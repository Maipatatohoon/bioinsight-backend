import os
import re
import sys
import time
import uuid
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import tempfile
import threading
import subprocess
import pandas as pd

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
BACKEND_URL = os.environ.get("BIOINSIGHT_BACKEND_URL", "http://localhost:8000")
WORKER_TOKEN = os.environ.get("BIOINSIGHT_WORKER_TOKEN", "dev-token-123")
WORKER_ID = f"worker-{uuid.uuid4().hex[:8]}"

HEADERS = {
    "Authorization": f"Bearer {WORKER_TOKEN}"
}

def print_log(msg: str):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] [{WORKER_ID}] {msg}", flush=True)


def _normalise_design(params: dict, meta_path: str) -> dict:
    """Match design/contrast factor names to the cleaned metadata columns.

    The API stores the grouping column as `condition` in the .clean.csv. Older
    callers still send the original name (e.g. `~ Group`); rewrite it so the R
    engines do not fail with "object 'Group' not found".
    """
    try:
        columns = list(pd.read_csv(meta_path, nrows=0).columns)
    except Exception:
        return params
    target = "condition" if "condition" in columns else next(
        (c for c in columns if c.lower() not in ("sample", "sample_id", "id")), None
    )
    if not target:
        return params

    params = dict(params or {})
    design = str(params.get("design") or f"~ {target}")
    rhs = design.split("~", 1)[1] if "~" in design else design
    rhs = re.sub(
        r"[A-Za-z_.][A-Za-z0-9_.]*",
        lambda match: match.group(0) if match.group(0) in columns else target,
        rhs,
    )
    params["design"] = "~ " + rhs.strip()

    contrast = params.get("contrast")
    if isinstance(contrast, str) and contrast:
        factor, sep, rest = contrast.partition("_")
        if sep and factor and factor not in columns:
            params["contrast"] = f"{target}_" + rest
    return params

def get_session():
    s = requests.Session()
    retry = Retry(
        total=5,
        backoff_factor=1,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["HEAD", "GET", "OPTIONS", "POST"]
    )
    adapter = HTTPAdapter(max_retries=retry)
    s.mount("http://", adapter)
    s.mount("https://", adapter)
    return s

http_session = get_session()

# ---------------------------------------------------------------------------
# Heartbeat Thread
# ---------------------------------------------------------------------------
class HeartbeatThread(threading.Thread):
    def __init__(self):
        super().__init__(daemon=True)
        self.current_job_id = None
        
    def run(self):
        while True:
            try:
                payload = {"worker_id": WORKER_ID, "current_job_id": self.current_job_id}
                http_session.post(f"{BACKEND_URL}/api/workers/heartbeat", json=payload, headers=HEADERS, timeout=10)
            except Exception as e:
                print_log(f"Heartbeat failed: {e}")
            time.sleep(15)

# ---------------------------------------------------------------------------
# Execution Logic
# ---------------------------------------------------------------------------
def run_r_script(script_path: str, counts_path: str, meta_path: str, params: dict) -> pd.DataFrame:
    with tempfile.TemporaryDirectory() as tmpdir:
        out_path = os.path.join(tmpdir, "results.csv")
        cmd = [
            "Rscript", script_path,
            "--counts", counts_path,
            "--meta", meta_path,
            "--out", out_path,
            "--design", params.get("design", "~ condition"),
            "--contrast", params.get("contrast", ""),
            "--numerator", str(params.get("treated_group", "")),
            "--denominator", str(params.get("control_group", "")),
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=900)
            if os.path.exists(out_path):
                df = pd.read_csv(out_path)
                missing = [c for c in ("Gene", "log2FoldChange", "pvalue", "padj") if c not in df.columns]
                if missing:
                    raise RuntimeError("R output is missing required column(s): " + ", ".join(missing))
                return df
            else:
                raise RuntimeError("R script completed but output file is missing.")
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"R script failed: {e.stderr}")

def _records(df: pd.DataFrame):
    """Convert a result frame to JSON-safe records (NaN/Inf -> null)."""
    clean = df.replace([float("inf"), float("-inf")], float("nan"))
    clean = clean.astype(object).where(pd.notnull(clean), None)
    return clean.to_dict(orient="records")


def _post_method(job_id: str, engine: str, status: str, results=None, error=None) -> bool:
    """Report one method's state. Returns False when the backend is too old."""
    payload = {"worker_id": WORKER_ID, "status": status}
    if results is not None:
        payload["results"] = results
    if error is not None:
        payload["error"] = error
    try:
        response = http_session.post(
            f"{BACKEND_URL}/api/workers/jobs/{job_id}/method/{engine}",
            json=payload,
            headers=HEADERS,
            timeout=120,
        )
        if response.status_code in (404, 405):
            return False
        response.raise_for_status()
        return True
    except requests.exceptions.RequestException as exc:
        print_log(f"Could not report {engine} {status}: {exc}")
        return False


def _post_finish(job_id: str, warnings: list) -> bool:
    try:
        response = http_session.post(
            f"{BACKEND_URL}/api/workers/jobs/{job_id}/finish",
            json={"worker_id": WORKER_ID, "warnings": warnings},
            headers=HEADERS,
            timeout=30,
        )
        if response.status_code in (404, 405):
            return False
        response.raise_for_status()
        return True
    except requests.exceptions.RequestException as exc:
        print_log(f"Could not finalize job: {exc}")
        return False


def _upload_legacy(job_id: str, results: dict, warnings: list):
    """Fallback for an older backend without per-method endpoints."""
    payload = {"worker_id": WORKER_ID, "results": results, "warnings": warnings}
    response = http_session.post(
        f"{BACKEND_URL}/api/workers/jobs/{job_id}/results",
        json=payload,
        headers=HEADERS,
        timeout=120,
    )
    response.raise_for_status()


def _is_cancelled(job_id: str) -> bool:
    try:
        response = http_session.get(f"{BACKEND_URL}/api/analysis/{job_id}", timeout=15)
        if response.status_code == 200:
            return response.json().get("state") == "CANCELLED"
    except Exception as exc:
        print_log(f"Status check failed: {exc}. Continuing.")
    return False


def process_job(job_id: str, params: dict):
    print_log(f"Processing job {job_id}")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        counts_path = os.path.join(tmpdir, "counts.csv")
        meta_path = os.path.join(tmpdir, "meta.csv")
        
        # Download datasets
        print_log("Downloading datasets...")
        r_counts = http_session.get(f"{BACKEND_URL}/api/workers/jobs/{job_id}/dataset/counts", headers=HEADERS, timeout=60)
        r_counts.raise_for_status()
        with open(counts_path, "wb") as f:
            f.write(r_counts.content)
            
        r_meta = http_session.get(f"{BACKEND_URL}/api/workers/jobs/{job_id}/dataset/metadata", headers=HEADERS, timeout=60)
        r_meta.raise_for_status()
        with open(meta_path, "wb") as f:
            f.write(r_meta.content)

        params = _normalise_design(params, meta_path)

        # Locate R scripts
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        deseq2_script = os.path.join(base_dir, "app", "deseq2", "deseq2_script.R")
        edger_script = os.path.join(base_dir, "app", "edger", "edger_script.R")
        limma_script = os.path.join(base_dir, "app", "limma", "limma_script.R")

        engines = [
            ("deseq2", deseq2_script),
            ("edger", edger_script),
            ("limma", limma_script),
        ]
        results_dict = {}
        warnings = []
        legacy = False

        for engine, script in engines:
            if _is_cancelled(job_id):
                print_log(f"Job {job_id} was cancelled. Discarding remaining work.")
                return
            if not legacy:
                if not _post_method(job_id, engine, "running"):
                    print_log("Per-method API unavailable; falling back to the legacy all-at-once upload.")
                    legacy = True
            try:
                print_log(f"Running {engine}...")
                df = run_r_script(script, counts_path, meta_path, params)
                table = _records(df)
                results_dict[engine] = table
                print_log(f"{engine} completed: {len(table)} genes.")
                if not legacy and not _post_method(job_id, engine, "completed", results=table):
                    # The backend accepted the per-method protocol but this
                    # table did not persist. Fall back to the legacy aggregate
                    # upload for this and every remaining engine so a network
                    # hiccup cannot silently lose a completed result.
                    print_log("Could not save this method incrementally; falling back to aggregate upload.")
                    legacy = True
            except Exception as e:
                warnings.append(f"{engine}: {str(e)}")
                print_log(f"{engine} failed: {e}")
                if not legacy:
                    _post_method(job_id, engine, "failed", error=str(e))

        if _is_cancelled(job_id):
            print_log(f"Job {job_id} was cancelled. Discarding results.")
            return

        if legacy:
            print_log("Uploading aggregate results to legacy endpoint...")
            _upload_legacy(job_id, results_dict, warnings)
            print_log(f"Job {job_id} finalized through the legacy endpoint.")
            return

        print_log("Finalizing job...")
        if not _post_finish(job_id, warnings):
            print_log("Finalize endpoint unavailable; sending aggregate results instead.")
            _upload_legacy(job_id, results_dict, warnings)
        print_log(f"Job {job_id} finalized.")

# ---------------------------------------------------------------------------
# Main Polling Loop
# ---------------------------------------------------------------------------
def main():
    print_log(f"Starting BioInsight 360 Headless Worker")
    print_log(f"Connecting to Backend: {BACKEND_URL}")
    
    # Register
    try:
        r = http_session.post(f"{BACKEND_URL}/api/workers/register", json={"worker_id": WORKER_ID, "version": "2.0.0"}, headers=HEADERS, timeout=15)
        r.raise_for_status()
        print_log("Successfully registered with backend.")
    except Exception as e:
        print_log(f"Failed to register with backend: {e}")
        sys.exit(1)
        
    heartbeat = HeartbeatThread()
    heartbeat.start()
    
    while True:
        try:
            # Poll for next job
            r = http_session.get(f"{BACKEND_URL}/api/workers/jobs/next", headers=HEADERS, timeout=15)
            r.raise_for_status()
            data = r.json()
            job_id = data.get("job_id")
            
            if job_id:
                # Attempt to claim
                claim_r = http_session.post(
                    f"{BACKEND_URL}/api/workers/jobs/{job_id}/claim", 
                    data={"worker_id": WORKER_ID}, 
                    headers=HEADERS,
                    timeout=15
                )
                if claim_r.status_code == 200:
                    heartbeat.current_job_id = job_id
                    try:
                        process_job(job_id, data.get("params", {}))
                    except Exception as e:
                        print_log(f"Job {job_id} failed: {e}")
                        http_session.post(
                            f"{BACKEND_URL}/api/workers/jobs/{job_id}/failed", 
                            data={"error": str(e)}, 
                            headers=HEADERS,
                            timeout=15
                        )
                    finally:
                        heartbeat.current_job_id = None
                else:
                    print_log(f"Failed to claim job {job_id} (might be claimed by another worker).")
            else:
                # No jobs available
                time.sleep(5)
                
        except requests.exceptions.RequestException as e:
            print_log(f"Connection error while polling: {e}")
            time.sleep(10)
        except Exception as e:
            print_log(f"Unexpected error: {e}")
            time.sleep(10)

if __name__ == "__main__":
    main()
