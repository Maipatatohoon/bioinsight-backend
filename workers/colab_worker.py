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
            "--contrast", params.get("contrast", "")
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
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

        results_dict = {}
        warnings = []

        # Run DESeq2
        try:
            print_log("Running DESeq2...")
            df = run_r_script(deseq2_script, counts_path, meta_path, params)
            df = df.where(pd.notnull(df), None)
            results_dict["deseq2"] = df.to_dict(orient="records")
        except Exception as e:
            warnings.append(f"DESeq2: {str(e)}")

        # Run edgeR
        try:
            print_log("Running edgeR...")
            df = run_r_script(edger_script, counts_path, meta_path, params)
            df = df.where(pd.notnull(df), None)
            results_dict["edger"] = df.to_dict(orient="records")
        except Exception as e:
            warnings.append(f"edgeR: {str(e)}")

        # Run limma
        try:
            print_log("Running limma...")
            df = run_r_script(limma_script, counts_path, meta_path, params)
            df = df.where(pd.notnull(df), None)
            results_dict["limma"] = df.to_dict(orient="records")
        except Exception as e:
            warnings.append(f"limma: {str(e)}")

        if not results_dict:
            raise RuntimeError("All R engines failed.")

        # Check cancellation before upload
        print_log("Checking job status before upload...")
        try:
            r_status = http_session.get(f"{BACKEND_URL}/api/analysis/{job_id}", timeout=15)
            if r_status.status_code == 200:
                if r_status.json().get("state") == "CANCELLED":
                    print_log(f"Job {job_id} was cancelled. Discarding results.")
                    return
        except Exception as e:
            print_log(f"Status check failed: {e}. Proceeding with upload.")

        # Upload results
        print_log("Uploading results...")
        payload = {
            "worker_id": WORKER_ID,
            "results": results_dict,
            "warnings": warnings
        }
        r_submit = http_session.post(f"{BACKEND_URL}/api/workers/jobs/{job_id}/results", json=payload, headers=HEADERS, timeout=60)
        r_submit.raise_for_status()
        print_log(f"Job {job_id} completed successfully.")

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
