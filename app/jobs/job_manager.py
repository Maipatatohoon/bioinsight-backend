# backend/app/jobs/job_manager.py
"""
Asynchronous Job Manager for BioInsight 360 Backend.
Prevents browser timeouts on long-running RNA-seq analyses by tracking job state.
"""
import uuid
import time
import threading
import contextlib
import sqlite3
import pandas as pd
import numpy as np
import os
import re
import json
from typing import Dict, Any, Optional, List
from app.analysis.validator import validate_analysis_input, DataValidationError
from app.analysis.runner import REngineRunner
from app.analysis.qc import compute_qc_metrics

# Engines the backend knows how to execute. DESeq2 and edgeR are the two
# required methods for the dissertation workflow; limma-voom remains available
# as an optional third method.
ENGINES = ("deseq2", "edger", "limma")
CORE_ENGINES = ("deseq2", "edger")


def _initial_methods() -> Dict[str, Dict[str, Any]]:
    return {
        engine: {
            "status": "pending",
            "error": None,
            "rows": 0,
            "updated_at": None,
        }
        for engine in ENGINES
    }


def ensure_methods(job: Optional[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    """Return a complete method-status map for old and new job records.

    Older jobs (and results uploaded by an older worker) have no ``methods``
    key. Derive it from the saved result tables so the API and UI remain
    backwards compatible without rewriting historical rows.
    """
    if not job:
        return _initial_methods()
    methods = _initial_methods()
    stored = job.get("methods") or {}
    for engine in ENGINES:
        if isinstance(stored.get(engine), dict):
            methods[engine].update(stored[engine])
    for engine in ENGINES:
        if engine in (job.get("results") or {}) and methods[engine]["status"] == "pending":
            methods[engine]["status"] = "completed"
            methods[engine]["rows"] = len(job["results"].get(engine) or [])
    return methods


def _terminal_state(methods: Dict[str, Dict[str, Any]], results: Dict[str, Any]) -> str:
    completed = [
        engine for engine in CORE_ENGINES
        if methods.get(engine, {}).get("status") == "completed" and results.get(engine)
    ]
    if len(completed) == len(CORE_ENGINES):
        return JobState.COMPLETED
    if completed:
        return JobState.PARTIALLY_COMPLETED
    return JobState.FAILED


def _normalise_design_params(params: Dict[str, Any], meta_path: str) -> Dict[str, Any]:
    """Make design/contrast match the cleaned metadata column names.

    validate_analysis_input() renames whichever grouping column was resolved to
    ``condition``. Callers that still send the original name (e.g. the older
    frontend sent ``~ Group``) would otherwise fail with
    "object 'Group' not found". Rewrite unknown factor names to ``condition``.
    """
    params = dict(params or {})
    try:
        columns = list(pd.read_csv(meta_path, nrows=0).columns)
    except Exception:
        return params
    if "condition" not in columns:
        return params

    design = str(params.get("design") or "~ condition")
    rhs = design.split("~", 1)[1] if "~" in design else design
    rhs = re.sub(
        r"[A-Za-z_.][A-Za-z0-9_.]*",
        lambda m: m.group(0) if m.group(0) in columns else "condition",
        rhs,
    )
    params["design"] = "~ " + rhs.strip()
    params["condition_col"] = "condition"

    contrast = params.get("contrast")
    if isinstance(contrast, str) and contrast:
        factor, sep, rest = contrast.partition("_")
        if sep and factor and factor not in columns:
            params["contrast"] = "condition_" + rest
    return params


class JobState:
    QUEUED = "QUEUED"
    VALIDATING = "VALIDATING"
    RUNNING = "RUNNING"
    ANALYZING = "ANALYZING"
    ENRICHMENT = "ENRICHMENT"
    GENERATING_RESULTS = "GENERATING_RESULTS"
    COMPLETED = "COMPLETED"
    PARTIALLY_COMPLETED = "PARTIALLY_COMPLETED"
    FAILED = "FAILED"
    CANCELLED = "CANCELLED"
    WORKER_DISCONNECTED = "WORKER_DISCONNECTED"


class JobStore:
    _instance = None
    _DEFAULT_DB_FILENAME = "bioinsight_jobs.db"

    def __new__(cls):
        if cls._instance is None:
            inst = super(JobStore, cls).__new__(cls)
            inst.lock = threading.Lock()
            inst._db_path = os.environ.get("BIOINSIGHT_DB_PATH") or cls._default_db_path()
            parent = os.path.dirname(inst._db_path)
            if parent:
                os.makedirs(parent, exist_ok=True)
            with inst._connect() as conn:
                conn.execute("PRAGMA journal_mode=WAL")
                conn.execute("CREATE TABLE IF NOT EXISTS jobs (job_id TEXT PRIMARY KEY, data TEXT)")
                conn.execute("CREATE TABLE IF NOT EXISTS workers (worker_id TEXT PRIMARY KEY, version TEXT, status TEXT, last_heartbeat REAL, current_job_id TEXT)")
            cls._instance = inst
        return cls._instance

    @classmethod
    def _default_db_path(cls) -> str:
        # app/jobs/job_manager.py -> app/jobs -> app -> backend
        backend_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        return os.path.join(backend_root, "data", cls._DEFAULT_DB_FILENAME)

    @contextlib.contextmanager
    def _connect(self, immediate: bool = False):
        """Open a SQLite connection, guarantee it is closed, and use a write lock
        (BEGIN IMMEDIATE) for read-modify-write sequences when requested."""
        conn = sqlite3.connect(self._db_path, timeout=30, isolation_level=None)
        try:
            conn.execute("PRAGMA busy_timeout=30000")
            if immediate:
                conn.execute("BEGIN IMMEDIATE")
            yield conn
            if immediate:
                conn.execute("COMMIT")
        except Exception:
            if immediate:
                try:
                    conn.execute("ROLLBACK")
                except Exception:
                    pass
            raise
        finally:
            conn.close()

    def _get_job_raw(self, job_id):
        with self._connect() as conn:
            cursor = conn.execute("SELECT data FROM jobs WHERE job_id = ?", (job_id,))
            row = cursor.fetchone()
            if row:
                return json.loads(row[0])
        return None

    def _save_job_raw(self, job_id, job_data):
        with self._connect() as conn:
            conn.execute("INSERT OR REPLACE INTO jobs (job_id, data) VALUES (?, ?)", (job_id, json.dumps(job_data)))

    def create_job(self, counts_path: str, metadata_path: str, params: Dict[str, Any]) -> str:
        job_id = str(uuid.uuid4())
        now = time.time()
        job_data = {
            "job_id": job_id,
            "state": JobState.VALIDATING, # Start in validating, will move to QUEUED
            "created_at": now,
            "updated_at": now,
            "progress": 0,
            "message": "Job queued for validation",
            "warnings": [],
            "error": None,
            "counts_path": counts_path,
            "metadata_path": metadata_path,
            "params": params,
            "results": None,
            "methods": _initial_methods(),
            "comparison": None,
            "qc_data": None,
            "worker_id": None
        }
        with self.lock:
            self._save_job_raw(job_id, job_data)
        return job_id

    def get_job(self, job_id: str) -> Optional[Dict[str, Any]]:
        with self.lock:
            return self._get_job_raw(job_id)

    def update_job(self, job_id: str, state: str, progress: int, message: str, warnings: List[str] = None, error: str = None, results: Dict[str, Any] = None):
        with self.lock:
            with self._connect(immediate=True) as conn:
                cursor = conn.execute("SELECT data FROM jobs WHERE job_id = ?", (job_id,))
                row = cursor.fetchone()
                if not row:
                    return
                job_data = json.loads(row[0])
                job_data["state"] = state
                job_data["progress"] = progress
                job_data["message"] = message
                job_data["updated_at"] = time.time()
                if warnings is not None:
                    job_data["warnings"].extend(warnings)
                if error is not None:
                    job_data["error"] = error
                if results is not None:
                    job_data["results"] = results
                conn.execute("UPDATE jobs SET data = ? WHERE job_id = ?", (json.dumps(job_data), job_id))

    def update_job_results_atomic(self, job_id: str, results: Dict[str, Any], warnings: List[str] = None, worker_id: str = None) -> bool:
        """Atomically save results and derive the terminal state.

        A job is COMPLETED only when both core methods produced a table. If
        exactly one did, it is PARTIALLY_COMPLETED so the UI does not report a
        two-method analysis as fully successful. ``results`` is retained for
        backwards compatibility with workers/API clients that upload every
        engine in one request.
        """
        with self.lock:
            with self._connect(immediate=True) as conn:
                cursor = conn.execute("SELECT data FROM jobs WHERE job_id = ?", (job_id,))
                row = cursor.fetchone()
                if not row:
                    return False
                job_data = json.loads(row[0])
                if job_data.get("state") not in [JobState.RUNNING, JobState.GENERATING_RESULTS]:
                    return False
                if worker_id is not None and job_data.get("worker_id") not in (None, worker_id):
                    return False

                now = time.time()
                methods = ensure_methods(job_data)
                merged = job_data.get("results") or {}
                for engine, table in (results or {}).items():
                    if engine not in ENGINES:
                        continue
                    merged[engine] = table
                    methods[engine].update({
                        "status": "completed",
                        "error": None,
                        "rows": len(table or []),
                        "updated_at": now,
                    })
                # The aggregate endpoint is final. Any core method that was not
                # in this payload cannot become available later, so report it as
                # failed instead of leaving a misleading "pending" status.
                for engine in CORE_ENGINES:
                    if (engine not in (results or {})
                            and methods[engine]["status"] not in ("completed", "failed")):
                        methods[engine].update({
                            "status": "failed",
                            "error": "No result table was uploaded for this method.",
                            "updated_at": now,
                        })
                job_data["results"] = merged
                job_data["methods"] = methods
                if warnings:
                    job_data["warnings"].extend(warnings)

                state = _terminal_state(methods, merged)
                job_data["state"] = state
                job_data["updated_at"] = now
                if state == JobState.COMPLETED:
                    job_data["progress"] = 100
                    job_data["message"] = "DESeq2 and edgeR completed successfully."
                    job_data["error"] = None
                elif state == JobState.PARTIALLY_COMPLETED:
                    job_data["progress"] = 100
                    job_data["message"] = "One differential-expression method completed; the other failed."
                    job_data["error"] = "; ".join(
                        f"{e}: {methods[e].get('error')}" for e in CORE_ENGINES
                        if methods[e]["status"] == "failed"
                    ) or job_data.get("error")
                else:
                    job_data["progress"] = 0
                    job_data["message"] = "No differential-expression method completed."
                    job_data["error"] = job_data.get("error") or "; ".join(
                        f"{e}: {methods[e].get('error')}" for e in CORE_ENGINES
                        if methods[e].get("error")
                    ) or "All engines failed."

                cursor = conn.execute(
                    "UPDATE jobs SET data = ? WHERE job_id = ? AND json_extract(data, '$.state') IN ('RUNNING', 'GENERATING_RESULTS')",
                    (json.dumps(job_data), job_id),
                )
                return cursor.rowcount == 1

    def update_job_method(self, job_id: str, engine: str, status: str,
                          results: Optional[List[Dict[str, Any]]] = None,
                          error: Optional[str] = None,
                          worker_id: Optional[str] = None) -> Optional[Dict[str, Any]]:
        """Persist one method's running/completed/failed state.

        Called as soon as an engine finishes so the frontend can display it
        while the other engine is still running. The job remains RUNNING until
        ``finish_job`` (or the legacy all-at-once upload) decides the terminal
        state.
        """
        if engine not in ENGINES or status not in ("running", "completed", "failed"):
            return None
        with self.lock:
            with self._connect(immediate=True) as conn:
                row = conn.execute("SELECT data FROM jobs WHERE job_id = ?", (job_id,)).fetchone()
                if not row:
                    return None
                job_data = json.loads(row[0])
                if job_data.get("state") in (
                    JobState.COMPLETED, JobState.PARTIALLY_COMPLETED,
                    JobState.FAILED, JobState.CANCELLED,
                ):
                    return None
                if worker_id is not None and job_data.get("worker_id") not in (None, worker_id):
                    return None

                now = time.time()
                methods = ensure_methods(job_data)
                method = methods[engine]
                if status == "completed" and not results:
                    status = "failed"
                    error = "The engine returned an empty result table."
                method["status"] = status
                method["updated_at"] = now
                if status == "completed":
                    method["error"] = None
                    method["rows"] = len(results or [])
                    method.setdefault("results", None)
                    job_data.setdefault("results", {})
                    if job_data["results"] is None:
                        job_data["results"] = {}
                    job_data["results"][engine] = results or []
                elif status == "failed":
                    method["error"] = error or "Engine failed without an error message."

                job_data["methods"] = methods
                job_data["state"] = JobState.RUNNING
                job_data["updated_at"] = now
                done = [e for e in CORE_ENGINES if methods[e]["status"] == "completed"]
                job_data["message"] = (
                    f"{engine} {status}."
                    if status != "running"
                    else f"Running {engine}..."
                )
                job_data["progress"] = 40 + int(45 * len(done) / len(CORE_ENGINES))
                conn.execute("UPDATE jobs SET data = ? WHERE job_id = ?", (json.dumps(job_data), job_id))
                return job_data

    def finish_job(self, job_id: str, worker_id: Optional[str] = None,
                   warnings: Optional[List[str]] = None) -> Optional[str]:
        """Finalise a job after the worker has reported every method.

        Unlike the legacy all-at-once endpoint this does not need the result
        tables again; they are already in the job record.
        """
        with self.lock:
            with self._connect(immediate=True) as conn:
                row = conn.execute("SELECT data FROM jobs WHERE job_id = ?", (job_id,)).fetchone()
                if not row:
                    return None
                job_data = json.loads(row[0])
                if job_data.get("state") in (
                    JobState.COMPLETED, JobState.PARTIALLY_COMPLETED,
                    JobState.FAILED, JobState.CANCELLED,
                ):
                    return job_data["state"]
                if worker_id is not None and job_data.get("worker_id") not in (None, worker_id):
                    return None

                now = time.time()
                methods = ensure_methods(job_data)
                # Core methods must report success or failure. limma-voom is
                # optional: if the worker never started it, leave it pending so
                # the UI says "not run" instead of inventing a failure.
                for engine in CORE_ENGINES:
                    if methods[engine]["status"] in ("pending", "running"):
                        methods[engine]["status"] = "failed"
                        methods[engine]["error"] = "No result was reported by the worker."
                        methods[engine]["updated_at"] = now
                if methods["limma"]["status"] == "running":
                    methods["limma"]["status"] = "failed"
                    methods["limma"]["error"] = "No result was reported by the worker."
                    methods["limma"]["updated_at"] = now
                results = job_data.get("results") or {}
                state = _terminal_state(methods, results)

                job_data["methods"] = methods
                job_data["state"] = state
                job_data["updated_at"] = now
                if warnings:
                    job_data["warnings"].extend(warnings)
                if state == JobState.COMPLETED:
                    job_data["progress"] = 100
                    job_data["message"] = "DESeq2 and edgeR completed successfully."
                elif state == JobState.PARTIALLY_COMPLETED:
                    job_data["progress"] = 100
                    job_data["message"] = "One differential-expression method completed; the other failed."
                    job_data["error"] = "; ".join(
                        f"{e}: {methods[e].get('error')}" for e in CORE_ENGINES
                        if methods[e]["status"] == "failed"
                    )
                else:
                    job_data["progress"] = 0
                    job_data["message"] = "No differential-expression method completed."
                    job_data["error"] = "; ".join(
                        f"{e}: {methods[e].get('error')}" for e in CORE_ENGINES
                    ) or "All engines failed."
                conn.execute("UPDATE jobs SET data = ? WHERE job_id = ?", (json.dumps(job_data), job_id))
                return state

    def save_comparison(self, job_id: str, comparison: Dict[str, Any]) -> bool:
        """Persist a user-triggered DESeq2-vs-edgeR comparison on the job."""
        with self.lock:
            with self._connect(immediate=True) as conn:
                row = conn.execute("SELECT data FROM jobs WHERE job_id = ?", (job_id,)).fetchone()
                if not row:
                    return False
                job_data = json.loads(row[0])
                job_data["comparison"] = comparison
                job_data["updated_at"] = time.time()
                conn.execute("UPDATE jobs SET data = ? WHERE job_id = ?", (json.dumps(job_data), job_id))
                return True

    def update_job_qc(self, job_id: str, qc_data: dict):
        with self.lock:
            with self._connect(immediate=True) as conn:
                cursor = conn.execute("SELECT data FROM jobs WHERE job_id = ?", (job_id,))
                row = cursor.fetchone()
                if not row:
                    return
                job_data = json.loads(row[0])
                job_data["qc_data"] = qc_data
                job_data["updated_at"] = time.time()
                conn.execute("UPDATE jobs SET data = ? WHERE job_id = ?", (json.dumps(job_data), job_id))

    def cancel_job(self, job_id: str) -> bool:
        with self.lock:
            with self._connect(immediate=True) as conn:
                cursor = conn.execute("SELECT data FROM jobs WHERE job_id = ?", (job_id,))
                row = cursor.fetchone()
                if not row:
                    return False
                job_data = json.loads(row[0])
                if job_data["state"] in [JobState.COMPLETED, JobState.FAILED, JobState.CANCELLED, JobState.WORKER_DISCONNECTED]:
                    return False
                job_data["state"] = JobState.CANCELLED
                job_data["message"] = "Job cancelled by user."
                job_data["updated_at"] = time.time()
                cursor = conn.execute(
                    "UPDATE jobs SET data = ? WHERE job_id = ? AND json_extract(data, '$.state') NOT IN ('COMPLETED','FAILED','CANCELLED','WORKER_DISCONNECTED')",
                    (json.dumps(job_data), job_id),
                )
                return cursor.rowcount == 1

    # --- Worker Management Methods ---
    def register_worker(self, worker_id: str, version: str):
        self.reap_stale_workers()
        now = time.time()
        with self.lock:
            with self._connect() as conn:
                conn.execute(
                    "INSERT OR REPLACE INTO workers (worker_id, version, status, last_heartbeat, current_job_id) VALUES (?, ?, 'ONLINE', ?, NULL)",
                    (worker_id, version, now)
                )
    
    def heartbeat_worker(self, worker_id: str, current_job_id: Optional[str] = None):
        now = time.time()
        with self.lock:
            with self._connect() as conn:
                conn.execute(
                    "UPDATE workers SET last_heartbeat = ?, current_job_id = ?, status = 'ONLINE' WHERE worker_id = ?",
                    (now, current_job_id, worker_id)
                )

    def reap_stale_workers(self, timeout_seconds=60):
        now = time.time()
        with self.lock:
            with self._connect(immediate=True) as conn:
                cursor = conn.execute(
                    "SELECT worker_id, current_job_id FROM workers WHERE status = 'ONLINE' AND (? - last_heartbeat) > ?",
                    (now, timeout_seconds)
                )
                stale_workers = cursor.fetchall()
                
                for worker_id, current_job_id in stale_workers:
                    conn.execute("UPDATE workers SET status = 'OFFLINE', current_job_id = NULL WHERE worker_id = ?", (worker_id,))
                    if current_job_id:
                        cursor2 = conn.execute("SELECT data FROM jobs WHERE job_id = ?", (current_job_id,))
                        row = cursor2.fetchone()
                        if row:
                            job_data = json.loads(row[0])
                            if job_data["state"] == JobState.RUNNING:
                                job_data["state"] = JobState.WORKER_DISCONNECTED
                                job_data["message"] = "Worker disconnected. Job failed."
                                job_data["error"] = "Worker heartbeat timeout"
                                job_data["updated_at"] = now
                                conn.execute("UPDATE jobs SET data = ? WHERE job_id = ?", (json.dumps(job_data), current_job_id))

    def get_next_available_job(self) -> Optional[Dict[str, Any]]:
        """Finds a job in QUEUED state."""
        self.reap_stale_workers()
        with self.lock:
            with self._connect() as conn:
                # We use JSON1 to extract state
                cursor = conn.execute("SELECT job_id, data FROM jobs WHERE json_extract(data, '$.state') = 'QUEUED' ORDER BY json_extract(data, '$.created_at') ASC LIMIT 1")
                row = cursor.fetchone()
                if row:
                    return json.loads(row[1])
        return None

    def claim_job(self, job_id: str, worker_id: str) -> bool:
        """Atomically claim a job. Exactly one caller wins even across processes,
        because the state check and the write happen in a single IMMEDIATE transaction."""
        with self.lock:
            with self._connect(immediate=True) as conn:
                cursor = conn.execute(
                    "SELECT data FROM jobs WHERE job_id = ? AND json_extract(data, '$.state') = 'QUEUED'",
                    (job_id,),
                )
                row = cursor.fetchone()
                if not row:
                    return False
                job_data = json.loads(row[0])
                job_data["state"] = JobState.RUNNING
                job_data["worker_id"] = worker_id
                job_data["updated_at"] = time.time()
                job_data["message"] = f"Job claimed by worker {worker_id}"
                cursor = conn.execute(
                    "UPDATE jobs SET data = ? WHERE job_id = ? AND json_extract(data, '$.state') = 'QUEUED'",
                    (json.dumps(job_data), job_id),
                )
                return cursor.rowcount == 1

class ComputeRuntime:
    """Base abstraction for compute providers."""
    def dispatch(self, store: JobStore, job_id: str):
        raise NotImplementedError()

class ColabRuntime(ComputeRuntime):
    """Dispatches job to the remote Colab worker by placing it in QUEUED state."""
    def dispatch(self, store: JobStore, job_id: str):
        store.update_job(job_id, JobState.QUEUED, 30, "Validation complete. Waiting for Colab compute worker to claim job.")

class LocalRuntime(ComputeRuntime):
    """Runs the job on this machine by invoking the R engines via Rscript."""
    def dispatch(self, store: JobStore, job_id: str):
        job = store.get_job(job_id)
        if not job:
            return
        counts_path = job["counts_path"] + ".clean.csv"
        meta_path = job["metadata_path"] + ".clean.csv"
        app_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        engines = [
            ("deseq2", os.path.join(app_dir, "deseq2", "deseq2_script.R")),
            ("edger", os.path.join(app_dir, "edger", "edger_script.R")),
            ("limma", os.path.join(app_dir, "limma", "limma_script.R")),
        ]
        params = _normalise_design_params(job.get("params", {}) or {}, meta_path)
        warnings = []
        for i, (name, script) in enumerate(engines):
            progress = 40 + int(40 * i / len(engines))
            store.update_job_method(job_id, name, "running", worker_id="local")
            store.update_job(job_id, JobState.RUNNING, progress, f"Running {name} locally...")
            try:
                df = REngineRunner.run_r_script(script, counts_path, meta_path, params)
                df = df.replace([float("inf"), float("-inf")], float("nan"))
                df = df.where(pd.notnull(df), None)
                table = df.to_dict(orient="records")
                store.update_job_method(job_id, name, "completed", results=table, worker_id="local")
            except Exception as e:  # one engine failing must not sink the others
                message = f"{name}: {e}"
                warnings.append(message)
                store.update_job_method(job_id, name, "failed", error=str(e), worker_id="local")

        store.update_job(job_id, JobState.RUNNING, 95, "Saving results...")
        state = store.finish_job(job_id, worker_id="local", warnings=warnings)
        if state is None:
            # The job disappeared or belongs to another worker. Never fabricate a
            # completed state in that case.
            return

class CloudRuntime(ComputeRuntime):
    """(Future) Runs the job on AWS/GCP Batch."""
    def dispatch(self, store: JobStore, job_id: str):
        pass

def process_job_async(job_id: str):
    """
    Locally performs validation and QC using Pandas. 
    If successful, places the job into QUEUED state so a remote worker can claim it.
    """
    store = JobStore()
    job = store.get_job(job_id)
    if not job:
        return

    try:
        store.update_job(job_id, JobState.VALIDATING, 10, "Validating input count matrix and metadata...")
        
        counts_sep = '\t' if job["counts_path"].endswith(('.tsv', '.txt')) else ','
        meta_sep = '\t' if job["metadata_path"].endswith(('.tsv', '.txt')) else ','
        
        counts_df = pd.read_csv(job["counts_path"], sep=counts_sep)
        meta_df = pd.read_csv(job["metadata_path"], sep=meta_sep)
        params = job["params"]

        clean_counts, clean_meta, warnings = validate_analysis_input(counts_df, meta_df, params)

        store.update_job(job_id, JobState.VALIDATING, 20, "Computing Quality Control metrics (PCA, Correlations)...", warnings=warnings)
        qc_data = compute_qc_metrics(clean_counts, clean_meta)
        store.update_job_qc(job_id, qc_data)
        
        clean_counts_path = job["counts_path"] + ".clean.csv"
        clean_meta_path = job["metadata_path"] + ".clean.csv"
        clean_counts.to_csv(clean_counts_path)
        clean_meta.to_csv(clean_meta_path)

        # Dispatch via the configured runtime
        runtime_name = os.environ.get("BIOINSIGHT_RUNTIME", "colab").strip().lower()
        runtime = LocalRuntime() if runtime_name == "local" else ColabRuntime()
        runtime.dispatch(store, job_id)

    except DataValidationError as e:
        store.update_job(job_id, JobState.FAILED, 0, f"Scientific Validation Error: {e.message}", error=e.message)
    except Exception as e:
        store.update_job(job_id, JobState.FAILED, 0, f"Pipeline Error: {str(e)}", error=str(e))
