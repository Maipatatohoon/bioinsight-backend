# backend/app/main.py
from fastapi import FastAPI, UploadFile, File, HTTPException, BackgroundTasks, Form, Depends, Security
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import Response, StreamingResponse, FileResponse
from fastapi.security import HTTPBearer, HTTPAuthorizationCredentials
from pydantic import BaseModel
from typing import Dict, Any, Optional, List
import pandas as pd
import io
import os
import hmac
import math
import tempfile

from app.analysis.profiler import profile_dataset
from app.analysis.models import ProfilingReport
from app.analysis.runner import REngineRunner
from app.analysis.enrichment import enrich_genes
from app.analysis.comparison import compute_comparison
from app.jobs.job_manager import (
    JobStore, JobState, process_job_async, ensure_methods, ENGINES, CORE_ENGINES,
)
from app.storage import save_upload
import json

app = FastAPI(
    title="BioInsight 360 Computational Engine",
    version="2.0.0",
    description="Scientific microservice API for RNA-seq differential expression (DESeq2, edgeR, limma-voom)."
)


def _json_safe(value):
    """Replace NaN/Infinity with null so Starlette can render strict JSON.

    DESeq2/edgeR/limma legitimately return NA for all-zero or filtered genes.
    pandas turns those into NaN (and occasionally +/-Inf for log2FC), which are
    valid Python floats but invalid JSON. FastAPI/Starlette then raises
    "Out of range float values are not JSON compliant" and the whole results
    request fails with a 500. Emitting null instead lets the UI skip the row.
    """
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value

DEFAULT_WORKER_TOKEN = "dev-token-123"
BIOINSIGHT_ENV = os.environ.get("BIOINSIGHT_ENV", "development")
WORKER_TOKEN = os.environ.get("BIOINSIGHT_WORKER_TOKEN", DEFAULT_WORKER_TOKEN)

if BIOINSIGHT_ENV == "production" and WORKER_TOKEN == DEFAULT_WORKER_TOKEN:
    raise RuntimeError(
        "BIOINSIGHT_WORKER_TOKEN must be set to a non-default value when "
        "BIOINSIGHT_ENV=production."
    )


def _allowed_origins():
    """Browser origins allowed by CORS. Configurable via BIOINSIGHT_ALLOWED_ORIGINS
    (comma-separated, or '*' to allow any). Defaults to local development."""
    raw = os.environ.get("BIOINSIGHT_ALLOWED_ORIGINS", "").strip()
    if raw == "*":
        return ["*"]
    origins = [o.strip() for o in raw.split(",") if o.strip()]
    return origins or ["http://localhost:3000", "http://localhost:8000"]


_CORS_ORIGINS = _allowed_origins()

app.add_middleware(
    CORSMiddleware,
    allow_origins=_CORS_ORIGINS,
    # '*' cannot be combined with credentials; we don't send credentials from the UI.
    allow_credentials=("*" not in _CORS_ORIGINS),
    allow_methods=["*"],
    allow_headers=["*"],
)

class AnalysisRequest(BaseModel):
    counts: Dict[str, Dict[str, Any]]
    metadata: Dict[str, Dict[str, Any]]
    params: Dict[str, Any]

@app.get("/")
def read_root():
    return {
        "title": "BioInsight 360 Backend",
        "status": "online",
        "supported_engines": ["DESeq2", "edgeR", "limma-voom"]
    }

@app.get("/health")
def health_check():
    return {"status": "healthy", "service": "bioinsight-backend"}

@app.post("/api/profile", response_model=ProfilingReport)
async def profile_upload(file: UploadFile = File(...)):
    if not file.filename.endswith(('.csv', '.tsv', '.txt')):
        raise HTTPException(status_code=400, detail="Only CSV or TSV files are supported.")
        
    try:
        contents = await file.read()
        sep = '\t' if file.filename.endswith(('.tsv', '.txt')) else ','
        df = pd.read_csv(io.StringIO(contents.decode('utf-8')), sep=sep)
        report = profile_dataset(df)
        return report
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error profiling file: {str(e)}")

# ==============================================================================
# Job-Oriented Workflow Endpoints
# ==============================================================================

@app.post("/api/analysis", status_code=202)
def submit_analysis_job(
    background_tasks: BackgroundTasks,
    counts_file: UploadFile = File(...),
    metadata_file: UploadFile = File(...),
    params: str = Form("{}")
):
    parsed_params = {}
    if params:
        try:
            parsed_params = json.loads(params)
        except json.JSONDecodeError:
            raise HTTPException(status_code=400, detail="Invalid JSON in params")
    if not isinstance(parsed_params, dict):
        raise HTTPException(status_code=400, detail="params must be a JSON object.")

    # Direction is part of the scientific result. Do not let an ambiguous
    # request fall back to alphabetical factor levels inside DESeq2/edgeR.
    control = str(parsed_params.get("control_group") or "").strip()
    treated = str(parsed_params.get("treated_group") or "").strip()
    if not control or not treated:
        raise HTTPException(
            status_code=400,
            detail="control_group and treated_group are required for raw-count analysis.",
        )
    if control == treated:
        raise HTTPException(status_code=400, detail="control_group and treated_group must differ.")
    parsed_params["control_group"] = control
    parsed_params["treated_group"] = treated

    # Save only after the request has passed the direction checks, so a rejected
    # submission cannot leave orphaned upload files.
    counts_path = save_upload(counts_file, "counts")
    metadata_path = save_upload(metadata_file, "metadata")

    store = JobStore()
    job_id = store.create_job(counts_path, metadata_path, parsed_params)
    background_tasks.add_task(process_job_async, job_id)
    created = store.get_job(job_id)
    return {
        "job_id": job_id,
        "state": created["state"] if created else JobState.VALIDATING,
        "message": "Analysis job submitted successfully. Poll GET /api/analysis/{job_id} for progress."
    }

@app.get("/api/analysis/{job_id}")
def get_job_status(job_id: str):
    store = JobStore()
    job = store.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job ID '{job_id}' not found.")

    results = job.get("results") or {}
    return _json_safe({
        "job_id": job["job_id"],
        "state": job["state"],
        "progress": job["progress"],
        "message": job["message"],
        "warnings": job["warnings"],
        "error": job["error"],
        "created_at": job["created_at"],
        "updated_at": job["updated_at"],
        "params": job.get("params") or {},
        "methods": ensure_methods(job),
        "results_available": sorted(results.keys()),
        "comparison_available": bool(job.get("comparison")),
        "qc_data": job.get("qc_data")
    })

@app.get("/api/analysis/{job_id}/results")
def get_job_results(job_id: str):
    """Return whatever method tables are saved, even while another is running.

    This is the key contract for individual-method results: DESeq2 must be
    viewable as soon as it finishes, without waiting for edgeR or a consensus.
    The response is marked ``partial`` until both core methods are complete.
    """
    store = JobStore()
    job = store.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job ID '{job_id}' not found.")

    results = job.get("results") or {}
    if not results:
        return _json_safe({
            "job_id": job_id,
            "status": "pending",
            "state": job["state"],
            "methods": ensure_methods(job),
            "warnings": job["warnings"],
            "error": job["error"],
            "results": {},
        })

    return _json_safe({
        "job_id": job_id,
        "status": "success" if job["state"] == JobState.COMPLETED else "partial",
        "state": job["state"],
        "methods": ensure_methods(job),
        "warnings": job["warnings"],
        "error": job["error"],
        "results": results,
    })


@app.get("/api/analysis/{job_id}/method/{engine}")
def get_method_results(job_id: str, engine: str):
    """Return one engine's saved table without requiring the other engine."""
    engine = engine.lower()
    if engine not in ENGINES:
        raise HTTPException(status_code=400, detail=f"Unknown engine '{engine}'. Choose one of: {', '.join(ENGINES)}.")
    store = JobStore()
    job = store.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job ID '{job_id}' not found.")
    results = job.get("results") or {}
    if engine not in results:
        method = ensure_methods(job).get(engine, {})
        raise HTTPException(
            status_code=404,
            detail=f"{engine} results are not available yet (status: {method.get('status', 'pending')}).",
        )
    return _json_safe({
        "job_id": job_id,
        "engine": engine,
        "state": job["state"],
        "methods": ensure_methods(job),
        "warnings": job["warnings"],
        "results": results[engine],
    })

@app.post("/api/analysis/{job_id}/cancel")
def cancel_job(job_id: str):
    store = JobStore()
    success = store.cancel_job(job_id)
    if not success:
        raise HTTPException(status_code=400, detail="Job cannot be cancelled (not found or already finished).")
    return {"status": "cancelled", "job_id": job_id}

@app.get("/api/analysis/{job_id}/files")
def download_job_results_csv(job_id: str, engine: Optional[str] = "deseq2"):
    store = JobStore()
    job = store.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job '{job_id}' not found.")

    engine_name = (engine or "deseq2").lower()
    if engine_name not in ENGINES:
        raise HTTPException(status_code=400, detail=f"Unknown engine '{engine_name}'.")
    results = job.get("results") or {}
    engine_data = results.get(engine_name)
    if not engine_data:
        raise HTTPException(status_code=404, detail=f"Engine '{engine_name}' results not available for job '{job_id}'.")

    df = pd.DataFrame(engine_data)
    stream = io.StringIO()
    df.to_csv(stream, index=False)

    return Response(
        content=stream.getvalue(),
        media_type="text/csv",
        headers={"Content-Disposition": f"attachment; filename=bioinsight_{engine_name}_results_{job_id[:8]}.csv"}
    )


# ==============================================================================
# Optional comparison (computed from the two saved result tables)
# ==============================================================================
class ComparisonRequest(BaseModel):
    padj_threshold: float = 0.05
    log2fc_threshold: float = 1.0


@app.post("/api/analysis/{job_id}/comparison")
def create_comparison(job_id: str, request: ComparisonRequest = ComparisonRequest()):
    store = JobStore()
    job = store.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job ID '{job_id}' not found.")
    results = job.get("results") or {}
    missing = [engine for engine in CORE_ENGINES if not results.get(engine)]
    if missing:
        raise HTTPException(
            status_code=409,
            detail=(
                "The comparison needs saved DESeq2 and edgeR results. Missing: "
                + ", ".join(missing)
            ),
        )
    try:
        comparison = compute_comparison(
            results["deseq2"],
            results["edger"],
            padj_threshold=request.padj_threshold,
            log2fc_threshold=request.log2fc_threshold,
        )
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    comparison["job_id"] = job_id
    if not store.save_comparison(job_id, comparison):
        raise HTTPException(status_code=500, detail="Could not save the comparison.")
    return _json_safe({"status": "success", "comparison": comparison})


@app.get("/api/analysis/{job_id}/comparison")
def get_comparison(job_id: str):
    store = JobStore()
    job = store.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job ID '{job_id}' not found.")
    comparison = job.get("comparison")
    if not comparison:
        raise HTTPException(status_code=404, detail="No comparison has been computed for this job yet.")
    return _json_safe({"status": "success", "comparison": comparison})


@app.get("/api/analysis/{job_id}/comparison.csv")
def download_comparison_csv(job_id: str):
    store = JobStore()
    job = store.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job ID '{job_id}' not found.")
    comparison = job.get("comparison")
    if not comparison:
        raise HTTPException(status_code=404, detail="No comparison has been computed for this job yet.")
    rows = comparison.get("rows") or []
    df = pd.DataFrame(rows)
    stream = io.StringIO()
    df.to_csv(stream, index=False)
    return Response(
        content=stream.getvalue(),
        media_type="text/csv",
        headers={"Content-Disposition": f"attachment; filename=bioinsight_comparison_{job_id[:8]}.csv"}
    )

# Synchronous endpoints for legacy backwards compatibility
def _run_sync_engine(script_parts: tuple, request: AnalysisRequest):
    """
    Runs a single R engine synchronously for a JSON request. The runner consumes
    file paths, so the in-memory counts/metadata are written to a temp dir first.
    """
    counts_df = pd.DataFrame(request.counts)
    meta_df = pd.DataFrame(request.metadata)
    script_path = os.path.join(os.path.dirname(__file__), *script_parts)
    with tempfile.TemporaryDirectory() as tmpdir:
        counts_path = os.path.join(tmpdir, "counts.csv")
        meta_path = os.path.join(tmpdir, "metadata.csv")
        counts_df.to_csv(counts_path)
        meta_df.to_csv(meta_path)
        result_df = REngineRunner.run_r_script(script_path, counts_path, meta_path, request.params)
    result_df = result_df.where(pd.notnull(result_df), None)
    return {"status": "success", "results": result_df.to_dict(orient="records")}

@app.post("/api/analyze/deseq2")
def run_deseq2_sync(request: AnalysisRequest):
    try:
        return _run_sync_engine(("deseq2", "deseq2_script.R"), request)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"DESeq2 analysis failed: {str(e)}")

@app.post("/api/analyze/edger")
def run_edger_sync(request: AnalysisRequest):
    try:
        return _run_sync_engine(("edger", "edger_script.R"), request)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"edgeR analysis failed: {str(e)}")

@app.post("/api/analyze/limma")
def run_limma_sync(request: AnalysisRequest):
    try:
        return _run_sync_engine(("limma", "limma_script.R"), request)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"limma-voom analysis failed: {str(e)}")

# ==============================================================================
# Enrichment / interpretation
# ==============================================================================
class EnrichmentRequest(BaseModel):
    genes: List[str]
    annotations: Dict[str, Any]
    universe: Optional[List[str]] = None
    min_term_size: int = 3
    max_term_size: int = 500

@app.post("/api/enrichment")
def run_enrichment(request: EnrichmentRequest):
    """Over-representation analysis of a gene list against supplied annotations."""
    results = enrich_genes(
        request.genes,
        request.annotations,
        universe=request.universe,
        min_term_size=request.min_term_size,
        max_term_size=request.max_term_size,
    )
    return {"count": len(results), "results": results}

# ==============================================================================
# Worker Endpoints
# ==============================================================================
security = HTTPBearer()

def verify_worker(credentials: HTTPAuthorizationCredentials = Security(security)):
    if not hmac.compare_digest(credentials.credentials, WORKER_TOKEN):
        raise HTTPException(status_code=401, detail="Invalid worker token")
    return credentials.credentials

class WorkerRegisterRequest(BaseModel):
    worker_id: str
    version: str

class WorkerHeartbeatRequest(BaseModel):
    worker_id: str
    current_job_id: Optional[str] = None

class WorkerResultsRequest(BaseModel):
    results: Dict[str, Any]
    warnings: Optional[List[str]] = []
    worker_id: str

class WorkerMethodRequest(BaseModel):
    worker_id: str
    status: str = "completed"
    results: Optional[List[Dict[str, Any]]] = None
    error: Optional[str] = None

class WorkerFinishRequest(BaseModel):
    worker_id: Optional[str] = None
    warnings: Optional[List[str]] = []

@app.post("/api/workers/register")
def register_worker(req: WorkerRegisterRequest, token: str = Depends(verify_worker)):
    store = JobStore()
    store.register_worker(req.worker_id, req.version)
    return {"status": "registered"}

@app.post("/api/workers/heartbeat")
def heartbeat_worker(req: WorkerHeartbeatRequest, token: str = Depends(verify_worker)):
    store = JobStore()
    store.heartbeat_worker(req.worker_id, req.current_job_id)
    return {"status": "ok"}

@app.get("/api/workers/jobs/next")
def get_next_job(token: str = Depends(verify_worker)):
    store = JobStore()
    job = store.get_next_available_job()
    if job:
        return {"job_id": job["job_id"], "params": job["params"]}
    return {"job_id": None}

@app.post("/api/workers/jobs/{job_id}/claim")
def claim_job(job_id: str, worker_id: str = Form(...), token: str = Depends(verify_worker)):
    store = JobStore()
    success = store.claim_job(job_id, worker_id)
    if not success:
        raise HTTPException(status_code=409, detail="Job already claimed or not queued")
    return {"status": "claimed"}

@app.get("/api/workers/jobs/{job_id}/dataset/{file_type}")
def download_dataset(job_id: str, file_type: str, token: str = Depends(verify_worker)):
    store = JobStore()
    job = store.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    
    if file_type == "counts":
        path = job["counts_path"] + ".clean.csv"
    elif file_type == "metadata":
        path = job["metadata_path"] + ".clean.csv"
    else:
        raise HTTPException(status_code=400, detail="Invalid file type")
        
    if not os.path.exists(path):
        raise HTTPException(status_code=404, detail="File not found")
        
    return FileResponse(path)

@app.post("/api/workers/jobs/{job_id}/results")
def upload_results(job_id: str, req: WorkerResultsRequest, token: str = Depends(verify_worker)):
    store = JobStore()
    job = store.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    if job["state"] == JobState.COMPLETED:
        raise HTTPException(status_code=409, detail="Results already submitted")
    if job["state"] == JobState.CANCELLED:
        raise HTTPException(status_code=409, detail="Job was cancelled; rejecting results")
    if job["state"] not in [JobState.RUNNING, JobState.GENERATING_RESULTS]:
        raise HTTPException(status_code=409, detail=f"Invalid state transition from {job['state']}")
    if job.get("worker_id") and req.worker_id != job["worker_id"]:
        raise HTTPException(status_code=409, detail="Worker does not own this job")

    success = store.update_job_results_atomic(
        job_id,
        results=req.results,
        warnings=req.warnings,
        worker_id=req.worker_id,
    )
    if not success:
        # State may have concurrently transitioned to CANCELLED or COMPLETED
        current_job = store.get_job(job_id)
        current_state = current_job.get("state") if current_job else "UNKNOWN"
        raise HTTPException(status_code=409, detail=f"Failed to submit results: job is now in {current_state} state")
        
    return {"status": "saved"}

@app.post("/api/workers/jobs/{job_id}/method/{engine}")
def upload_method_result(job_id: str, engine: str, req: WorkerMethodRequest, token: str = Depends(verify_worker)):
    """Persist one engine's status/result as soon as it finishes.

    This lets the browser show DESeq2 while edgeR is still running, and keeps a
    genuine edgeR failure separate from a successful DESeq2 result.
    """
    engine = engine.lower()
    if engine not in ENGINES:
        raise HTTPException(status_code=400, detail=f"Unknown engine '{engine}'.")
    if req.status not in ("running", "completed", "failed"):
        raise HTTPException(status_code=400, detail="status must be running, completed, or failed.")
    if req.status == "completed" and req.results is None:
        raise HTTPException(status_code=400, detail="Completed method uploads require a results table.")
    if req.status == "failed" and not req.error:
        raise HTTPException(status_code=400, detail="Failed method uploads require an error message.")

    job = JobStore().get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    updated = JobStore().update_job_method(
        job_id,
        engine,
        req.status,
        results=req.results,
        error=req.error,
        worker_id=req.worker_id,
    )
    if updated is None:
        raise HTTPException(status_code=409, detail="Job is finished, cancelled, or owned by another worker.")
    return {"status": req.status, "engine": engine, "job_state": updated["state"]}


@app.post("/api/workers/jobs/{job_id}/finish")
def finish_job(job_id: str, req: WorkerFinishRequest, token: str = Depends(verify_worker)):
    """Finalise per-method uploads without re-sending the large result tables."""
    job = JobStore().get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    state = JobStore().finish_job(job_id, worker_id=req.worker_id, warnings=req.warnings)
    if state is None:
        raise HTTPException(status_code=409, detail="Job is already finished or owned by another worker.")
    return {"status": "finished", "job_state": state}


@app.post("/api/workers/jobs/{job_id}/failed")
def mark_job_failed(job_id: str, error: str = Form(...), token: str = Depends(verify_worker)):
    store = JobStore()
    store.update_job(
        job_id, 
        JobState.FAILED, 
        0, 
        f"Worker failed: {error}", 
        error=error
    )
    return {"status": "updated"}
