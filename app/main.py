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
import tempfile

from app.analysis.profiler import profile_dataset
from app.analysis.models import ProfilingReport
from app.analysis.runner import REngineRunner
from app.analysis.enrichment import enrich_genes
from app.jobs.job_manager import JobStore, JobState, process_job_async
from app.storage import save_upload
import json

app = FastAPI(
    title="BioInsight 360 Computational Engine",
    version="2.0.0",
    description="Scientific microservice API for RNA-seq differential expression (DESeq2, edgeR, limma-voom)."
)

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
    (comma-separated); defaults to local development, never '*'."""
    raw = os.environ.get("BIOINSIGHT_ALLOWED_ORIGINS", "")
    origins = [o.strip() for o in raw.split(",") if o.strip()]
    return origins or ["http://localhost:3000", "http://localhost:8000"]


app.add_middleware(
    CORSMiddleware,
    allow_origins=_allowed_origins(),
    allow_credentials=True,
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
    counts_path = save_upload(counts_file, "counts")
    metadata_path = save_upload(metadata_file, "metadata")
    
    parsed_params = {}
    if params:
        try:
            parsed_params = json.loads(params)
        except json.JSONDecodeError:
            raise HTTPException(status_code=400, detail="Invalid JSON in params")
            
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
    
    return {
        "job_id": job["job_id"],
        "state": job["state"],
        "progress": job["progress"],
        "message": job["message"],
        "warnings": job["warnings"],
        "error": job["error"],
        "created_at": job["created_at"],
        "updated_at": job["updated_at"],
        "qc_data": job.get("qc_data")
    }

@app.get("/api/analysis/{job_id}/results")
def get_job_results(job_id: str):
    store = JobStore()
    job = store.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job ID '{job_id}' not found.")
    if job["state"] != JobState.COMPLETED:
        raise HTTPException(status_code=400, detail=f"Job '{job_id}' is not complete yet. Current state: {job['state']}")
    
    return {
        "job_id": job_id,
        "status": "success",
        "warnings": job["warnings"],
        "results": job["results"]
    }

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
    if not job or job["state"] != JobState.COMPLETED:
        raise HTTPException(status_code=404, detail=f"Completed job '{job_id}' not found.")
    
    results = job.get("results", {})
    engine_data = results.get(engine.lower())
    if not engine_data:
        raise HTTPException(status_code=404, detail=f"Engine '{engine}' results not available for job '{job_id}'.")
    
    df = pd.DataFrame(engine_data)
    stream = io.StringIO()
    df.to_csv(stream, index=False)
    
    return Response(
        content=stream.getvalue(),
        media_type="text/csv",
        headers={"Content-Disposition": f"attachment; filename=bioinsight_{engine}_results_{job_id[:8]}.csv"}
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
