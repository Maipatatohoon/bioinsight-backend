# backend/app/main.py
from fastapi import FastAPI, UploadFile, File, HTTPException, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import Response, StreamingResponse
from pydantic import BaseModel
from typing import Dict, Any, Optional
import pandas as pd
import io
import os

from app.analysis.profiler import profile_dataset
from app.analysis.models import ProfilingReport
from app.analysis.runner import REngineRunner
from app.jobs.job_manager import JobStore, JobState, process_job_async

app = FastAPI(
    title="BioInsight 360 Computational Engine",
    version="2.0.0",
    description="Scientific microservice API for RNA-seq differential expression (DESeq2, edgeR, limma-voom)."
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
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
def submit_analysis_job(request: AnalysisRequest, background_tasks: BackgroundTasks):
    store = JobStore()
    job_id = store.create_job(request.counts, request.metadata, request.params)
    background_tasks.add_task(process_job_async, job_id)
    return {
        "job_id": job_id,
        "state": JobState.QUEUED,
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
        "updated_at": job["updated_at"]
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
@app.post("/api/analyze/deseq2")
def run_deseq2_sync(request: AnalysisRequest):
    try:
        counts_df = pd.DataFrame(request.counts)
        meta_df = pd.DataFrame(request.metadata)
        script_path = os.path.join(os.path.dirname(__file__), "deseq2", "deseq2_script.R")
        result_df = REngineRunner.run_r_script(script_path, counts_df, meta_df, request.params)
        result_df = result_df.where(pd.notnull(result_df), None)
        return {"status": "success", "results": result_df.to_dict(orient="records")}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"DESeq2 analysis failed: {str(e)}")

@app.post("/api/analyze/edger")
def run_edger_sync(request: AnalysisRequest):
    try:
        counts_df = pd.DataFrame(request.counts)
        meta_df = pd.DataFrame(request.metadata)
        script_path = os.path.join(os.path.dirname(__file__), "edger", "edger_script.R")
        result_df = REngineRunner.run_r_script(script_path, counts_df, meta_df, request.params)
        result_df = result_df.where(pd.notnull(result_df), None)
        return {"status": "success", "results": result_df.to_dict(orient="records")}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"edgeR analysis failed: {str(e)}")

@app.post("/api/analyze/limma")
def run_limma_sync(request: AnalysisRequest):
    try:
        counts_df = pd.DataFrame(request.counts)
        meta_df = pd.DataFrame(request.metadata)
        script_path = os.path.join(os.path.dirname(__file__), "limma", "limma_script.R")
        result_df = REngineRunner.run_r_script(script_path, counts_df, meta_df, request.params)
        result_df = result_df.where(pd.notnull(result_df), None)
        return {"status": "success", "results": result_df.to_dict(orient="records")}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"limma-voom analysis failed: {str(e)}")
