from fastapi import FastAPI, UploadFile, File, HTTPException
from fastapi.middleware.cors import CORSMiddleware
import pandas as pd
import io
from app.analysis.profiler import profile_dataset
from app.analysis.models import ProfilingReport

app = FastAPI(title="BioInsight 360 Backend", version="1.0.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"], # Update this in production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/")
def read_root():
    return {"message": "Welcome to BioInsight 360 Backend API"}

@app.get("/health")
def health_check():
    return {"status": "healthy"}

@app.post("/api/profile", response_model=ProfilingReport)
async def profile_upload(file: UploadFile = File(...)):
    if not file.filename.endswith(('.csv', '.tsv', '.txt')):
        raise HTTPException(status_code=400, detail="Only CSV or TSV files are supported.")
        
    try:
        contents = await file.read()
        
        # Determine separator
        sep = '\t' if file.filename.endswith(('.tsv', '.txt')) else ','
        
        df = pd.read_csv(io.StringIO(contents.decode('utf-8')), sep=sep)
        
        # Profile the dataset
        report = profile_dataset(df)
        return report
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error processing file: {str(e)}")

from pydantic import BaseModel
from typing import Dict, Any
from app.analysis.runner import REngineRunner
import os

class AnalysisRequest(BaseModel):
    counts: Dict[str, Dict[str, int]]
    metadata: Dict[str, Dict[str, str]]
    params: Dict[str, Any]

@app.post("/api/analyze/deseq2")
def run_deseq2(request: AnalysisRequest):
    try:
        counts_df = pd.DataFrame(request.counts)
        meta_df = pd.DataFrame(request.metadata)
        
        script_path = os.path.join(os.path.dirname(__file__), "deseq2", "deseq2_script.R")
        
        result_df = REngineRunner.run_r_script(script_path, counts_df, meta_df, request.params)
        
        # Replace NaN with None for JSON serialization
        result_df = result_df.where(pd.notnull(result_df), None)
        
        return {"status": "success", "results": result_df.to_dict(orient="index")}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"DESeq2 analysis failed: {str(e)}")

