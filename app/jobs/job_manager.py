# backend/app/jobs/job_manager.py
"""
Asynchronous Job Manager for BioInsight 360 Backend.
Prevents browser timeouts on long-running RNA-seq analyses by tracking job state.
"""
import uuid
import time
import threading
import pandas as pd
import numpy as np
import os
import json
from typing import Dict, Any, Optional, List
from app.analysis.validator import validate_analysis_input, DataValidationError
from app.analysis.runner import REngineRunner

class JobState:
    QUEUED = "QUEUED"
    VALIDATING = "VALIDATING"
    ANALYZING = "ANALYZING"
    ENRICHMENT = "ENRICHMENT"
    GENERATING_RESULTS = "GENERATING_RESULTS"
    COMPLETED = "COMPLETED"
    FAILED = "FAILED"

class JobStore:
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(JobStore, cls).__new__(cls)
            cls._instance.jobs = {}
            cls._instance.lock = threading.Lock()
        return cls._instance

    def create_job(self, counts: Dict[str, Dict[str, Any]], metadata: Dict[str, Dict[str, Any]], params: Dict[str, Any]) -> str:
        job_id = str(uuid.uuid4())
        now = time.time()
        job_data = {
            "job_id": job_id,
            "state": JobState.QUEUED,
            "created_at": now,
            "updated_at": now,
            "progress": 0,
            "message": "Job queued for processing",
            "warnings": [],
            "error": None,
            "counts": counts,
            "metadata": metadata,
            "params": params,
            "results": None
        }
        with self.lock:
            self.jobs[job_id] = job_data
        return job_id

    def get_job(self, job_id: str) -> Optional[Dict[str, Any]]:
        with self.lock:
            return self.jobs.get(job_id)

    def update_job(self, job_id: str, state: str, progress: int, message: str, warnings: List[str] = None, error: str = None, results: Dict[str, Any] = None):
        with self.lock:
            if job_id in self.jobs:
                self.jobs[job_id]["state"] = state
                self.jobs[job_id]["progress"] = progress
                self.jobs[job_id]["message"] = message
                self.jobs[job_id]["updated_at"] = time.time()
                if warnings is not None:
                    self.jobs[job_id]["warnings"].extend(warnings)
                if error is not None:
                    self.jobs[job_id]["error"] = error
                if results is not None:
                    self.jobs[job_id]["results"] = results

def process_job_async(job_id: str):
    store = JobStore()
    job = store.get_job(job_id)
    if not job:
        return

    try:
        # Step 1: Validating
        store.update_job(job_id, JobState.VALIDATING, 10, "Validating input count matrix and metadata...")
        
        counts_df = pd.DataFrame(job["counts"])
        meta_df = pd.DataFrame(job["metadata"])
        params = job["params"]

        clean_counts, clean_meta, warnings = validate_analysis_input(counts_df, meta_df, params)

        # Step 2: Analyzing (DESeq2, edgeR, limma)
        store.update_job(job_id, JobState.ANALYZING, 30, "Running R-Bioconductor differential expression engines (DESeq2, edgeR, limma)...", warnings=warnings)

        base_dir = os.path.dirname(os.path.dirname(__file__))
        deseq2_script = os.path.join(base_dir, "deseq2", "deseq2_script.R")
        edger_script = os.path.join(base_dir, "edger", "edger_script.R")
        limma_script = os.path.join(base_dir, "limma", "limma_script.R")

        results_dict = {}

        # Run DESeq2
        try:
            store.update_job(job_id, JobState.ANALYZING, 45, "Executing DESeq2 dispersion estimation & Wald test...")
            deseq2_df = REngineRunner.run_r_script(deseq2_script, clean_counts, clean_meta, params)
            deseq2_df = deseq2_df.where(pd.notnull(deseq2_df), None)
            results_dict["deseq2"] = deseq2_df.to_dict(orient="records")
        except Exception as e:
            store.update_job(job_id, JobState.ANALYZING, 45, f"DESeq2 failed: {str(e)}", warnings=[f"DESeq2: {str(e)}"])

        # Run edgeR
        try:
            store.update_job(job_id, JobState.ANALYZING, 60, "Executing edgeR TMM normalization & Quasi-Likelihood F-test...")
            edger_df = REngineRunner.run_r_script(edger_script, clean_counts, clean_meta, params)
            edger_df = edger_df.where(pd.notnull(edger_df), None)
            results_dict["edger"] = edger_df.to_dict(orient="records")
        except Exception as e:
            store.update_job(job_id, JobState.ANALYZING, 60, f"edgeR failed: {str(e)}", warnings=[f"edgeR: {str(e)}"])

        # Run limma-voom
        try:
            store.update_job(job_id, JobState.ANALYZING, 75, "Executing limma-voom precision weighting & empirical Bayes...")
            limma_df = REngineRunner.run_r_script(limma_script, clean_counts, clean_meta, params)
            limma_df = limma_df.where(pd.notnull(limma_df), None)
            results_dict["limma"] = limma_df.to_dict(orient="records")
        except Exception as e:
            store.update_job(job_id, JobState.ANALYZING, 75, f"limma failed: {str(e)}", warnings=[f"limma: {str(e)}"])

        if not results_dict:
            raise RuntimeError("All differential expression engines (DESeq2, edgeR, limma) failed.")

        # Step 3: Generating Results & Consensus
        store.update_job(job_id, JobState.GENERATING_RESULTS, 90, "Building multi-tool consensus tables and summary statistics...")

        # Step 4: Completed
        store.update_job(job_id, JobState.COMPLETED, 100, "Analysis completed successfully!", results=results_dict)

    except DataValidationError as e:
        store.update_job(job_id, JobState.FAILED, 0, f"Scientific Validation Error: {e.message}", error=e.message)
    except Exception as e:
        store.update_job(job_id, JobState.FAILED, 0, f"Pipeline Error: {str(e)}", error=str(e))
