# backend/tests/test_fixes.py
"""
Regression tests for the audit fixes. These are hermetic (temp SQLite DB via
BIOINSIGHT_DB_PATH) and must fail against the pre-fix code.
"""
import os
import sys
import gc
import shutil
import tempfile
import warnings
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.jobs import job_manager
from app.jobs.job_manager import JobStore, JobState


class HermeticJobStoreTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="bioinsight_test_")
        os.environ["BIOINSIGHT_DB_PATH"] = os.path.join(self.tmp, "jobs.db")
        JobStore._instance = None

    def tearDown(self):
        JobStore._instance = None
        os.environ.pop("BIOINSIGHT_DB_PATH", None)
        shutil.rmtree(self.tmp, ignore_errors=True)

    # --- Fix 2: configurable DB path -------------------------------------
    def test_db_path_is_configurable(self):
        store = JobStore()
        self.assertEqual(store._db_path, os.environ["BIOINSIGHT_DB_PATH"])
        self.assertTrue(os.path.exists(store._db_path))

    # --- Fix 1: CANCELLED exists, cancel works ---------------------------
    def test_cancel_running_job_sets_cancelled(self):
        store = JobStore()
        job_id = store.create_job("c.csv", "m.csv", {})
        store.update_job(job_id, JobState.QUEUED, 30, "queued")
        self.assertTrue(store.cancel_job(job_id))
        self.assertEqual(store.get_job(job_id)["state"], JobState.CANCELLED)

    def test_cancel_completed_job_returns_false(self):
        store = JobStore()
        job_id = store.create_job("c.csv", "m.csv", {})
        store.update_job(job_id, JobState.COMPLETED, 100, "done")
        self.assertFalse(store.cancel_job(job_id))

    # --- Fix 3: connections are closed -----------------------------------
    def test_operations_do_not_leak_connections(self):
        store = JobStore()
        store.create_job("c.csv", "m.csv", {})  # warm up
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            for _ in range(60):
                job_id = store.create_job("c.csv", "m.csv", {})
                store.get_job(job_id)
            gc.collect()
        leaked = [w for w in caught if issubclass(w.category, ResourceWarning)]
        self.assertEqual(leaked, [], f"unclosed sqlite connections: {leaked[:2]}")

    # --- Fix 4: claim is exclusive ---------------------------------------
    def test_second_claim_returns_false(self):
        store = JobStore()
        job_id = store.create_job("c.csv", "m.csv", {})
        store.update_job(job_id, JobState.QUEUED, 30, "queued")
        self.assertTrue(store.claim_job(job_id, "worker-a"))
        self.assertFalse(store.claim_job(job_id, "worker-b"))
        self.assertEqual(store.get_job(job_id)["worker_id"], "worker-a")


class ApiWorkflowTest(unittest.TestCase):
    """End-to-end HTTP flow that does NOT require R: submit -> validate/QC -> QUEUED."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="bioinsight_api_")
        os.environ["BIOINSIGHT_DB_PATH"] = os.path.join(self.tmp, "jobs.db")
        JobStore._instance = None
        from app.main import app
        from starlette.testclient import TestClient
        self.client = TestClient(app)
        self.counts = (
            "gene,S1,S2,S3,S4\n"
            "A,100,120,110,90\n"
            "B,50,60,70,80\n"
            "C,200,190,220,210\n"
            "D,10,12,11,9\n"
        )
        self.meta = "sample,condition\nS1,Control\nS2,Control\nS3,Treated\nS4,Treated\n"

    def tearDown(self):
        JobStore._instance = None
        os.environ.pop("BIOINSIGHT_DB_PATH", None)
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _submit(self):
        return self.client.post(
            "/api/analysis",
            files={
                "counts_file": ("counts.csv", self.counts, "text/csv"),
                "metadata_file": ("meta.csv", self.meta, "text/csv"),
            },
            data={"params": "{}"},
        )

    # --- Fix 6: response state matches stored state ----------------------
    def test_submit_reports_real_state(self):
        resp = self._submit()
        self.assertEqual(resp.status_code, 202)
        self.assertEqual(resp.json()["state"], JobState.VALIDATING)

    def test_async_pipeline_reaches_queued_with_qc(self):
        resp = self._submit()
        job_id = resp.json()["job_id"]
        status = self.client.get(f"/api/analysis/{job_id}").json()
        self.assertEqual(status["state"], JobState.QUEUED)
        self.assertIsNotNone(status["qc_data"])

    # --- Fix 5: sync routes hand file paths (not DataFrames) to runner ---
    def test_sync_route_passes_file_paths_to_runner(self):
        from unittest import mock
        import pandas as pd
        from app.analysis.runner import REngineRunner

        captured = {}

        def fake_run(script_path, counts_path, meta_path, params):
            captured["counts_path"] = counts_path
            captured["meta_path"] = meta_path
            captured["existed"] = os.path.exists(counts_path) and os.path.exists(meta_path)
            return pd.DataFrame([{"Gene": "A", "log2FoldChange": 1.0}])

        payload = {
            "counts": {"A": {"S1": 10, "S2": 12, "S3": 11, "S4": 9}},
            "metadata": {
                "S1": {"condition": "Control"},
                "S2": {"condition": "Control"},
                "S3": {"condition": "Treated"},
                "S4": {"condition": "Treated"},
            },
            "params": {},
        }
        with mock.patch.object(REngineRunner, "run_r_script", side_effect=fake_run):
            resp = self.client.post("/api/analyze/deseq2", json=payload)

        self.assertEqual(resp.status_code, 200, resp.text)
        self.assertIsInstance(captured["counts_path"], str)
        self.assertTrue(captured["existed"], "runner must receive live file paths, not DataFrames")

    # --- Fix 7: result upload is bound to the claiming worker ------------
    def test_results_rejected_from_non_owning_worker(self):
        store = JobStore()
        job_id = store.create_job("c.csv", "m.csv", {})
        store.update_job(job_id, JobState.QUEUED, 30, "queued")
        store.claim_job(job_id, "worker-a")
        headers = {"Authorization": "Bearer dev-token-123"}
        resp = self.client.post(
            f"/api/workers/jobs/{job_id}/results",
            json={"worker_id": "worker-b", "results": {"deseq2": []}, "warnings": []},
            headers=headers,
        )
        self.assertEqual(resp.status_code, 409, resp.text)


if __name__ == "__main__":
    unittest.main()
