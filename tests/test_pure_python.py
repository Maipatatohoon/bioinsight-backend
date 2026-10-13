# backend/tests/test_pure_python.py
import unittest
import sys
import os
import shutil
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.jobs.job_manager import JobStore, JobState

class TestJobManager(unittest.TestCase):

    def setUp(self):
        self._tmp = tempfile.mkdtemp(prefix="bioinsight_pure_")
        os.environ["BIOINSIGHT_DB_PATH"] = os.path.join(self._tmp, "jobs.db")
        JobStore._instance = None

    def tearDown(self):
        JobStore._instance = None
        os.environ.pop("BIOINSIGHT_DB_PATH", None)
        shutil.rmtree(self._tmp, ignore_errors=True)

    def test_job_lifecycle(self):
        store = JobStore()
        job_id = store.create_job(
            counts_path="/tmp/fake_counts.csv",
            metadata_path="/tmp/fake_meta.csv",
            params={"condition_col": "condition"}
        )
        self.assertIsNotNone(job_id)
        
        job = store.get_job(job_id)
        self.assertEqual(job["state"], JobState.VALIDATING)
        self.assertEqual(job["progress"], 0)
        
        # Transition to QUEUED as ColabRuntime does
        store.update_job(job_id, JobState.QUEUED, 30, "Ready for worker")

        # Test transition to RUNNING via claim
        claimed_success = store.claim_job(job_id, "worker-1")
        self.assertTrue(claimed_success)
        claimed = store.get_job(job_id)
        self.assertEqual(claimed["state"], JobState.RUNNING)
        self.assertEqual(claimed["worker_id"], "worker-1")
        
        # Test atomic result completion
        results_payload = {"deseq2": [{"Gene": "Gene_1", "log2FoldChange": 1.5}]}
        success = store.update_job_results_atomic(job_id, results_payload, warnings=["Low count warning"])
        self.assertTrue(success)
        
        completed = store.get_job(job_id)
        self.assertEqual(completed["state"], JobState.COMPLETED)
        self.assertEqual(completed["progress"], 100)
        self.assertIn("deseq2", completed["results"])
        self.assertIn("Low count warning", completed["warnings"])

        # Test duplicate result submission is rejected
        dup_success = store.update_job_results_atomic(job_id, results_payload)
        self.assertFalse(dup_success)

if __name__ == "__main__":
    unittest.main()
