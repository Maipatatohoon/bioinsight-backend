# backend/tests/test_pure_python.py
import unittest
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.jobs.job_manager import JobStore, JobState

class TestJobManager(unittest.TestCase):

    def test_job_lifecycle(self):
        store = JobStore()
        job_id = store.create_job(
            counts={"Gene_1": {"Sample_1": 10}},
            metadata={"Sample_1": {"condition": "Control"}},
            params={"condition_col": "condition"}
        )
        self.assertIsNotNone(job_id)
        
        job = store.get_job(job_id)
        self.assertEqual(job["state"], JobState.QUEUED)
        self.assertEqual(job["progress"], 0)
        
        store.update_job(
            job_id,
            state=JobState.ANALYZING,
            progress=50,
            message="Processing DESeq2...",
            warnings=["Low count warning"]
        )
        
        updated = store.get_job(job_id)
        self.assertEqual(updated["state"], JobState.ANALYZING)
        self.assertEqual(updated["progress"], 50)
        self.assertIn("Low count warning", updated["warnings"])
        
        store.update_job(
            job_id,
            state=JobState.COMPLETED,
            progress=100,
            message="Done",
            results={"deseq2": [{"Gene": "Gene_1", "log2FoldChange": 1.5}]}
        )
        
        completed = store.get_job(job_id)
        self.assertEqual(completed["state"], JobState.COMPLETED)
        self.assertIn("deseq2", completed["results"])

if __name__ == "__main__":
    unittest.main()
