# backend/tests/test_local_runtime.py
"""
LocalRuntime runs the R engines in-process (subprocess Rscript) so the full job
pipeline can be exercised without a Colab worker. R itself is mocked here.
"""
import os
import sys
import shutil
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pandas as pd

from app.jobs.job_manager import JobStore, JobState, LocalRuntime


class TestLocalRuntime(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="bioinsight_local_")
        os.environ["BIOINSIGHT_DB_PATH"] = os.path.join(self.tmp, "jobs.db")
        JobStore._instance = None
        self.store = JobStore()

        self.counts = os.path.join(self.tmp, "counts.csv")
        self.meta = os.path.join(self.tmp, "meta.csv")
        with open(self.counts + ".clean.csv", "w") as f:
            f.write("gene,S1,S2,S3,S4\nA,10,12,11,9\nB,5,6,7,8\n")
        with open(self.meta + ".clean.csv", "w") as f:
            f.write("sample,condition\nS1,Control\nS2,Control\nS3,Treated\nS4,Treated\n")

        self.job_id = self.store.create_job(self.counts, self.meta, {"condition_col": "condition"})
        self.store.update_job(self.job_id, JobState.QUEUED, 30, "queued")

    def tearDown(self):
        JobStore._instance = None
        os.environ.pop("BIOINSIGHT_DB_PATH", None)
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_dispatch_runs_engines_and_completes(self):
        df = pd.DataFrame([{"Gene": "A", "log2FoldChange": 1.0, "padj": 0.01}])
        with mock.patch("app.analysis.runner.REngineRunner.run_r_script", return_value=df) as m:
            LocalRuntime().dispatch(self.store, self.job_id)
        job = self.store.get_job(self.job_id)
        self.assertEqual(job["state"], JobState.COMPLETED)
        self.assertEqual(m.call_count, 3)
        self.assertIn("deseq2", job["results"])
        self.assertIn("edger", job["results"])
        self.assertIn("limma", job["results"])

    def test_dispatch_keeps_partial_results(self):
        good = pd.DataFrame([{"Gene": "A", "log2FoldChange": 1.0, "padj": 0.01}])
        with mock.patch(
            "app.analysis.runner.REngineRunner.run_r_script",
            side_effect=[RuntimeError("DESeq2 failed"), good, good],
        ):
            LocalRuntime().dispatch(self.store, self.job_id)
        job = self.store.get_job(self.job_id)
        self.assertEqual(job["state"], JobState.PARTIALLY_COMPLETED)
        self.assertNotIn("deseq2", job["results"])
        self.assertIn("edger", job["results"])
        self.assertTrue(any("deseq2" in w.lower() for w in job["warnings"]))

    def test_dispatch_fails_when_all_engines_fail(self):
        with mock.patch("app.analysis.runner.REngineRunner.run_r_script", side_effect=RuntimeError("boom")):
            LocalRuntime().dispatch(self.store, self.job_id)
        job = self.store.get_job(self.job_id)
        self.assertEqual(job["state"], JobState.FAILED)
        self.assertIsNotNone(job["error"])


if __name__ == "__main__":
    unittest.main()
