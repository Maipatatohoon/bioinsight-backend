# backend/tests/test_backend.py
import unittest
import sys
import os
import shutil
import tempfile
import pandas as pd
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.analysis.validator import validate_analysis_input, DataValidationError
from app.jobs.job_manager import JobStore, JobState
from app.analysis.profiler import profile_dataset

class TestBioInsightBackend(unittest.TestCase):

    def setUp(self):
        self._tmp = tempfile.mkdtemp(prefix="bioinsight_be_")
        os.environ["BIOINSIGHT_DB_PATH"] = os.path.join(self._tmp, "jobs.db")
        JobStore._instance = None
        np.random.seed(42)
        genes = [f"Gene_{i}" for i in range(1, 21)]
        samples = ["Control_1", "Control_2", "Treated_1", "Treated_2"]
        data = np.random.negative_binomial(n=10, p=0.3, size=(20, 4))
        self.counts_df = pd.DataFrame(data, index=genes, columns=samples)
        self.meta_df = pd.DataFrame({
            "condition": ["Control", "Control", "Treated", "Treated"],
            "organism": ["Homo sapiens"] * 4
        }, index=samples)

    def tearDown(self):
        JobStore._instance = None
        os.environ.pop("BIOINSIGHT_DB_PATH", None)
        shutil.rmtree(self._tmp, ignore_errors=True)

    def test_valid_input_validation(self):
        clean_counts, clean_meta, warnings = validate_analysis_input(self.counts_df, self.meta_df, {"condition_col": "condition"})
        self.assertEqual(len(clean_counts), 20)
        self.assertEqual(len(clean_counts.columns), 4)

    def test_negative_counts_validation(self):
        self.counts_df.iloc[0, 0] = -5
        with self.assertRaises(DataValidationError) as context:
            validate_analysis_input(self.counts_df, self.meta_df, {"condition_col": "condition"})
        self.assertIn("contains negative values", str(context.exception))

    def test_insufficient_replicates_validation(self):
        self.meta_df.iloc[3, 0] = "OtherGroup"
        with self.assertRaises(DataValidationError) as context:
            validate_analysis_input(self.counts_df, self.meta_df, {"condition_col": "condition"})
        self.assertIn("require at least 2 biological replicates", str(context.exception))

    def test_job_store_creation(self):
        store = JobStore()
        job_id = store.create_job("/tmp/c.csv", "/tmp/m.csv", {})
        self.assertIsNotNone(job_id)
        job = store.get_job(job_id)
        self.assertEqual(job["state"], JobState.VALIDATING)
        
        store.update_job(job_id, JobState.COMPLETED, 100, "Finished", results={"status": "ok"})
        job_updated = store.get_job(job_id)
        self.assertEqual(job_updated["state"], JobState.COMPLETED)
        self.assertEqual(job_updated["results"], {"status": "ok"})

    def test_dataset_profiler(self):
        df = pd.DataFrame({
            "Gene": ["ACTB", "GAPDH", "TP53"],
            "log2FoldChange": [1.5, -0.8, 2.1],
            "pvalue": [0.001, 0.04, 0.0001],
            "padj": [0.005, 0.09, 0.001]
        })
        report = profile_dataset(df)
        self.assertGreater(len(report.analysis_units), 0)
        self.assertFalse(report.analysis_units[0].is_raw_counts)

if __name__ == "__main__":
    unittest.main()
