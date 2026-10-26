"""Per-method result API and saved comparison tests (no R required)."""
import json
import os
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.analysis.comparison import compute_comparison
from app.jobs.job_manager import JobStore, JobState, ensure_methods


DESEQ2 = [
    {"Gene": "A", "log2FoldChange": 2.0, "baseMean": 100, "stat": 5, "pvalue": 0.001, "padj": 0.01},
    {"Gene": "B", "log2FoldChange": -1.5, "baseMean": 80, "stat": -4, "pvalue": 0.002, "padj": 0.02},
    {"Gene": "C", "log2FoldChange": 0.2, "baseMean": 50, "stat": 1, "pvalue": 0.5, "padj": 0.8},
    {"Gene": "D", "log2FoldChange": 1.2, "baseMean": 30, "stat": 3, "pvalue": 0.01, "padj": 0.04},
]
EDGER = [
    {"Gene": "A", "log2FoldChange": 1.8, "baseMean": 90, "stat": 4, "pvalue": 0.002, "padj": 0.02},
    {"Gene": "B", "log2FoldChange": 1.4, "baseMean": 70, "stat": 4, "pvalue": 0.003, "padj": 0.03},
    {"Gene": "C", "log2FoldChange": 0.1, "baseMean": 40, "stat": 0, "pvalue": 0.9, "padj": 0.95},
    {"Gene": "E", "log2FoldChange": 1.1, "baseMean": 20, "stat": 2, "pvalue": 0.02, "padj": 0.05},
]


class TestComparisonMath(unittest.TestCase):
    def test_shared_only_direction_and_consensus(self):
        result = compute_comparison(DESEQ2, EDGER, padj_threshold=0.05, log2fc_threshold=1.0)
        summary = result["summary"]
        self.assertEqual(summary["genes_compared"], 5)
        self.assertEqual(summary["shared_significant"], 2)   # A, B
        self.assertEqual(summary["direction_agree"], 1)      # A
        self.assertEqual(summary["direction_disagree"], 1)   # B (opposite signs)
        self.assertEqual(summary["consensus"], 1)
        self.assertEqual(summary["consensus_up"], 1)
        self.assertEqual(summary["deseq2_only"], 1)          # D
        self.assertEqual(summary["edger_only"], 1)           # E
        by_gene = {row["gene"]: row for row in result["rows"]}
        self.assertTrue(by_gene["A"]["consensus"])
        self.assertEqual(by_gene["B"]["direction_agreement"], "disagree")
        self.assertEqual(by_gene["D"]["direction_agreement"], "one_method")

    def test_zero_fold_change_is_not_a_direction_call(self):
        d = [{"Gene": "Z", "log2FoldChange": 0, "pvalue": 1e-9, "padj": 1e-8}]
        e = [{"Gene": "Z", "log2FoldChange": 0, "pvalue": 1e-9, "padj": 1e-8}]
        result = compute_comparison(d, e, padj_threshold=0.05, log2fc_threshold=0)
        row = result["rows"][0]
        self.assertFalse(row["deseq2_significant"])
        self.assertFalse(row["edger_significant"])
        self.assertFalse(row["consensus"])
        self.assertEqual(result["summary"]["direction_agree"], 0)

    def test_invalid_thresholds_rejected(self):
        with self.assertRaises(ValueError):
            compute_comparison(DESEQ2, EDGER, padj_threshold=0)
        with self.assertRaises(ValueError):
            compute_comparison(DESEQ2, EDGER, log2fc_threshold=-1)


class TestMethodApi(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="bioinsight_methods_")
        os.environ["BIOINSIGHT_DB_PATH"] = os.path.join(self.tmp, "jobs.db")
        JobStore._instance = None
        from app.main import app
        from starlette.testclient import TestClient
        self.client = TestClient(app)
        self.store = JobStore()
        self.job_id = self.store.create_job("c.csv", "m.csv", {"control_group": "C", "treated_group": "T"})
        self.store.update_job(self.job_id, JobState.QUEUED, 30, "queued")
        self.assertTrue(self.store.claim_job(self.job_id, "worker-a"))
        self.headers = {"Authorization": "Bearer dev-token-123"}

    def tearDown(self):
        JobStore._instance = None
        os.environ.pop("BIOINSIGHT_DB_PATH", None)
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _post_method(self, engine, status, results=None, error=None):
        payload = {"worker_id": "worker-a", "status": status}
        if results is not None:
            payload["results"] = results
        if error is not None:
            payload["error"] = error
        return self.client.post(
            f"/api/workers/jobs/{self.job_id}/method/{engine}",
            json=payload,
            headers=self.headers,
        )

    def test_partial_results_visible_while_edgeR_runs(self):
        self.assertEqual(self._post_method("deseq2", "running").status_code, 200)
        self.assertEqual(self._post_method("deseq2", "completed", DESEQ2).status_code, 200)
        self.assertEqual(self._post_method("edger", "running").status_code, 200)

        status = self.client.get(f"/api/analysis/{self.job_id}")
        self.assertEqual(status.status_code, 200)
        self.assertEqual(status.json()["state"], JobState.RUNNING)
        self.assertEqual(status.json()["methods"]["deseq2"]["status"], "completed")
        self.assertEqual(status.json()["methods"]["edger"]["status"], "running")

        partial = self.client.get(f"/api/analysis/{self.job_id}/results")
        self.assertEqual(partial.status_code, 200)
        self.assertEqual(partial.json()["status"], "partial")
        self.assertIn("deseq2", partial.json()["results"])
        self.assertNotIn("edger", partial.json()["results"])

        one = self.client.get(f"/api/analysis/{self.job_id}/method/deseq2")
        self.assertEqual(one.status_code, 200)
        self.assertEqual(len(one.json()["results"]), len(DESEQ2))
        self.assertEqual(self.client.get(f"/api/analysis/{self.job_id}/method/edger").status_code, 404)

        csv = self.client.get(f"/api/analysis/{self.job_id}/files?engine=deseq2")
        self.assertEqual(csv.status_code, 200)
        self.assertIn("Gene,log2FoldChange", csv.text)

    def test_empty_method_table_is_reported_as_failed(self):
        resp = self._post_method("deseq2", "completed", [])
        self.assertEqual(resp.status_code, 200, resp.text)
        status = self.client.get(f"/api/analysis/{self.job_id}").json()
        self.assertEqual(status["methods"]["deseq2"]["status"], "failed")
        self.assertIn("empty", status["methods"]["deseq2"]["error"].lower())

    def test_aggregate_upload_does_not_overwrite_a_saved_method_as_failed(self):
        self._post_method("deseq2", "completed", DESEQ2)
        aggregate = self.client.post(
            f"/api/workers/jobs/{self.job_id}/results",
            json={"worker_id": "worker-a", "results": {"edger": EDGER}, "warnings": []},
            headers=self.headers,
        )
        self.assertEqual(aggregate.status_code, 200, aggregate.text)
        status = self.client.get(f"/api/analysis/{self.job_id}").json()
        self.assertEqual(status["state"], JobState.COMPLETED)
        self.assertEqual(status["methods"]["deseq2"]["status"], "completed")
        self.assertEqual(status["methods"]["edger"]["status"], "completed")

    def test_finish_marks_partial_not_completed(self):
        self._post_method("deseq2", "completed", DESEQ2)
        self._post_method("edger", "failed", error="edgeR dispersion failure")
        finish = self.client.post(
            f"/api/workers/jobs/{self.job_id}/finish",
            json={"worker_id": "worker-a", "warnings": []},
            headers=self.headers,
        )
        self.assertEqual(finish.status_code, 200)
        self.assertEqual(finish.json()["job_state"], JobState.PARTIALLY_COMPLETED)

        status = self.client.get(f"/api/analysis/{self.job_id}").json()
        self.assertEqual(status["state"], JobState.PARTIALLY_COMPLETED)
        self.assertEqual(status["methods"]["edger"]["status"], "failed")
        self.assertIn("dispersion", status["error"])
        self.assertIn("deseq2", self.client.get(f"/api/analysis/{self.job_id}/results").json()["results"])

    def test_comparison_uses_saved_tables_and_persists(self):
        self._post_method("deseq2", "completed", DESEQ2)
        self._post_method("edger", "completed", EDGER)
        self.client.post(
            f"/api/workers/jobs/{self.job_id}/finish",
            json={"worker_id": "worker-a"},
            headers=self.headers,
        )
        created = self.client.post(
            f"/api/analysis/{self.job_id}/comparison",
            json={"padj_threshold": 0.05, "log2fc_threshold": 1.0},
        )
        self.assertEqual(created.status_code, 200, created.text)
        comparison = created.json()["comparison"]
        self.assertEqual(comparison["summary"]["consensus"], 1)
        self.assertEqual(comparison["summary"]["shared_significant"], 2)
        saved = self.client.get(f"/api/analysis/{self.job_id}/comparison")
        self.assertEqual(saved.status_code, 200)
        self.assertEqual(saved.json()["comparison"]["summary"]["consensus"], 1)
        csv = self.client.get(f"/api/analysis/{self.job_id}/comparison.csv")
        self.assertEqual(csv.status_code, 200)
        self.assertIn("deseq2_log2fc", csv.text)

    def test_finish_does_not_invent_a_limma_failure(self):
        self._post_method("deseq2", "completed", DESEQ2)
        self._post_method("edger", "completed", EDGER)
        finish = self.client.post(
            f"/api/workers/jobs/{self.job_id}/finish",
            json={"worker_id": "worker-a", "warnings": []},
            headers=self.headers,
        )
        self.assertEqual(finish.status_code, 200)
        self.assertEqual(finish.json()["job_state"], JobState.COMPLETED)
        status = self.client.get(f"/api/analysis/{self.job_id}").json()
        self.assertEqual(status["methods"]["limma"]["status"], "pending")
        self.assertEqual(status["methods"]["deseq2"]["status"], "completed")
        self.assertEqual(status["methods"]["edger"]["status"], "completed")

    def test_comparison_does_not_fall_back_to_raw_pvalue(self):
        raw_only = [{"Gene": "A", "log2FoldChange": 2.0, "pvalue": 0.001}]
        result = compute_comparison(raw_only, raw_only, padj_threshold=0.05, log2fc_threshold=1.0)
        self.assertEqual(result["summary"]["consensus"], 0)
        self.assertEqual(result["summary"]["shared_significant"], 0)
        self.assertFalse(result["rows"][0]["deseq2_significant"])

    def test_comparison_requires_both_methods(self):
        self._post_method("deseq2", "completed", DESEQ2)
        response = self.client.post(
            f"/api/analysis/{self.job_id}/comparison",
            json={"padj_threshold": 0.05, "log2fc_threshold": 1.0},
        )
        self.assertEqual(response.status_code, 409)
        self.assertIn("edger", response.json()["detail"].lower())

    def test_method_upload_rejects_other_worker(self):
        response = self.client.post(
            f"/api/workers/jobs/{self.job_id}/method/deseq2",
            json={"worker_id": "worker-b", "status": "completed", "results": DESEQ2},
            headers=self.headers,
        )
        self.assertEqual(response.status_code, 409)

    def test_legacy_job_without_methods_is_normalised(self):
        self.store.update_job(
            self.job_id, JobState.COMPLETED, 100, "legacy",
            results={"deseq2": DESEQ2, "edger": EDGER},
        )
        methods = ensure_methods(self.store.get_job(self.job_id))
        self.assertEqual(methods["deseq2"]["status"], "completed")
        self.assertEqual(methods["edger"]["rows"], len(EDGER))
        self.assertEqual(methods["limma"]["status"], "pending")


if __name__ == "__main__":
    unittest.main()
