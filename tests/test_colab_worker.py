import os
import sys
import unittest
from unittest.mock import patch, MagicMock

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from workers.colab_worker import process_job


class TestColabWorker(unittest.TestCase):
    def _session(self, state="RUNNING"):
        counts_resp = MagicMock(content=b"counts")
        meta_resp = MagicMock(content=b"meta")
        status_resp = MagicMock(status_code=200)
        status_resp.json.return_value = {"state": state}

        def get(url, **kwargs):
            if "dataset/counts" in url:
                return counts_resp
            if "dataset/metadata" in url:
                return meta_resp
            if "/api/analysis/" in url:
                return status_resp
            raise ValueError(f"Unexpected GET {url}")

        session = MagicMock()
        session.get.side_effect = get
        session.post.return_value = MagicMock(status_code=200)
        return session

    def _posts(self, session):
        return [call.args[0] for call in session.post.call_args_list]

    @patch("workers.colab_worker.run_r_script")
    @patch("workers.colab_worker.http_session")
    def test_process_job_uploads_each_method_then_finishes(self, mock_session, mock_run_r):
        session = self._session()
        mock_session.get = session.get
        mock_session.post = session.post
        mock_run_r.return_value = pd.DataFrame([
            {"Gene": "A", "log2FoldChange": 1.2, "pvalue": 0.01, "padj": 0.02}
        ])

        process_job("job123", {"design": "~ condition"})

        self.assertEqual(mock_run_r.call_count, 3)
        urls = self._posts(mock_session)
        for engine in ("deseq2", "edger", "limma"):
            self.assertTrue(any(f"/method/{engine}" in u for u in urls), urls)
        self.assertTrue(any(u.endswith("/finish") for u in urls), urls)
        self.assertFalse(any(u.endswith("/results") for u in urls), urls)

        completed = [
            c.kwargs["json"] for c in mock_session.post.call_args_list
            if c.args[0].endswith("/method/deseq2") and c.kwargs["json"]["status"] == "completed"
        ]
        self.assertEqual(completed[0]["results"][0]["Gene"], "A")

    @patch("workers.colab_worker.run_r_script")
    @patch("workers.colab_worker.http_session")
    def test_process_job_falls_back_to_aggregate_when_incremental_upload_fails(
        self, mock_session, mock_run_r
    ):
        session = self._session()
        mock_session.get = session.get

        def post(url, **kwargs):
            if url.endswith("/method/deseq2") and kwargs.get("json", {}).get("status") == "completed":
                return MagicMock(status_code=404)
            return MagicMock(status_code=200)

        mock_session.post.side_effect = post
        mock_run_r.return_value = pd.DataFrame([
            {"Gene": "A", "log2FoldChange": 1.2, "pvalue": 0.01, "padj": 0.02}
        ])

        process_job("job789", {"design": "~ condition"})

        urls = self._posts(mock_session)
        self.assertTrue(any(u.endswith("/results") for u in urls), urls)
        self.assertFalse(any(u.endswith("/finish") for u in urls), urls)

    @patch("workers.colab_worker.run_r_script")
    @patch("workers.colab_worker.http_session")
    def test_process_job_stops_when_cancelled(self, mock_session, mock_run_r):
        session = self._session(state="CANCELLED")
        mock_session.get = session.get
        mock_session.post = session.post

        process_job("job123", {"design": "~ condition"})

        self.assertEqual(mock_run_r.call_count, 0)
        mock_session.post.assert_not_called()

    @patch("workers.colab_worker.run_r_script")
    @patch("workers.colab_worker.http_session")
    def test_process_job_reports_partial_failure(self, mock_session, mock_run_r):
        session = self._session()
        mock_session.get = session.get
        mock_session.post = session.post
        good = pd.DataFrame([{"Gene": "A", "log2FoldChange": 1.0, "pvalue": 0.01, "padj": 0.02}])
        mock_run_r.side_effect = [RuntimeError("DESeq2 dispersion error"), good, good]

        process_job("job456", {"design": "~ condition"})

        calls = mock_session.post.call_args_list
        failed = [
            c.kwargs["json"] for c in calls
            if c.args[0].endswith("/method/deseq2") and c.kwargs["json"]["status"] == "failed"
        ]
        self.assertTrue(failed)
        self.assertIn("dispersion", failed[0]["error"])
        self.assertTrue(any(c.args[0].endswith("/finish") for c in calls))


if __name__ == "__main__":
    unittest.main()
