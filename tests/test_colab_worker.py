import unittest
from unittest.mock import patch, MagicMock
import sys
import os

# Add backend to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from workers.colab_worker import process_job, main, HeartbeatThread

class TestColabWorker(unittest.TestCase):
    @patch("workers.colab_worker.http_session")
    @patch("workers.colab_worker.run_r_script")
    def test_process_job_success(self, mock_run_r, mock_session):
        # Mock GET dataset counts/meta
        mock_counts_resp = MagicMock()
        mock_counts_resp.content = b"counts"
        mock_meta_resp = MagicMock()
        mock_meta_resp.content = b"meta"
        
        # Mock job status
        mock_status_resp = MagicMock()
        mock_status_resp.status_code = 200
        mock_status_resp.json.return_value = {"state": "RUNNING"}
        
        # Mock UPLOAD result
        mock_upload_resp = MagicMock()
        
        def mock_get(url, **kwargs):
            if "dataset/counts" in url:
                return mock_counts_resp
            elif "dataset/metadata" in url:
                return mock_meta_resp
            elif "/api/analysis/" in url:
                return mock_status_resp
            raise ValueError(f"Unexpected GET {url}")
            
        mock_session.get.side_effect = mock_get
        mock_session.post.return_value = mock_upload_resp
        
        # Mock R scripts
        mock_df = MagicMock()
        mock_df.where.return_value.to_dict.return_value = [{"gene": "A"}]
        mock_run_r.return_value = mock_df
        
        process_job("job123", {"design": "~ condition"})
        
        self.assertEqual(mock_run_r.call_count, 3) # DESeq2, edgeR, limma
        mock_session.post.assert_called_once()
        args, kwargs = mock_session.post.call_args
        self.assertIn("job123/results", args[0])
        self.assertIn("deseq2", kwargs["json"]["results"])
        
    @patch("workers.colab_worker.http_session")
    @patch("workers.colab_worker.run_r_script")
    def test_process_job_cancelled(self, mock_run_r, mock_session):
        # Mock GET dataset counts/meta
        mock_counts_resp = MagicMock()
        mock_counts_resp.content = b"counts"
        mock_meta_resp = MagicMock()
        mock_meta_resp.content = b"meta"
        
        # Mock job status (CANCELLED)
        mock_status_resp = MagicMock()
        mock_status_resp.status_code = 200
        mock_status_resp.json.return_value = {"state": "CANCELLED"}
        
        def mock_get(url, **kwargs):
            if "dataset/counts" in url:
                return mock_counts_resp
            elif "dataset/metadata" in url:
                return mock_meta_resp
            elif "/api/analysis/" in url:
                return mock_status_resp
            raise ValueError(f"Unexpected GET {url}")
            
        mock_session.get.side_effect = mock_get
        
        # Mock R scripts
        mock_df = MagicMock()
        mock_df.where.return_value.to_dict.return_value = [{"gene": "A"}]
        mock_run_r.return_value = mock_df
        
        process_job("job123", {"design": "~ condition"})
        
        self.assertEqual(mock_run_r.call_count, 3) # DESeq2, edgeR, limma
        mock_session.post.assert_not_called() # Should NOT upload results if cancelled

    @patch("workers.colab_worker.http_session")
    @patch("workers.colab_worker.run_r_script")
    def test_process_job_partial_failure(self, mock_run_r, mock_session):
        # Test that if one R engine fails, warning is recorded and remaining engines are uploaded
        mock_counts_resp = MagicMock()
        mock_counts_resp.content = b"counts"
        mock_meta_resp = MagicMock()
        mock_meta_resp.content = b"meta"
        mock_status_resp = MagicMock()
        mock_status_resp.status_code = 200
        mock_status_resp.json.return_value = {"state": "RUNNING"}
        
        def mock_get(url, **kwargs):
            if "dataset/counts" in url:
                return mock_counts_resp
            elif "dataset/metadata" in url:
                return mock_meta_resp
            elif "/api/analysis/" in url:
                return mock_status_resp
            raise ValueError(f"Unexpected GET {url}")
            
        mock_session.get.side_effect = mock_get
        
        # DESeq2 fails, edgeR succeeds, limma succeeds
        mock_df = MagicMock()
        mock_df.where.return_value.to_dict.return_value = [{"gene": "A"}]
        mock_run_r.side_effect = [RuntimeError("DESeq2 dispersion error"), mock_df, mock_df]
        
        process_job("job456", {"design": "~ condition"})
        
        self.assertEqual(mock_run_r.call_count, 3)
        mock_session.post.assert_called_once()
        args, kwargs = mock_session.post.call_args
        self.assertIn("job456/results", args[0])
        self.assertNotIn("deseq2", kwargs["json"]["results"])
        self.assertIn("edger", kwargs["json"]["results"])
        self.assertIn("limma", kwargs["json"]["results"])
        self.assertTrue(any("DESeq2" in w for w in kwargs["json"]["warnings"]))

if __name__ == "__main__":
    unittest.main()
