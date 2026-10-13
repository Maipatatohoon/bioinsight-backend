# backend/tests/test_concurrency.py
"""
Proves claim_job is exclusive even across separate processes (the real risk is
two workers on two machines, not two threads in one process).
"""
import os
import sys
import shutil
import tempfile
import multiprocessing as mp
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.jobs.job_manager import JobStore, JobState


def _try_claim(args):
    db_path, job_id, worker_id = args
    os.environ["BIOINSIGHT_DB_PATH"] = db_path
    JobStore._instance = None
    store = JobStore()
    return store.claim_job(job_id, worker_id)


@unittest.skipUnless("fork" in mp.get_all_start_methods(), "requires fork")
class TestConcurrentClaim(unittest.TestCase):
    def test_exactly_one_worker_claims(self):
        tmp = tempfile.mkdtemp(prefix="bioinsight_conc_")
        db_path = os.path.join(tmp, "jobs.db")
        os.environ["BIOINSIGHT_DB_PATH"] = db_path
        JobStore._instance = None
        try:
            store = JobStore()
            job_id = store.create_job("c.csv", "m.csv", {})
            store.update_job(job_id, JobState.QUEUED, 30, "queued")

            ctx = mp.get_context("fork")
            with ctx.Pool(6) as pool:
                results = pool.map(_try_claim, [(db_path, job_id, f"worker-{i}") for i in range(6)])

            self.assertEqual(sum(bool(r) for r in results), 1,
                             f"exactly one worker must win the claim, got {results}")
        finally:
            JobStore._instance = None
            os.environ.pop("BIOINSIGHT_DB_PATH", None)
            shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
