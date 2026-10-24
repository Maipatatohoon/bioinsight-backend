# backend/tests/test_default_db_path.py
"""Regression: JobStore must not depend on BIOINSIGHT_DB_PATH being set.

`_default_db_path` was a @staticmethod that referenced `cls`, so instantiating
JobStore without the env var raised NameError -> HTTP 500 on every
POST /api/analysis. Colab and Docker set the variable, which hid the bug; a
plain `uvicorn app.main:app` did not.
"""
import os
import shutil
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.jobs.job_manager import JobStore


class DefaultDbPathTest(unittest.TestCase):
    def setUp(self):
        self._tmp = None
        JobStore._instance = None
        os.environ.pop("BIOINSIGHT_DB_PATH", None)

    def tearDown(self):
        JobStore._instance = None
        os.environ.pop("BIOINSIGHT_DB_PATH", None)
        # Remove the DB this test creates under backend/data, if it was not
        # there before, so the working tree stays clean.
        if self._tmp and os.path.exists(self._tmp):
            os.remove(self._tmp)

    def test_job_store_builds_a_default_path_without_env_var(self):
        store = JobStore()
        self.assertTrue(store._db_path.endswith("bioinsight_jobs.db"))
        self.assertTrue(os.path.isabs(store._db_path))
        self._tmp = store._db_path

    def test_env_var_still_wins(self):
        import tempfile
        tmpdir = tempfile.mkdtemp(prefix="bioinsight_dbpath_")
        try:
            custom = os.path.join(tmpdir, "custom.db")
            os.environ["BIOINSIGHT_DB_PATH"] = custom
            store = JobStore()
            self.assertEqual(store._db_path, custom)
        finally:
            JobStore._instance = None
            os.environ.pop("BIOINSIGHT_DB_PATH", None)
            shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
