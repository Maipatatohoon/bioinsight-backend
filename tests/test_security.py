# backend/tests/test_security.py
"""
Config-level security tests: configurable CORS, production-safe token handling,
constant-time comparison. These reload app.main because configuration is read at
import time.
"""
import os
import sys
import importlib
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

_ENV_KEYS = ("BIOINSIGHT_ALLOWED_ORIGINS", "BIOINSIGHT_ENV", "BIOINSIGHT_WORKER_TOKEN")


class BackendConfigTest(unittest.TestCase):
    def setUp(self):
        self._saved = {k: os.environ.get(k) for k in _ENV_KEYS}

    def tearDown(self):
        for k, v in self._saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
        # Restore the default configuration for the rest of the suite.
        self._reload()

    def _reload(self):
        import app.main
        return importlib.reload(app.main)

    def test_allowed_origins_parsed_from_env(self):
        os.environ["BIOINSIGHT_ALLOWED_ORIGINS"] = "https://a.example, https://b.example ,"
        m = self._reload()
        self.assertEqual(m._allowed_origins(), ["https://a.example", "https://b.example"])

    def test_allowed_origins_default_is_not_wildcard(self):
        os.environ.pop("BIOINSIGHT_ALLOWED_ORIGINS", None)
        m = self._reload()
        origins = m._allowed_origins()
        self.assertIn("http://localhost:3000", origins)
        self.assertNotIn("*", origins)

    def test_allowed_origins_star(self):
        os.environ["BIOINSIGHT_ALLOWED_ORIGINS"] = "*"
        m = self._reload()
        self.assertEqual(m._allowed_origins(), ["*"])
        self.assertEqual(m._CORS_ORIGINS, ["*"])

    def test_production_refuses_default_token(self):
        os.environ["BIOINSIGHT_ENV"] = "production"
        os.environ.pop("BIOINSIGHT_WORKER_TOKEN", None)
        with self.assertRaises(RuntimeError):
            self._reload()

    def test_production_accepts_real_token(self):
        os.environ["BIOINSIGHT_ENV"] = "production"
        os.environ["BIOINSIGHT_WORKER_TOKEN"] = "s3cret-token"
        m = self._reload()
        self.assertEqual(m.WORKER_TOKEN, "s3cret-token")

    def test_token_verification_accepts_and_rejects(self):
        os.environ["BIOINSIGHT_ENV"] = "development"
        os.environ["BIOINSIGHT_WORKER_TOKEN"] = "abc123"
        m = self._reload()
        creds_ok = m.HTTPAuthorizationCredentials(scheme="Bearer", credentials="abc123")
        creds_bad = m.HTTPAuthorizationCredentials(scheme="Bearer", credentials="nope")
        self.assertEqual(m.verify_worker(creds_ok), "abc123")
        with self.assertRaises(m.HTTPException):
            m.verify_worker(creds_bad)


if __name__ == "__main__":
    unittest.main()
