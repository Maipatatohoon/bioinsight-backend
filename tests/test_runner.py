# backend/tests/test_runner.py
"""
Exercises REngineRunner against a fake `Rscript` on PATH, so the subprocess
handoff, output parsing and error handling can be tested without R installed.
"""
import os
import sys
import stat
import shutil
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.analysis.runner import REngineRunner

FAKE_OK = """#!/usr/bin/env bash
out=""
while [ $# -gt 0 ]; do case "$1" in --out) out="$2"; shift 2;; *) shift;; esac; done
printf 'Gene,log2FoldChange,baseMean,stat,pvalue,padj\\nA,1.2,10,3.1,0.01,0.05\\n' > "$out"
"""

FAKE_BAD = """#!/usr/bin/env bash
out=""
while [ $# -gt 0 ]; do case "$1" in --out) out="$2"; shift 2;; *) shift;; esac; done
printf 'nope\\n' > "$out"
"""

FAKE_EXIT = "#!/usr/bin/env bash\necho 'boom' 1>&2\nexit 3\n"


class TestREngineRunner(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="bioinsight_runner_")
        self.counts = os.path.join(self.tmp, "c.csv")
        self.meta = os.path.join(self.tmp, "m.csv")
        open(self.counts, "w").write("gene,S1,S2\nA,1,2\n")
        open(self.meta, "w").write("sample,condition\nS1,C\nS2,T\n")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _fake_rscript(self, body):
        path = os.path.join(self.tmp, "Rscript")
        with open(path, "w") as f:
            f.write(body)
        os.chmod(path, os.stat(path).st_mode | stat.S_IEXEC)
        return path

    def _run_with(self, body):
        self._fake_rscript(body)
        env = dict(os.environ, PATH=self.tmp + os.pathsep + os.environ.get("PATH", ""))
        with mock.patch.dict(os.environ, env, clear=True):
            return REngineRunner.run_r_script("script.R", self.counts, self.meta, {})

    def test_parses_output_and_keeps_gene_column(self):
        df = self._run_with(FAKE_OK)
        self.assertIn("Gene", df.columns)          # gene ids must survive
        self.assertEqual(df.iloc[0]["Gene"], "A")
        self.assertIn("padj", df.columns)

    def test_missing_columns_raise(self):
        with self.assertRaises(RuntimeError):
            self._run_with(FAKE_BAD)

    def test_nonzero_exit_raises_with_message(self):
        with self.assertRaises(RuntimeError) as ctx:
            self._run_with(FAKE_EXIT)
        self.assertIn("boom", str(ctx.exception))

    def test_missing_rscript_raises_clear_error(self):
        empty = os.path.join(self.tmp, "empty")
        os.makedirs(empty, exist_ok=True)
        with mock.patch.dict(os.environ, {"PATH": empty}, clear=True):
            with self.assertRaises(RuntimeError) as ctx:
                REngineRunner.run_r_script("script.R", self.counts, self.meta, {})
        self.assertIn("Rscript", str(ctx.exception))


if __name__ == "__main__":
    unittest.main()
