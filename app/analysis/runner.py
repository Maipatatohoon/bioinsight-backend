import subprocess
import shutil
import tempfile
import os
from typing import Dict, Any

import pandas as pd

# Every engine must emit these so the frontend can consume the result uniformly.
REQUIRED_COLUMNS = ["Gene", "log2FoldChange", "pvalue", "padj"]


class REngineRunner:
    """
    Orchestrates the execution of R scripts via subprocess. Subprocess + temp
    files is used instead of rpy2 — it is robust in containers and avoids
    rpy2's memory behaviour.
    """

    @staticmethod
    def run_r_script(script_path: str, counts_path: str, meta_path: str,
                     params: Dict[str, Any], timeout: int = 900) -> pd.DataFrame:
        """Run an R engine on the given count/metadata files and return its table.

        Raises RuntimeError with an actionable message on any failure so callers
        (the local runtime or the worker) can record it as a warning.
        """
        if shutil.which("Rscript") is None:
            raise RuntimeError(
                "Rscript was not found on PATH. Install R + Bioconductor, or run "
                "this engine on the worker (Google Colab)."
            )

        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, "results.csv")

            cmd = [
                "Rscript", script_path,
                "--counts", counts_path,
                "--meta", meta_path,
                "--out", out_path,
                "--design", params.get("design", "~ condition"),
                "--contrast", params.get("contrast", ""),
                "--numerator", str(params.get("treated_group", "")),
                "--denominator", str(params.get("control_group", "")),
            ]

            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=timeout)
            except subprocess.TimeoutExpired:
                raise RuntimeError(f"R script timed out after {timeout}s: {script_path}")
            except subprocess.CalledProcessError as e:
                detail = (e.stderr or e.stdout or "").strip()[-800:]
                raise RuntimeError(f"R script failed: {detail}")

            if not os.path.exists(out_path):
                raise RuntimeError("R script completed but output file is missing.")

            df = pd.read_csv(out_path)
            missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
            if missing:
                raise RuntimeError(
                    "R output is missing required column(s): " + ", ".join(missing)
                )
            return df
