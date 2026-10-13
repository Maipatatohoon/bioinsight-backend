import subprocess
import pandas as pd
import tempfile
import os
from typing import Dict, Any

class REngineRunner:
    """
    Orchestrates the execution of R scripts via subprocess, avoiding heavy rpy2 memory leaks
    if used indiscriminately. For deep integration, rpy2 is also available, but subprocess
    with temporary files is extremely robust for containerized pipelines.
    """
    
    @staticmethod
    def run_r_script(script_path: str, counts_path: str, meta_path: str, params: Dict[str, Any]) -> pd.DataFrame:
        """
        1. Executes the R script directly on the uploaded files.
        2. Reads the resulting CSV.
        3. Cleans up temp output file.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, "results.csv")
            
            # Construct R command
            # Assuming the R script takes arguments: --counts, --meta, --out, etc.
            cmd = [
                "Rscript", script_path,
                "--counts", counts_path,
                "--meta", meta_path,
                "--out", out_path,
                "--design", params.get("design", "~ condition"),
                "--contrast", params.get("contrast", "")
            ]
            
            try:
                result = subprocess.run(cmd, check=True, capture_output=True, text=True)
                # print(result.stdout)
                
                if os.path.exists(out_path):
                    return pd.read_csv(out_path, index_col=0)
                else:
                    raise RuntimeError("R script completed but output file is missing.")
                    
            except subprocess.CalledProcessError as e:
                raise RuntimeError(f"R script failed: {e.stderr}")
