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
    def run_r_script(script_path: str, count_matrix: pd.DataFrame, metadata: pd.DataFrame, params: Dict[str, Any]) -> pd.DataFrame:
        """
        1. Writes count_matrix and metadata to temp CSVs.
        2. Executes the R script.
        3. Reads the resulting CSV.
        4. Cleans up temp files.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            counts_path = os.path.join(tmpdir, "counts.csv")
            meta_path = os.path.join(tmpdir, "meta.csv")
            out_path = os.path.join(tmpdir, "results.csv")
            
            count_matrix.to_csv(counts_path)
            metadata.to_csv(meta_path)
            
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
