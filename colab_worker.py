# backend/colab_worker.py
"""
BioInsight 360 - Google Colab Computational Worker Script
Runs inside Google Colab (or any remote GPU/CPU node) as a dedicated worker.

Features:
1. Environment Setup (Installs R-base, DESeq2, edgeR, limma-voom, pyngrok, nest_asyncio)
2. Starts FastAPI backend on port 8000
3. Exposes secure public ngrok tunnel
4. Exposes job-oriented endpoints (/api/analysis, /api/analysis/{job_id}/results)
"""

import os
import sys
import subprocess

def setup_environment():
    print("=================================================================")
    print("🚀 BioInsight 360 Colab Worker Initialization")
    print("=================================================================")
    
    # 1. Install System & R packages
    print("📦 [1/3] Installing R-base & Bioconductor packages (DESeq2, edgeR, limma)...")
    cmd_r = (
        "apt-get update -qq && apt-get install -y -qq r-base && "
        "R -e \"install.packages('BiocManager', repos='http://cran.us.r-project.org')\" && "
        "R -e \"BiocManager::install(c('DESeq2', 'edgeR', 'limma', 'argparse'))\""
    )
    subprocess.run(cmd_r, shell=True, check=True)

    # 2. Install Python requirements
    print("📦 [2/3] Installing Python dependencies & pyngrok...")
    cmd_py = "pip install -q fastapi uvicorn pandas numpy pydantic pyngrok nest-asyncio"
    subprocess.run(cmd_py, shell=True, check=True)

    print("✅ Environment setup complete!")

def run_worker(ngrok_token: str, port: int = 8000):
    setup_environment()

    import nest_asyncio
    from pyngrok import ngrok
    import uvicorn

    print("🔒 [3/3] Authenticating ngrok and opening secure tunnel...")
    ngrok.set_auth_token(ngrok_token)

    public_url = ngrok.connect(port).public_url
    print("\n" + "="*66)
    print(f"🎉 BIOINSIGHT 360 WORKER IS LIVE!")
    print(f"🔗 Public API URL: {public_url}")
    print(f"👉 Copy this URL into the BioInsight 360 Frontend ('⚡ Colab Server')")
    print("="*66 + "\n")

    from app.main import app
    nest_asyncio.apply()
    uvicorn.run(app, host="0.0.0.0", port=port)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python colab_worker.py <YOUR_NGROK_AUTHTOKEN>")
        sys.exit(1)
    
    token = sys.argv[1]
    run_worker(token)
