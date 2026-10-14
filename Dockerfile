# BioInsight 360 backend — FastAPI + R/Bioconductor (DESeq2, edgeR, limma)
FROM python:3.10-slim

LABEL org.opencontainers.image.title="BioInsight 360 backend" \
      org.opencontainers.image.description="Multi-method RNA-Seq consensus engine (DESeq2 / edgeR / limma)" \
      org.opencontainers.image.source="https://github.com/Maipatatohoon/bioinsight-backend"

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PORT=7860 \
    BIOINSIGHT_DB_PATH=/data/bioinsight_jobs.db

# System libraries for R and Bioconductor compilation
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libglpk-dev \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# R + Bioconductor packages (slow layer; cached)
RUN Rscript -e 'install.packages(c("BiocManager","argparse"), repos="https://cloud.r-project.org")' \
    && Rscript -e 'BiocManager::install(c("DESeq2","edgeR","limma"), ask=FALSE, update=FALSE)'

WORKDIR /app

COPY requirements.txt .
RUN pip install -r requirements.txt

COPY . .

# Run as a non-root user with a writable data volume
RUN mkdir -p /data && useradd -m -u 10001 appuser && chown -R appuser /app /data
USER appuser

EXPOSE 7860

HEALTHCHECK --interval=30s --timeout=5s --start-period=25s --retries=3 \
  CMD python -c "import urllib.request,sys; sys.exit(0 if urllib.request.urlopen('http://127.0.0.1:7860/health', timeout=3).status == 200 else 1)"

CMD ["sh", "-c", "uvicorn app.main:app --host 0.0.0.0 --port ${PORT:-7860}"]
