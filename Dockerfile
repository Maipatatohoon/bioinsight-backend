FROM python:3.10-slim

# Install system dependencies for R and C/C++ compilation
RUN apt-get update && apt-get install -y \
    r-base \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install BiocManager and core Bioconductor packages
RUN Rscript -e 'install.packages(c("BiocManager", "argparse"), repos="http://cran.rstudio.com/")' \
    && Rscript -e 'BiocManager::install(c("DESeq2", "edgeR", "limma"))'

# Set working directory
WORKDIR /app

# Install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of the application
COPY . .

# Expose port (Hugging Face Spaces default is 7860)
EXPOSE 7860

# Start FastAPI server
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "7860"]
