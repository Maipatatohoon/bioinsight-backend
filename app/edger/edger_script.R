#!/usr/bin/env Rscript
# ==============================================================================
# edgeR Execution Script for BioInsight 360 Backend
# ==============================================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(argparse)
})

# 1. Parse Arguments
parser <- ArgumentParser(description="Run edgeR differential expression analysis.")
parser$add_argument("--counts", type="character", required=TRUE, help="Path to raw count matrix CSV")
parser$add_argument("--meta", type="character", required=TRUE, help="Path to metadata CSV")
parser$add_argument("--out", type="character", required=TRUE, help="Path to save output results CSV")
parser$add_argument("--design", type="character", default="~ condition", help="Design formula")
parser$add_argument("--contrast", type="character", default="", help="Contrast string")

args <- parser$parse_args()

# 2. Load Input Files
cts <- read.csv(args$counts, row.names=1, check.names=FALSE)
coldata <- read.csv(args$meta, row.names=1, check.names=FALSE)

# Ensure sample names match
common_samples <- intersect(colnames(cts), rownames(coldata))
if (length(common_samples) < 2) {
  stop("Fewer than 2 overlapping samples between count matrix and metadata.")
}
cts <- cts[, common_samples, drop=FALSE]
coldata <- coldata[common_samples, , drop=FALSE]

for (col in colnames(coldata)) {
  coldata[[col]] <- as.factor(coldata[[col]])
}

# 3. Create DGEList & Filter
dge <- DGEList(counts=as.matrix(cts), samples=coldata)
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# 4. TMM Normalization & Dispersion Estimation
dge <- calcNormFactors(dge, method="TMM")
design_formula <- as.formula(args$design)
design_mat <- model.matrix(design_formula, data=coldata)

dge <- estimateDisp(dge, design_mat)
fit <- glmQLFit(dge, design_mat)

# 5. Contrast Testing
if (nchar(args$contrast) > 0) {
  parts <- unlist(strsplit(args$contrast, "_vs_"))
  if (length(parts) == 2) {
    sub_parts <- unlist(strsplit(parts[1], "_"))
    factor_name <- sub_parts[1]
    num_level <- paste(sub_parts[-1], collapse="_")
    den_level <- parts[2]
    
    coef_name <- paste0(factor_name, num_level)
    if (coef_name %in% colnames(design_mat)) {
      qlf <- glmQLFTest(fit, coef=coef_name)
    } else {
      qlf <- glmQLFTest(fit, coef=2)
    }
  } else {
    qlf <- glmQLFTest(fit, coef=2)
  }
} else {
  qlf <- glmQLFTest(fit, coef=2)
}

# 6. Format Output Table
res_table <- topTags(qlf, n=Inf)$table
res_table$Gene <- rownames(res_table)

# Map edgeR columns to standard BioInsight 360 backend schema
res_df <- data.frame(
  Gene = res_table$Gene,
  log2FoldChange = res_table$logFC,
  baseMean = 2^(res_table$logCPM),
  stat = res_table$F,
  pvalue = res_table$PValue,
  padj = res_table$FDR,
  stringsAsFactors = FALSE
)

write.csv(res_df, file=args$out, row.names=FALSE)
cat("edgeR execution completed successfully.\n")
