#!/usr/bin/env Rscript
# ==============================================================================
# limma-voom Execution Script for BioInsight 360 Backend
# ==============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(argparse)
})

# 1. Parse Arguments
parser <- ArgumentParser(description="Run limma-voom differential expression analysis.")
parser$add_argument("--counts", type="character", required=TRUE, help="Path to raw count matrix CSV")
parser$add_argument("--meta", type="character", required=TRUE, help="Path to metadata CSV")
parser$add_argument("--out", type="character", required=TRUE, help="Path to save output results CSV")
parser$add_argument("--design", type="character", default="~ condition", help="Design formula")
parser$add_argument("--contrast", type="character", default="", help="Contrast string")

args <- parser$parse_args()

# 2. Load Input Files
cts <- read.csv(args$counts, row.names=1, check.names=FALSE)
coldata <- read.csv(args$meta, row.names=1, check.names=FALSE)

common_samples <- intersect(colnames(cts), rownames(coldata))
if (length(common_samples) < 2) {
  stop("Fewer than 2 overlapping samples between count matrix and metadata.")
}
cts <- cts[, common_samples, drop=FALSE]
coldata <- coldata[common_samples, , drop=FALSE]

for (col in colnames(coldata)) {
  coldata[[col]] <- as.factor(coldata[[col]])
}

# 3. Create DGEList, Filter & Normalize
dge <- DGEList(counts=as.matrix(cts), samples=coldata)
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge, method="TMM")

# 4. Apply Voom Transformation & Fit Linear Model
design_formula <- as.formula(args$design)
design_mat <- model.matrix(design_formula, data=coldata)

v <- voom(dge, design_mat, plot=FALSE)
fit <- lmFit(v, design_mat)

# 5. Contrast Testing & Empirical Bayes Moderation
coef_idx <- 2
if (nchar(args$contrast) > 0) {
  parts <- unlist(strsplit(args$contrast, "_vs_"))
  if (length(parts) == 2) {
    sub_parts <- unlist(strsplit(parts[1], "_"))
    factor_name <- sub_parts[1]
    num_level <- paste(sub_parts[-1], collapse="_")
    coef_name <- paste0(factor_name, num_level)
    if (coef_name %in% colnames(design_mat)) {
      coef_idx <- coef_name
    }
  }
}

fit2 <- eBayes(fit)
res_table <- topTable(fit2, number=Inf, coef=coef_idx)
res_table$Gene <- rownames(res_table)

# Map limma columns to standard BioInsight 360 backend schema
res_df <- data.frame(
  Gene = res_table$Gene,
  log2FoldChange = res_table$logFC,
  baseMean = 2^(res_table$AveExpr),
  stat = res_table$t,
  pvalue = res_table$P.Value,
  padj = res_table$adj.P.Val,
  stringsAsFactors = FALSE
)

write.csv(res_df, file=args$out, row.names=FALSE)
cat("limma-voom execution completed successfully.\n")
