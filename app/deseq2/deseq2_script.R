#!/usr/bin/env Rscript
# ==============================================================================
# DESeq2 Execution Script for BioInsight 360 Backend
# ==============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(argparse)
})

# 1. Parse Arguments
parser <- ArgumentParser(description="Run DESeq2 differential expression analysis.")
parser$add_argument("--counts", type="character", required=TRUE, help="Path to raw count matrix CSV")
parser$add_argument("--meta", type="character", required=TRUE, help="Path to metadata CSV")
parser$add_argument("--out", type="character", required=TRUE, help="Path to save output results CSV")
parser$add_argument("--design", type="character", default="~ condition", help="Design formula (e.g., ~ condition)")
parser$add_argument("--contrast", type="character", default="", help="Contrast string (e.g., condition_Treated_vs_Control)")

args <- parser$parse_args()

# 2. Load Input Files
cts <- read.csv(args$counts, row.names=1, check.names=FALSE)
coldata <- read.csv(args$meta, row.names=1, check.names=FALSE)

# Ensure matching sample order
common_samples <- intersect(colnames(cts), rownames(coldata))
if (length(common_samples) < 2) {
  stop("Fewer than 2 overlapping samples between count matrix and metadata.")
}
cts <- cts[, common_samples, drop=FALSE]
coldata <- coldata[common_samples, , drop=FALSE]

# Convert all factor columns in coldata to factors
for (col in colnames(coldata)) {
  coldata[[col]] <- as.factor(coldata[[col]])
}

# 3. Construct DESeqDataSet
design_formula <- as.formula(args$design)
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(cts)),
  colData = coldata,
  design = design_formula
)

# 4. Low-count Pre-filtering (at least 10 reads across samples)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# 5. Execute DESeq Pipeline
dds <- DESeq(dds, quiet=TRUE)

# 6. Extract Results
if (nchar(args$contrast) > 0) {
  # Parse contrast string format: factor_num_vs_den
  parts <- unlist(strsplit(args$contrast, "_vs_"))
  if (length(parts) == 2) {
    sub_parts <- unlist(strsplit(parts[1], "_"))
    factor_name <- sub_parts[1]
    num_level <- paste(sub_parts[-1], collapse="_")
    den_level <- parts[2]
    res <- results(dds, contrast=c(factor_name, num_level, den_level))
  } else {
    res <- results(dds)
  }
} else {
  res <- results(dds)
}

# 7. Format Output DataFrame
res_df <- as.data.frame(res)
res_df$Gene <- rownames(res_df)

# Standardize output column names
colnames(res_df)[colnames(res_df) == "log2FoldChange"] <- "log2FoldChange"
colnames(res_df)[colnames(res_df) == "pvalue"] <- "pvalue"
colnames(res_df)[colnames(res_df) == "padj"] <- "padj"
colnames(res_df)[colnames(res_df) == "baseMean"] <- "baseMean"
colnames(res_df)[colnames(res_df) == "stat"] <- "stat"

# Reorder columns cleanly
cols_order <- c("Gene", "log2FoldChange", "baseMean", "stat", "pvalue", "padj")
existing_cols <- intersect(cols_order, colnames(res_df))
res_df <- res_df[, c(existing_cols, setdiff(colnames(res_df), existing_cols))]

write.csv(res_df, file=args$out, row.names=FALSE)
cat("DESeq2 execution completed successfully.\n")
