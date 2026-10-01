#!/usr/bin/env Rscript
# ==============================================================================
# edgeR Execution Script for BioInsight 360 Backend
# ==============================================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(argparse)
})

# 1. Setup Argument Parser
parser <- ArgumentParser(description="Run edgeR analysis on count data.")
parser$add_argument("--counts", type="character", required=TRUE, help="Path to input counts CSV")
parser$add_argument("--metadata", type="character", required=TRUE, help="Path to input metadata CSV")
parser$add_argument("--formula", type="character", required=TRUE, help="Design formula (e.g., ~ condition)")
parser$add_argument("--contrast", type="character", required=TRUE, help="Contrast string (e.g., condition_Treated_vs_Control)")
parser$add_argument("--output", type="character", required=TRUE, help="Path to save output results CSV")
args <- parser$parse_args()

# 2. Load Data
counts <- read.csv(args$counts, row.names=1, check.names=FALSE)
metadata <- read.csv(args$metadata, row.names=1, check.names=FALSE)

# Ensure sample names match
common_samples <- intersect(colnames(counts), rownames(metadata))
counts <- counts[, common_samples, drop=FALSE]
metadata <- metadata[common_samples, , drop=FALSE]

# 3. Parse Formula and Contrast
form <- as.formula(args$formula)
contrast_parts <- unlist(strsplit(args$contrast, "_"))
if (length(contrast_parts) != 4 || contrast_parts[3] != "vs") {
    stop(paste("Invalid contrast format:", args$contrast, ". Expected format: Factor_LevelA_vs_LevelB"))
}

factor_name <- contrast_parts[1]
level_numerator <- contrast_parts[2]
level_denominator <- contrast_parts[4]

if (!factor_name %in% colnames(metadata)) {
    stop(paste("Factor", factor_name, "not found in metadata columns."))
}

# Ensure the factor is treated as a factor variable
metadata[[factor_name]] <- factor(metadata[[factor_name]])

# 4. Construct DGEList and Filter
y <- DGEList(counts=counts, samples=metadata)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

# 5. Normalization
y <- calcNormFactors(y)

# 6. Design Matrix
design <- model.matrix(form, data=metadata)

# 7. Dispersion Estimation
y <- estimateDisp(y, design)

# 8. Exact Test (or GLM if complex design)
fit <- glmQLFit(y, design)

# Identify the coefficient for the contrast
coef_name <- paste0(factor_name, level_numerator)
if (!coef_name %in% colnames(design)) {
    # If standard naming fails, create a custom contrast matrix
    design <- model.matrix(~ 0 + metadata[[factor_name]])
    colnames(design) <- levels(metadata[[factor_name]])
    
    my_contrast <- makeContrasts(
        contrasts = paste(level_numerator, "-", level_denominator),
        levels = design
    )
    qlf <- glmQLFTest(fit, contrast=my_contrast)
} else {
    qlf <- glmQLFTest(fit, coef=coef_name)
}

# 9. Format Results
res <- topTags(qlf, n=Inf)$table
res$Gene <- rownames(res)
res <- res[, c("Gene", "logFC", "logCPM", "F", "PValue", "FDR")]

# Rename columns to match standard backend output expectations (DESeq2 format style for UI compatibility)
colnames(res) <- c("Gene", "log2FoldChange", "baseMean", "stat", "pvalue", "padj")

# 10. Save Output
write.csv(res, args$output, row.names=FALSE, quote=FALSE)
cat("edgeR analysis completed successfully.\n")
