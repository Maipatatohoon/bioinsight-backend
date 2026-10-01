#!/usr/bin/env Rscript
# ==============================================================================
# limma-voom Execution Script for BioInsight 360 Backend
# ==============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(argparse)
})

# 1. Setup Argument Parser
parser <- ArgumentParser(description="Run limma-voom analysis on count data.")
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

metadata[[factor_name]] <- factor(metadata[[factor_name]])

# 4. Filter and Normalize using DGEList + TMM
dge <- DGEList(counts=counts)
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

# 5. Design Matrix
design <- model.matrix(form, data=metadata)

# 6. Apply Voom Transformation
v <- voom(dge, design, plot=FALSE)

# 7. Fit Linear Model
fit <- lmFit(v, design)

# 8. Setup Contrast Matrix & eBayes
coef_name <- paste0(factor_name, level_numerator)
if (!coef_name %in% colnames(design)) {
    design <- model.matrix(~ 0 + metadata[[factor_name]])
    colnames(design) <- levels(metadata[[factor_name]])
    v <- voom(dge, design, plot=FALSE)
    fit <- lmFit(v, design)
    
    contrast_matrix <- makeContrasts(
        contrasts = paste(level_numerator, "-", level_denominator),
        levels = design
    )
    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2)
} else {
    fit2 <- eBayes(fit)
}

# 9. Format Results
res <- topTable(fit2, number=Inf, coef=if (coef_name %in% colnames(design)) coef_name else 1)
res$Gene <- rownames(res)
res <- res[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]

# Rename columns for unified backend schema (DESeq2 standard keys)
colnames(res) <- c("Gene", "log2FoldChange", "baseMean", "stat", "pvalue", "padj")

# 10. Save Output
write.csv(res, args$output, row.names=FALSE, quote=FALSE)
cat("limma-voom analysis completed successfully.\n")
