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
parser$add_argument("--numerator", type="character", default="", help="Treatment level for the contrast")
parser$add_argument("--denominator", type="character", default="", help="Control level for the contrast")

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
keep <- filterByExpr(dge, group=coldata[["condition"]])
dge <- dge[keep, , keep.lib.sizes=FALSE]

# 4. TMM Normalization & Dispersion Estimation
dge <- calcNormFactors(dge, method="TMM")
design_formula <- as.formula(args$design)
design_mat <- model.matrix(design_formula, data=coldata)

dge <- estimateDisp(dge, design_mat)
fit <- glmQLFit(dge, design_mat)

# 5. Contrast Testing
if (nchar(args$numerator) > 0 && nchar(args$denominator) > 0) {
  factor_name <- "condition"
  if (!(factor_name %in% colnames(coldata))) {
    stop("The --numerator/--denominator contrast needs a metadata column named 'condition'.")
  }
  levels_seen <- levels(coldata[[factor_name]])
  if (!(args$numerator %in% levels_seen) || !(args$denominator %in% levels_seen)) {
    stop(sprintf(
      "Contrast levels not found. Numerator '%s', denominator '%s', available: %s",
      args$numerator, args$denominator, paste(levels_seen, collapse=", ")
    ))
  }
  # Build the contrast from R's actual design contrasts instead of assuming
  # how model.matrix sanitises level names. This supports group labels with
  # spaces, underscores, punctuation, or non-ASCII characters.
  condition_contrasts <- contrasts(coldata[[factor_name]])
  if (is.null(condition_contrasts) || ncol(design_mat) != ncol(condition_contrasts) + 1) {
    stop("The numerator/denominator contrast requires a single-factor design such as '~ condition'.")
  }
  if (!(args$numerator %in% rownames(condition_contrasts)) ||
      !(args$denominator %in% rownames(condition_contrasts))) {
    stop(sprintf(
      "Contrast levels do not match the design: numerator '%s', denominator '%s', available: %s",
      args$numerator, args$denominator, paste(rownames(condition_contrasts), collapse=", ")
    ))
  }
  contrast_values <- condition_contrasts[args$numerator, , drop=TRUE] -
                     condition_contrasts[args$denominator, , drop=TRUE]
  cv <- setNames(rep(0, ncol(design_mat)), colnames(design_mat))
  cv[2:ncol(design_mat)] <- as.numeric(contrast_values)
  qlf <- glmQLFTest(fit, contrast=cv)
} else if (nchar(args$contrast) > 0) {
  parts <- unlist(strsplit(args$contrast, "_vs_"))
  if (length(parts) != 2) {
    stop(sprintf("Unrecognised contrast '%s' (expected factor_level_vs_level).", args$contrast))
  }
  sub_parts <- unlist(strsplit(parts[1], "_"))
  factor_name <- sub_parts[1]
  num_level <- paste(sub_parts[-1], collapse="_")
  den_level <- parts[2]
  num_coef <- paste0(factor_name, num_level)
  den_coef <- paste0(factor_name, den_level)

  cv <- setNames(rep(0, ncol(design_mat)), colnames(design_mat))
  if (num_coef %in% names(cv)) cv[num_coef] <- cv[num_coef] + 1
  if (den_coef %in% names(cv)) cv[den_coef] <- cv[den_coef] - 1
  if (!(num_coef %in% names(cv)) || !(den_coef %in% names(cv))) {
    stop(sprintf("Contrast '%s' does not match the design matrix columns: %s",
                 args$contrast, paste(colnames(design_mat), collapse=", ")))
  }
  qlf <- glmQLFTest(fit, contrast=cv)
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
  logCPM = res_table$logCPM,
  stat = res_table$F,
  pvalue = res_table$PValue,
  padj = res_table$FDR,
  stringsAsFactors = FALSE
)

write.csv(res_df, file=args$out, row.names=FALSE)
cat("edgeR execution completed successfully.\n")
