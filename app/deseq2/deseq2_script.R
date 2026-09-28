#!/usr/bin/env Rscript
suppressMessages(library(argparse))
suppressMessages(library(DESeq2))

parser <- ArgumentParser()
parser$add_argument("--counts", type="character", help="Path to raw count matrix CSV")
parser$add_argument("--meta", type="character", help="Path to metadata CSV")
parser$add_argument("--out", type="character", help="Path to save results CSV")
parser$add_argument("--design", type="character", default="~ condition", help="Design formula")
parser$add_argument("--contrast", type="character", help="Contrast string e.g. condition_Treatment_vs_Control")

args <- parser$parse_args()

# Load data
cts <- read.csv(args$counts, row.names=1)
coldata <- read.csv(args$meta, row.names=1)

# Ensure sample names match
cts <- cts[, rownames(coldata)]

# Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = as.formula(args$design))

# Optional: Pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- DESeq(dds)

# Get Results
if (!is.null(args$contrast) && args$contrast != "") {
    # Parse contrast string assuming format "factor_numerator_vs_denominator"
    parts <- unlist(strsplit(args$contrast, "_vs_"))
    if(length(parts) == 2) {
        factor_name <- strsplit(parts[1], "_")[[1]][1]
        num <- strsplit(parts[1], "_")[[1]][2]
        den <- parts[2]
        res <- results(dds, contrast=c(factor_name, num, den))
    } else {
        res <- results(dds)
    }
} else {
    res <- results(dds)
}

# Write out
resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), file=args$out)
