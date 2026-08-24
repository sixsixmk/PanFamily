#!/usr/bin/env Rscript

suppressMessages(library(optparse))

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, 
              help = "Input PAV matrix TSV file (e.g., final_pav_matrix.tsv)", metavar = "FILE"),
  make_option(c("-o", "--output"), type = "character", default = "gene_categories.tsv", 
              help = "Output file for gene classification [default: %default]", metavar = "FILE"),
  make_option(c("--core_threshold"), type = "double", default = 0.95, 
              help = "Threshold for core genes (default: 95%%)", metavar = "FLOAT")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$input)) {
  stop("Error: Please provide an input PAV matrix file using -i or --input. Use --help for more details.")
}

pav_matrix <- read.table(opt$input, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

gene_presence_counts <- rowSums(pav_matrix)
num_strains <- ncol(pav_matrix)
core_threshold <- opt$core_threshold * num_strains   
specific_threshold <- 1                              

gene_categories <- data.frame(
  Gene = rownames(pav_matrix),
  PresenceCount = gene_presence_counts,
  Category = ifelse(
    gene_presence_counts >= core_threshold, "Core Gene",
    ifelse(gene_presence_counts <= specific_threshold, "Private Gene", "Dispenable Gene")
  )
)

write.table(gene_categories, opt$output, sep = "\t", row.names = FALSE, quote = FALSE)


cat("Classification results:\n")
print(table(gene_categories$Category))

cat("\nOutput saved to:", opt$output, "\n")
