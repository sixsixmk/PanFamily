#!/usr/bin/env Rscript

library(dplyr)
library(argparse)

parser <- ArgumentParser(description = "Merge multiple PAV matrix files")
parser$add_argument("-dir", "--directory", type = "character", required = TRUE,
                    help = "Directory path containing PAV matrix files")
parser$add_argument("-o", "--output", type = "character", default = "final_merged_pav_matrix.tsv",
                    help = "Path to the final merged PAV matrix output file (default: final_merged_pav_matrix.tsv)")
args <- parser$parse_args()
directory <- args$directory
pav_files <- list.files(directory, pattern = "_pav\\.(csv|tsv)$", full.names = TRUE)
pav_files <- pav_files[!grepl("merge_pav\\.R|final_merged_pav_matrix\\.(csv|tsv)", pav_files)]
final_pav_matrix <- NULL
for (file in pav_files) {
  cat("Processing file:", file, "\n")
  pav_matrix <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  if (is.null(final_pav_matrix)) {
    rownames(pav_matrix) <- pav_matrix[, 1]
    final_pav_matrix <- pav_matrix[, -1, drop = FALSE]
  } else {
    rownames(pav_matrix) <- pav_matrix[, 1]
    pav_matrix <- pav_matrix[, -1, drop = FALSE]
    all_genes <- union(rownames(final_pav_matrix), rownames(pav_matrix))
    temp_final <- matrix(0, nrow = length(all_genes), ncol = ncol(final_pav_matrix))
    rownames(temp_final) <- all_genes
    colnames(temp_final) <- colnames(final_pav_matrix)
    temp_final[rownames(final_pav_matrix), ] <- as.matrix(final_pav_matrix)
    
    temp_current <- matrix(0, nrow = length(all_genes), ncol = ncol(pav_matrix))
    rownames(temp_current) <- all_genes
    colnames(temp_current) <- colnames(pav_matrix)
    temp_current[rownames(pav_matrix), ] <- as.matrix(pav_matrix)
    final_pav_matrix <- cbind(temp_final, temp_current)
  }
}
final_pav_matrix[is.na(final_pav_matrix)] <- 0
write.table(final_pav_matrix, args$output, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

cat("Merging completed. Output file:", args$output, "\n")
