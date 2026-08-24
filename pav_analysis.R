#!/usr/bin/env Rscript

# PAV (Presence/Absence Variation) Analysis Script

if (!require("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!require("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}
if (!require("argparse", quietly = TRUE)) {
  install.packages("argparse")
}

library(dplyr)
library(tidyr)
library(argparse)

parser <- ArgumentParser(description = "Analyze gene presence/absence variation (PAV) based on BLAST results")

parser$add_argument("-i", "--input", type = "character", required = TRUE,
                    help = "Input BLAST result file (tab-delimited format)")

parser$add_argument("-o", "--output", type = "character", required = TRUE,
                    help = "Output PAV matrix file (.tsv format)")


parser$add_argument("-s", "--strain", type = "character", required = TRUE,
                    help = "Strain/cultivar identifier for this analysis")

parser$add_argument("-t", "--threshold", type = "double", default = 90.0,
                    help = "Identity threshold for presence determination [default: %(default)s]")

parser$add_argument("--evalue", type = "double", default = 1e-5,
                    help = "E-value threshold for filtering BLAST hits [default: %(default)s]")

parser$add_argument("--length", type = "integer", default = 50,
                    help = "Minimum alignment length threshold [default: %(default)s]")


parser$add_argument("-q", "--query_prefix", type = "character", default = "",
                    help = "Prefix to filter query sequences (e.g., 'TaTPS' for TaTPS genes)")

parser$add_argument("--query_list", type = "character", default = "",
                    help = "File containing list of query genes to analyze (one per line)")


parser$add_argument("-v", "--verbose", action = "store_true", default = FALSE,
                    help = "Enable verbose output")

parser$add_argument("--summary", action = "store_true", default = FALSE,
                    help = "Generate summary statistics")


args <- parser$parse_args()


log_info <- function(message) {
  if (args$verbose) {
    cat("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] INFO: ", message, "\n", sep = "")
  }
}

log_info("Starting PAV analysis...")
log_info(paste("Input file:", args$input))
log_info(paste("Strain:", args$strain))
log_info(paste("Identity threshold:", args$threshold, "%"))

# Read BLAST results
log_info("Reading BLAST results...")
blast_result <- read.table(args$input, header = FALSE, sep = "\t", stringsAsFactors = FALSE)


if (ncol(blast_result) < 12) {
  stop("BLAST file should have at least 12 columns (standard -outfmt 6 format)")
}


colnames(blast_result) <- c("query", "subject", "identity", "length", "mismatch", "gapopen", 
                           "qstart", "qend", "sstart", "send", "evalue", "bitscore")

log_info(paste("Read", nrow(blast_result), "BLAST hits"))


log_info("Applying filters...")
original_count <- nrow(blast_result)


blast_result <- blast_result %>% filter(evalue <= args$evalue)
log_info(paste("After e-value filter (<=", args$evalue, "):", nrow(blast_result), "hits"))


blast_result <- blast_result %>% filter(length >= args$length)
log_info(paste("After length filter (>=", args$length, "):", nrow(blast_result), "hits"))

if (args$query_prefix != "") {
  blast_result <- blast_result %>% filter(grepl(paste0("^", args$query_prefix), query))
  log_info(paste("After query prefix filter (", args$query_prefix, "):", nrow(blast_result), "hits"))
}

if (args$query_list != "") {
  if (file.exists(args$query_list)) {
    query_genes <- read.table(args$query_list, header = FALSE, stringsAsFactors = FALSE)$V1
    blast_result <- blast_result %>% filter(query %in% query_genes)
    log_info(paste("After query list filter:", nrow(blast_result), "hits"))
  } else {
    warning(paste("Query list file not found:", args$query_list))
  }
}

if (nrow(blast_result) == 0) {
  stop("No BLAST hits remaining after filtering!")
}

blast_result$strain <- args$strain

log_info("Processing best hits per query...")
best_hits <- blast_result %>%
  group_by(query, strain) %>%
  slice_max(identity, n = 1, with_ties = FALSE) %>% 
  ungroup()

log_info(paste("Selected", nrow(best_hits), "best hits"))
best_hits$present <- ifelse(best_hits$identity >= args$threshold, 1, 0)
all_queries <- unique(best_hits$query)
log_info(paste("Analyzing", length(all_queries), "query genes"))
log_info("Creating PAV matrix...")

strain_column <- sapply(all_queries, function(gene) {
  hit <- best_hits[best_hits$query == gene, ]
  if (nrow(hit) > 0) {
    return(hit$present[1])
  } else {
    return(0)  
  }
})


pav_matrix <- data.frame(
  Gene = all_queries,
  stringsAsFactors = FALSE
)

pav_matrix[[args$strain]] <- strain_column


if (args$summary) {
  log_info("Generating summary statistics...")
  
  total_genes <- nrow(pav_matrix)
  present_genes <- sum(pav_matrix[[args$strain]])
  absent_genes <- total_genes - present_genes
  presence_rate <- round(present_genes / total_genes * 100, 2)
  
  cat("\n=== PAV Analysis Summary ===\n")
  cat("Strain:", args$strain, "\n")
  cat("Total query genes:", total_genes, "\n")
  cat("Present genes:", present_genes, "\n")
  cat("Absent genes:", absent_genes, "\n")
  cat("Presence rate:", presence_rate, "%\n")
  
  if (args$verbose) {
    cat("\n=== Present genes ===\n")
    present_gene_names <- pav_matrix$Gene[pav_matrix[[args$strain]] == 1]
    cat(paste(present_gene_names, collapse = ", "), "\n")
    
    cat("\n=== Absent genes ===\n")
    absent_gene_names <- pav_matrix$Gene[pav_matrix[[args$strain]] == 0]
    cat(paste(absent_gene_names, collapse = ", "), "\n")
  }
  

  summary_file <- gsub("\\.(csv|tsv|txt)$", "_summary.txt", args$output)
  summary_content <- paste(
    "PAV Analysis Summary",
    "===================",
    paste("Analysis date:", Sys.time()),
    paste("Strain:", args$strain),
    paste("Input file:", args$input),
    paste("Identity threshold:", args$threshold, "%"),
    paste("E-value threshold:", args$evalue),
    paste("Length threshold:", args$length),
    "",
    "Results:",
    paste("Total query genes:", total_genes),
    paste("Present genes:", present_genes),
    paste("Absent genes:", absent_genes),
    paste("Presence rate:", presence_rate, "%"),
    "",
    "Present genes:",
    paste(present_gene_names, collapse = ", "),
    "",
    "Absent genes:",
    paste(absent_gene_names, collapse = ", "),
    sep = "\n"
  )
  
  write(summary_content, summary_file)
  log_info(paste("Summary saved to:", summary_file))
}

log_info(paste("Saving PAV matrix to:", args$output))
write.table(pav_matrix, args$output, sep = "\t", row.names = FALSE, quote = FALSE)

log_info("PAV analysis completed successfully!")

if (args$verbose) {
  cat("\n=== PAV Matrix Preview ===\n")
  print(head(pav_matrix, 10))
  if (nrow(pav_matrix) > 10) {
    cat("... and", nrow(pav_matrix) - 10, "more genes\n")
  }
}