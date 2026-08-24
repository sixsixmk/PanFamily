# Title: Gene Clustering from BLAST Results
# Description: This script takes BLASTP results, builds a network graph based on identity,
#              finds clusters (gene families), and assigns a group name to each gene.

# --- 1. Setup and Argument Parsing ---

# Install required packages if they are not already installed
if (!require("igraph", quietly = TRUE)) {
  install.packages("igraph")
}
if (!require("argparse", quietly = TRUE)) {
  install.packages("argparse")
}

library(igraph)
library(argparse)
parser <- ArgumentParser(description="Cluster genes from BLAST results and assign group names.")
parser$add_argument("-i", "--input", type="character", required=TRUE,
                    help="Path to the input BLAST results CSV file. Must have 3 columns: Gene1, Gene2, Identity.")

parser$add_argument("-o", "--output", type="character", required=TRUE,
                    help="Path for the output CSV file for the gene-to-group mapping.")

parser$add_argument("-t", "--identity_threshold", type="double", default=90.0,
                    help="Minimum identity percentage to consider a link [default: %(default)s].")

parser$add_argument("-p", "--prefix", type="character", required=TRUE,
                    help="Prefix for the generated cluster names.")

parser$add_argument("-s", "--start_number", type="integer", required=TRUE,
                    help="Starting number for the generated cluster names.")

# Parse the arguments provided by the user
args <- parser$parse_args()


# --- 2. Core Logic ---

cat("Starting gene clustering process...\n")

# Read the BLASTP results table
cat("--> Reading BLAST results from:", args$input, "\n")
# Modified to handle tab-separated BLAST output format
blastp_results <- read.table(args$input, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

blastp_results <- blastp_results[, 1:3]
colnames(blastp_results) <- c("Gene1", "Gene2", "Identity") # Rename columns for consistency

# Filter gene pairs with identity >= threshold
cat("--> Filtering pairs with identity >=", args$identity_threshold, "%\n")
filtered_data <- subset(blastp_results, Identity >= args$identity_threshold)

if (nrow(filtered_data) == 0) {
  cat("Warning: No gene pairs passed the identity threshold. Output file will be empty.\n")
  # Create an empty file and exit
  file.create(args$output)
  quit(save = "no", status = 0)
}

# Build a gene network graph
cat("--> Building gene network...\n")
edges <- filtered_data[, c("Gene1", "Gene2")]
g <- graph_from_data_frame(edges, directed = FALSE)

# Find connected components (gene clusters)
cat("--> Identifying gene clusters...\n")
clusters <- clusters(g)$membership
gene_groups <- split(names(clusters), clusters)
cat("--> Found", length(gene_groups), "gene clusters.\n")

# Create new column names for the clusters
column_names <- paste0(args$prefix, args$start_number:(args$start_number - 1 + length(gene_groups)))

# Convert the list of groups to a two-column data frame (GeneID, GroupName)
gene_to_group_mapping <- data.frame(
  GeneID = unlist(gene_groups, use.names = FALSE),
  Group = rep(column_names, sapply(gene_groups, length))
)


# --- 3. Save Results ---
cat("\n=== Gene Clustering Results Summary ===\n")
cat("Total genes clustered:", nrow(gene_to_group_mapping), "\n")
cat("Number of gene families:", length(gene_groups), "\n")
cat("Identity threshold used:", args$identity_threshold, "%\n\n")
cluster_sizes <- table(gene_to_group_mapping$Group)
cluster_summary <- data.frame(
  GeneFamily = names(cluster_sizes),
  Size = as.numeric(cluster_sizes)
)
cluster_summary <- cluster_summary[order(cluster_summary$Size, decreasing = TRUE), ]
cat("=== Gene Family Size Distribution ===\n")
print(cluster_summary, row.names = FALSE)
cat("\n=== Sample Gene-to-Family Mappings ===\n")
sample_size <- min(10, nrow(gene_to_group_mapping))
print(head(gene_to_group_mapping, sample_size), row.names = FALSE)

if (nrow(gene_to_group_mapping) > sample_size) {
  cat("... and", nrow(gene_to_group_mapping) - sample_size, "more entries\n")
}
write.table(gene_to_group_mapping, args$output, 
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n--> Process complete. Full mapping result saved to:", args$output, "\n")
summary_output <- gsub("\\.(csv|tsv)$", "_summary.tsv", args$output)
write.table(cluster_summary, summary_output, 
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("--> Cluster summary saved to:", summary_output, "\n")

