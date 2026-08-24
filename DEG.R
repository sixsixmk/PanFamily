if (!require("DESeq2"))
    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    BiocManager::install("DESeq2")

library(DESeq2)

deseq2 <- function(subcountMatrix, sampleInfo, pvalue_cutoff, log2fc_cutoff) {
    condition <- sampleInfo$condition
    dds <- DESeqDataSetFromMatrix(countData = subcountMatrix, colData = sampleInfo, design = ~condition)
    dds$condition <- relevel(dds$condition, ref="ctrl")
    keep <- rowSums(counts(dds)) >= 1.5 * ncol(subcountMatrix)
    dds <- dds[keep, ]
    dds <- DESeq(dds, quiet = FALSE)
    dds <- na.omit(dds)
    res <- results(dds, contrast = c("condition", "trait", "ctrl"))
    return(res)
}

DEGanalysis <- function(countMatrix, ctrl, trait, pvalue_cutoff, log2fc_cutoff) {
    sampleInfo <- data.frame(
        condition = c(rep("ctrl", length(ctrl)), rep("trait", length(trait))),
        row.names = c(ctrl, trait)
    )
    sampleInfo$condition <- factor(sampleInfo$condition)
    res <- deseq2(countMatrix, sampleInfo, pvalue_cutoff, log2fc_cutoff)
    res <- na.omit(res)
    filter_res <- res[res$padj < pvalue_cutoff, ]
    
   
    up <- filter_res[filter_res$log2FoldChange > log2fc_cutoff, ]    
    down <- filter_res[filter_res$log2FoldChange < -log2fc_cutoff, ]  
    mid <- filter_res[abs(filter_res$log2FoldChange) <= log2fc_cutoff, ] 
    
    return(list(upregulated = up, downregulated = down, intermediate = mid, all = res))
}


parse_group_file <- function(file) {
    # Read file: Control Treatment
    df <- read.table(file, header = TRUE, sep = "\t", comment.char = "", stringsAsFactors = FALSE)
    
    if (ncol(df) != 2) {
        stop("The group file must contain exactly two columns: Control and Treatment")
    }
    
    colnames(df) <- c("Control", "Treatment")

    ctrl_samples <- trimws(df$Control)
    ctrl_samples <- ctrl_samples[ctrl_samples != "" & !is.na(ctrl_samples)]
    
    treat_samples <- trimws(df$Treatment)
    treat_samples <- treat_samples[treat_samples != "" & !is.na(treat_samples)]
    
    return(list(control = ctrl_samples, treatment = treat_samples))
}

args <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    return(args)
}

main <- function() {
    args <- args()
    
    if (length(args) < 3) {
        stop("Usage: Rscript DEG.R <count_matrix> <group_file> <output_dir> [pvalue_cutoff] [log2fc_cutoff]")
    }
    
    countMatrix_path <- args[1]
    group_file_path <- args[2]
    output_dir <- args[3]
    pvalue_cutoff <- ifelse(length(args) >= 4, as.numeric(args[4]), 0.05)
    log2fc_cutoff <- ifelse(length(args) >= 5, as.numeric(args[5]), 1)
    
    cat("========================================\n")
    cat("DESeq2 Differential Expression Analysis\n")
    cat("========================================\n")
    cat(paste("Count matrix:", countMatrix_path, "\n"))
    cat(paste("Group file:", group_file_path, "\n"))
    cat(paste("Output directory:", output_dir, "\n"))
    cat(paste("P-value cutoff:", pvalue_cutoff, "\n"))
    cat(paste("Log2FC cutoff:", log2fc_cutoff, "\n"))
    cat("========================================\n\n")
    
    # Create output directory
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    
    # Read group file
    cat("Reading group file...\n")
    groups <- parse_group_file(group_file_path)
    
    ctrl_samples <- groups$control
    treat_samples <- groups$treatment
    
    cat(paste("Number of control samples:", length(ctrl_samples), "\n"))
    cat(paste("  Samples:", paste(ctrl_samples, collapse=", "), "\n"))
    cat(paste("Number of treatment samples:", length(treat_samples), "\n"))
    cat(paste("  Samples:", paste(treat_samples, collapse=", "), "\n\n"))
    
    # Read count matrix
    cat("Reading count matrix...\n")
    countMatrix <- read.csv(countMatrix_path, header = TRUE, sep = ",", row.names = 1, check.names = FALSE)
    
    cat(paste("Matrix dimensions:", nrow(countMatrix), "genes ×", ncol(countMatrix), "samples\n\n"))
    
    # Check if samples exist in matrix
    all_samples <- c(ctrl_samples, treat_samples)
    missing_samples <- setdiff(all_samples, colnames(countMatrix))
    if (length(missing_samples) > 0) {
        cat("Warning: The following samples do not exist in the count matrix:\n")
        cat(paste("  -", missing_samples, "\n"))
        cat("These samples will be ignored\n\n")
    }
    
    # Filter existing samples
    ctrl_samples <- intersect(ctrl_samples, colnames(countMatrix))
    treat_samples <- intersect(treat_samples, colnames(countMatrix))
    
    if (length(ctrl_samples) < 2) {
        stop("At least 2 control samples are required")
    }
    
    if (length(treat_samples) < 2) {
        stop("At least 2 treatment samples are required")
    }
    
    # Extract submatrix
    cat("Extracting submatrix...\n")
    subCountMatrix <- countMatrix[, c(ctrl_samples, treat_samples)]
    
    cat(paste("  Control group:", length(ctrl_samples), "samples\n"))
    cat(paste("  Treatment group:", length(treat_samples), "samples\n\n"))
    
    # Run DESeq2
    cat("Running DESeq2 differential analysis...\n")
    res <- DEGanalysis(subCountMatrix, ctrl_samples, treat_samples, pvalue_cutoff, log2fc_cutoff)
    
    # Output files
    all_file <- file.path(output_dir, "DEG_all.csv")
    up_file <- file.path(output_dir, "DEG_up.csv")
    down_file <- file.path(output_dir, "DEG_down.csv")
    mid_file <- file.path(output_dir, "DEG_intermediate.csv")
    
    write.csv(res$all, file = all_file, row.names = TRUE)
    write.csv(res$upregulated, file = up_file, row.names = TRUE)
    write.csv(res$downregulated, file = down_file, row.names = TRUE)
    write.csv(res$intermediate, file = mid_file, row.names = TRUE)
    
    cat("\n========================================\n")
    cat("Differential Expression Analysis Results:\n")
    cat(paste("  Total genes:", nrow(res$all[!is.na(res$all$padj), ]), "\n"))
    cat(paste("  Significant DEGs (padj <", pvalue_cutoff, "):", 
              nrow(res$all[!is.na(res$all$padj) & res$all$padj < pvalue_cutoff, ]), "\n"))
    cat(paste("  Upregulated genes (log2FC >", log2fc_cutoff, "):", nrow(res$upregulated), "\n"))
    cat(paste("  Downregulated genes (log2FC <", -log2fc_cutoff, "):", nrow(res$downregulated), "\n"))
    cat(paste("  Intermediate expression (|log2FC| ≤", log2fc_cutoff, "):", nrow(res$intermediate), "\n"))
    cat("\nOutput files:\n")
    cat(paste("  - All results:", all_file, "\n"))
    cat(paste("  - Upregulated genes:", up_file, "\n"))
    cat(paste("  - Downregulated genes:", down_file, "\n"))
    cat(paste("  - Intermediate expression:", mid_file, "\n"))
    cat("========================================\n")
    cat("Analysis completed!\n")
    cat("========================================\n")
}

main()
