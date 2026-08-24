#!/usr/bin/env Rscript

# ============================================================================
# PanFamily Single Cell RNA-seq Analysis - Core R Script
# ============================================================================

  library(optparse)
  library(yaml)
  library(Seurat)
  library(scMayoMap)
  library(dplyr)
  library(ggrepel)
  library(ggplot2)
  library(openxlsx)
  library(patchwork)
  library(monocle) 
  library(stringr)
  library(Matrix)
  library(tidyverse)
  library(AnnotationDbi)
  library(org.At.tair.db)
  library(ComplexHeatmap)
  library(igraph)
  library(HGNChelper)
  library(ggsci)
 
option_list <- list(
  make_option(c("-c", "--config"), type="character", default=NULL,
              help="Configuration file path", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=".",
              help="Output directory", metavar="character"),
  make_option(c("-s", "--step"), type="character", default="all",
              help="Execution step", metavar="character"),
  make_option(c("-t", "--threads"), type="integer", default=4,
              help="Number of threads", metavar="integer"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Verbose output")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate parameters
if (is.null(opt$config)) {
  stop("Configuration file must be provided (--config)")
}

# ============================================================================
# Utility Functions
# ============================================================================

# Get nested configuration value
get_nested_config <- function(config, key_path, default = NULL) {
  keys <- strsplit(key_path, "\\.")[[1]]
  value <- config
  
  for (key in keys) {
    if (key %in% names(value)) {
      value <- value[[key]]
    } else {
      return(default)
    }
  }
  
  return(value)
}

# Log functions
log_info <- function(message) {
  cat(sprintf("[INFO] %s - %s\n", Sys.time(), message))
}

log_error <- function(message) {
  cat(sprintf("[ERROR] %s - %s\n", Sys.time(), message), file = stderr())
}

log_step <- function(step, message) {
  cat(sprintf("[STEP %s] %s - %s\n", step, Sys.time(), message))
}

# Create safe output filename
safe_filename <- function(name, ext = "png") {
  clean_name <- gsub("[^A-Za-z0-9_-]", "_", name)
  return(paste0(clean_name, ".", ext))
}

# General function for saving plots
save_plot <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  tryCatch({
    ggsave(filename, plot = plot, width = width, height = height, dpi = dpi, limitsize = FALSE)
    log_info(sprintf("Plot saved: %s", filename))
  }, error = function(e) {
    log_error(sprintf("Failed to save plot %s: %s", filename, e$message))
  })
}

# ============================================================================
# Interactive Helper Functions
# ============================================================================

# Check if in interactive mode
is_interactive_mode <- function(config, step_name) {
  execution_mode <- get_nested_config(config, "execution.mode", default = "interactive")
  if (execution_mode == "batch") return(FALSE)
  step_interactive <- get_nested_config(config, paste0("execution.interactive_settings.steps.", step_name), default = FALSE)
  return(step_interactive)
}

# Simplified user input function with timeout
get_user_input_with_timeout <- function(prompt, timeout_seconds = 120) {
  # Display prompt
  cat(prompt)
  flush.console()  # Ensure prompt is displayed immediately
  
  # Use the most reliable method to read user input
  tryCatch({
    # Use scan to read user input, reliable in most environments
    result <- scan(file = "stdin", what = character(), nlines = 1, quiet = TRUE)
    
    if (length(result) == 0) {
      return("")
    } else {
      return(trimws(result[1]))
    }
  }, error = function(e) {
    # Backup method: use readline
    cat("\n[Backup input method] Please enter: ")
    flush.console()
    tryCatch({
      result <- readline()
      return(trimws(result))
    }, error = function(e2) {
      # Last resort: return empty string, use default values
      cat("\n[Input read failed, using default values]\n")
      return("")
    })
  })
}

# Checkpoint detection function
detect_checkpoint <- function(output_dir) {
  checkpoints <- list(
    list(step = "branch_completed", file = "monocle_cds_after_branch.rds", name = "Monocle object after branch analysis completion"),
    list(step = "branch_analysis", file = "monocle_cds.rds", name = "Pseudotime Monocle object"),
    list(step = "pseudotime", file = "annotated_seurat_object.rds", name = "Annotated Seurat object"),
    list(step = "annotation", file = "clustered_seurat_object.rds", name = "Clustered Seurat object")
  )
  
  for (checkpoint in checkpoints) {
    file_path <- file.path(output_dir, checkpoint$file)
    if (file.exists(file_path)) {
      log_info(sprintf("Checkpoint detected: %s (%s)", checkpoint$name, checkpoint$file))
      return(list(step = checkpoint$step, file = file_path, object_file = checkpoint$file))
    }
  }
  
  log_info("No valid checkpoint detected, will start analysis from beginning")
  return(list(step = "start", file = NULL, object_file = NULL))
}

# Interactive plot sizing with timeout mechanism
interactive_plot_sizing <- function(plot_func, plot_name, default_width, default_height, output_dir, config) {
  # Check if in interactive mode
  execution_mode <- get_nested_config(config, "execution.mode", default = "interactive")
  if (execution_mode == "batch") {
    # Batch mode, use default dimensions directly
    plot <- plot_func()
    save_plot(plot, file.path(output_dir, paste0(plot_name, ".png")), 
              width = default_width, height = default_height)
    return(plot)
  }
  
  # Interactive mode
  current_width <- default_width
  current_height <- default_height
  
  while (TRUE) {
    # Generate plot with current dimensions
    plot <- plot_func()
    temp_file <- file.path(output_dir, paste0(plot_name, "_temp.png"))
    save_plot(plot, temp_file, width = current_width, height = current_height)
    
    # Display interactive information
    cat("\n=== Interactive Plot Sizing ===\n")
    cat("Plot name:", plot_name, "\n")
    cat("Current dimensions:", current_width, "x", current_height, "\n")
    cat("Temporary file:", temp_file, "\n")
    cat("\nPlease check the plot, are you satisfied with current dimensions?\n")
    cat("Y - Satisfied, save and continue\n")
    cat("N - Not satisfied, adjust dimensions\n")
    
    # User input with timeout (120 seconds)
    choice <- get_user_input_with_timeout("Select (Y/N, will use default dimensions if no response in 120 seconds): ", 120)
    
    if (choice == "TIMEOUT" || toupper(choice) == "Y" || choice == "") {
      # Timeout or user satisfied, save final plot
      final_file <- file.path(output_dir, paste0(plot_name, ".png"))
      save_plot(plot, final_file, width = current_width, height = current_height)
      if (choice == "TIMEOUT") {
        log_info(sprintf("Interactive timeout, plot saved with default dimensions: %s (dimensions: %dx%d)", final_file, current_width, current_height))
      } else {
        log_info(sprintf("Plot saved: %s (dimensions: %dx%d)", final_file, current_width, current_height))
      }
      
      # Delete temporary file
      if (file.exists(temp_file)) {
        file.remove(temp_file)
      }
      
      return(plot)
    } else if (toupper(choice) == "N") {
      # User not satisfied, adjust dimensions
      cat("\n=== Adjust Plot Dimensions ===\n")
      cat("Current width:", current_width, "\n")
      cat("Current height:", current_height, "\n")
      
      # Get new width (with timeout)
      input_width <- get_user_input_with_timeout(sprintf("Enter new width [current:%d, will keep current value if no response in 120 seconds]: ", current_width), 120)
      if (input_width != "TIMEOUT" && input_width != "") {
        new_width <- as.numeric(input_width)
        if (!is.na(new_width) && new_width > 0) {
          current_width <- new_width
        }
      }
      
      # Get new height (with timeout)
      input_height <- get_user_input_with_timeout(sprintf("Enter new height [current:%d, will keep current value if no response in 120 seconds]: ", current_height), 120)
      if (input_height != "TIMEOUT" && input_height != "") {
        new_height <- as.numeric(input_height)
        if (!is.na(new_height) && new_height > 0) {
          current_height <- new_height
        }
      }
      
      cat("New dimensions:", current_width, "x", current_height, "\n")
      cat("Regenerating plot...\n")
      
    } else {
      cat("Invalid choice, please enter Y or N\n")
    }
  }
}

# Generate data statistics summary
generate_qc_summary <- function(seurat_obj, output_dir) {
  tryCatch({
    stats <- data.frame(
      metric = c("Cell count", "Gene count", "nFeature_RNA_min", "nFeature_RNA_max", "nFeature_RNA_median",
                "nCount_RNA_min", "nCount_RNA_max", "nCount_RNA_median",
                "percent.mt_min", "percent.mt_max", "percent.mt_median"),
      value = c(
        ncol(seurat_obj), nrow(seurat_obj),
        min(seurat_obj$nFeature_RNA), max(seurat_obj$nFeature_RNA), median(seurat_obj$nFeature_RNA),
        min(seurat_obj$nCount_RNA), max(seurat_obj$nCount_RNA), median(seurat_obj$nCount_RNA),
        round(min(seurat_obj$percent.mt), 2), round(max(seurat_obj$percent.mt), 2), round(median(seurat_obj$percent.mt), 2)
      )
    )
    
    # Save statistics file
    stats_file <- file.path(output_dir, "qc_data_stats.txt")
    writeLines(c(
      sprintf("Total cells: %d", ncol(seurat_obj)),
      sprintf("Total genes: %d", nrow(seurat_obj)),
      sprintf("nFeature_RNA: min=%d, max=%d, median=%.0f", 
              min(seurat_obj$nFeature_RNA), max(seurat_obj$nFeature_RNA), median(seurat_obj$nFeature_RNA)),
      sprintf("nCount_RNA: min=%d, max=%d, median=%.0f",
              min(seurat_obj$nCount_RNA), max(seurat_obj$nCount_RNA), median(seurat_obj$nCount_RNA)),
      sprintf("percent.mt: min=%.2f%%, max=%.2f%%, median=%.2f%%",
              min(seurat_obj$percent.mt), max(seurat_obj$percent.mt), median(seurat_obj$percent.mt))
    ), stats_file)
    
    log_info(sprintf("Data statistics saved: %s", stats_file))
    return(stats_file)
  }, error = function(e) {
    log_error(sprintf("Failed to generate data statistics: %s", e$message))
    return(NULL)
  })
}

# Get default QC parameters
get_default_qc_params <- function() {
  return(list(
    min_features_per_cell = 200,
    max_features_per_cell = 2500,
    max_mt_percent = 5.0,
    max_count_per_cell = 20000
  ))
}

# Ensure configuration parameters are complete, use defaults for missing values
ensure_config_complete <- function(config) {
  defaults <- get_default_qc_params()
  
  # Check and supplement QC parameters
  if (is.null(config$quality_control$min_features_per_cell)) {
    config$quality_control$min_features_per_cell <- defaults$min_features_per_cell
    log_info("Using default value: min_features_per_cell = 200")
  }
  if (is.null(config$quality_control$max_features_per_cell)) {
    config$quality_control$max_features_per_cell <- defaults$max_features_per_cell
    log_info("Using default value: max_features_per_cell = 2500")
  }
  if (is.null(config$quality_control$max_mt_percent)) {
    config$quality_control$max_mt_percent <- defaults$max_mt_percent
    log_info("Using default value: max_mt_percent = 5.0")
  }
  if (is.null(config$quality_control$max_count_per_cell)) {
    config$quality_control$max_count_per_cell <- defaults$max_count_per_cell
    log_info("Using default value: max_count_per_cell = 20000")
  }
  
  return(config)
}

# Wait for user confirmation to continue
wait_for_user_interaction <- function(step_name, image_file = NULL, config = NULL) {
  if (!is.null(config)) {
    execution_mode <- get_nested_config(config, "execution.mode", default = "interactive")
    if (execution_mode == "batch") return(TRUE)
  }
  
  # Create interaction signal file
  signal_file <- file.path(opt$output, paste0(".interaction_", step_name))
  
  # Write interaction information
  interaction_info <- list(
    step = step_name,
    timestamp = Sys.time(),
    image_file = image_file,
    status = "waiting"
  )
  
  writeLines(paste(names(interaction_info), interaction_info, sep = ": "), signal_file)
  
  # Wait for bash script to handle interaction and update signal file
  log_info(sprintf("Waiting for user interaction completion: %s", step_name))
  
  # Simple waiting mechanism - in practice, bash script handles interaction and deletes signal file
  max_wait <- 300  # Wait up to 5 minutes
  wait_time <- 0
  
  while (file.exists(signal_file) && wait_time < max_wait) {
    Sys.sleep(1)
    wait_time <- wait_time + 1
  }
  
  if (wait_time >= max_wait) {
    log_error("User interaction timeout")
    return(FALSE)
  }
  
  return(TRUE)
}

# ============================================================================
# Configuration Loading and Validation
# ============================================================================

load_config <- function(config_file) {
  if (!file.exists(config_file)) {
    stop(sprintf("Configuration file does not exist: %s", config_file))
  }
  
  config <- yaml.load_file(config_file)
  log_info("Configuration file loaded successfully")
  return(config)
}

# ============================================================================
# CellRanger Integration Module
# ============================================================================

# Run cellranger mkref to create reference genome
run_cellranger_mkref <- function(config, output_dir) {
  log_step("0.1", "CellRanger mkref - Create reference genome")
  
  # Get configuration parameters
  cellranger_path <- get_nested_config(config, "cellranger.cellranger_path", default = "cellranger")
  genome_name <- get_nested_config(config, "cellranger.mkref.genome_name", default = "custom_genome")
  fasta_file <- get_nested_config(config, "cellranger.mkref.fasta_file", default = NULL)
  gtf_file <- get_nested_config(config, "cellranger.mkref.gtf_file", default = NULL)
  
  if (is.null(fasta_file) || is.null(gtf_file)) {
    log_error("Missing required parameters: fasta_file or gtf_file")
    stop("CellRanger mkref requires fasta_file and gtf_file")
  }
  
  # Check if files exist
  if (!file.exists(fasta_file)) {
    stop(sprintf("FASTA file does not exist: %s", fasta_file))
  }
  if (!file.exists(gtf_file)) {
    stop(sprintf("GTF file does not exist: %s", gtf_file))
  }
  
  # Set output directory
  ref_output_dir <- file.path(output_dir, "cellranger_reference")
  
  # Build cellranger mkref command
  cmd <- sprintf("%s mkref --genome=%s --fasta=%s --genes=%s",
                 cellranger_path, genome_name, fasta_file, gtf_file)
  
  # If additional thread configuration provided
  nthreads <- get_nested_config(config, "cellranger.mkref.nthreads", default = NULL)
  if (!is.null(nthreads)) {
    cmd <- paste(cmd, sprintf("--nthreads=%d", nthreads))
  }
  
  log_info(sprintf("Running CellRanger mkref command: %s", cmd))
  log_info("This may take a long time, please be patient...")
  
  # Execute command
  result <- system(cmd, intern = FALSE)
  
  if (result != 0) {
    log_error("CellRanger mkref execution failed")
    stop("CellRanger mkref failed")
  }
  
  # Return reference genome path
  ref_genome_path <- file.path(getwd(), genome_name)
  log_info(sprintf("Reference genome created successfully: %s", ref_genome_path))
  
  return(ref_genome_path)
}

# Run cellranger count for quantification analysis
run_cellranger_count <- function(config, output_dir, transcriptome_path = NULL) {
  log_step("0.2", "CellRanger count - Gene expression quantification")
  
  # Get configuration parameters
  cellranger_path <- get_nested_config(config, "cellranger.cellranger_path", default = "cellranger")
  sample_id <- get_nested_config(config, "cellranger.count.sample_id", default = "sample")
  fastqs_dir <- get_nested_config(config, "cellranger.count.fastqs_dir", default = NULL)
  sample_name <- get_nested_config(config, "cellranger.count.sample_name", default = NULL)
  
  if (is.null(fastqs_dir)) {
    log_error("Missing required parameter: fastqs_dir")
    stop("CellRanger count requires fastqs_dir")
  }
  
  # Check if FASTQ directory exists
  if (!dir.exists(fastqs_dir)) {
    stop(sprintf("FASTQ directory does not exist: %s", fastqs_dir))
  }
  
  # Get or use provided transcriptome path
  if (is.null(transcriptome_path)) {
    transcriptome_path <- get_nested_config(config, "cellranger.count.transcriptome", default = NULL)
    if (is.null(transcriptome_path)) {
      log_error("Missing required parameter: transcriptome")
      stop("CellRanger count requires transcriptome path")
    }
  }
  
  # Check if transcriptome path exists
  if (!dir.exists(transcriptome_path)) {
    stop(sprintf("Transcriptome path does not exist: %s", transcriptome_path))
  }
  
  # Build cellranger count command
  cmd <- sprintf("%s count --id=%s --transcriptome=%s --fastqs=%s",
                 cellranger_path, sample_id, transcriptome_path, fastqs_dir)
  
  # Add sample name (if provided)
  if (!is.null(sample_name)) {
    cmd <- paste(cmd, sprintf("--sample=%s", sample_name))
  }
  
  # Add optional parameters
  force_cells <- get_nested_config(config, "cellranger.count.force_cells", default = NULL)
  if (!is.null(force_cells)) {
    cmd <- paste(cmd, sprintf("--force-cells=%d", force_cells))
  }
  
  expect_cells <- get_nested_config(config, "cellranger.count.expect_cells", default = NULL)
  if (!is.null(expect_cells)) {
    cmd <- paste(cmd, sprintf("--expect-cells=%d", expect_cells))
  }
  
  localcores <- get_nested_config(config, "cellranger.count.localcores", default = NULL)
  if (!is.null(localcores)) {
    cmd <- paste(cmd, sprintf("--localcores=%d", localcores))
  }
  
  localmem <- get_nested_config(config, "cellranger.count.localmem", default = NULL)
  if (!is.null(localmem)) {
    cmd <- paste(cmd, sprintf("--localmem=%d", localmem))
  }
  
  create_bam <- get_nested_config(config, "cellranger.count.create_bam", default = TRUE)
  if (!create_bam) {
    cmd <- paste(cmd, "--create-bam=false")
  }
  
  log_info(sprintf("Running CellRanger count command: %s", cmd))
  log_info("This may take a long time (several hours), please be patient...")
  
  # Execute command
  result <- system(cmd, intern = FALSE)
  
  if (result != 0) {
    log_error("CellRanger count execution failed")
    stop("CellRanger count failed")
  }
  
  # Return output matrix path
  matrix_path <- file.path(getwd(), sample_id, "outs", "filtered_feature_bc_matrix")
  
  if (!dir.exists(matrix_path)) {
    log_error(sprintf("CellRanger output directory does not exist: %s", matrix_path))
    stop("CellRanger count output directory not found")
  }
  
  log_info(sprintf("CellRanger count completed, output matrix path: %s", matrix_path))
  
  return(matrix_path)
}

# Complete CellRanger pipeline (mkref + count)
run_cellranger_pipeline <- function(config, output_dir) {
  log_info("Starting complete CellRanger pipeline")
  
  # Check if CellRanger is enabled
  cellranger_enabled <- get_nested_config(config, "cellranger.enabled", default = FALSE)
  if (!cellranger_enabled) {
    log_info("CellRanger not enabled, skipping")
    return(NULL)
  }
  
  # Check if need to run mkref
  run_mkref <- get_nested_config(config, "cellranger.mkref.enabled", default = FALSE)
  transcriptome_path <- NULL
  
  if (run_mkref) {
    log_info("Running CellRanger mkref...")
    transcriptome_path <- run_cellranger_mkref(config, output_dir)
  } else {
    log_info("Skipping CellRanger mkref, using existing reference genome")
    transcriptome_path <- get_nested_config(config, "cellranger.count.transcriptome", default = NULL)
  }
  
  # Check if need to run count
  run_count <- get_nested_config(config, "cellranger.count.enabled", default = FALSE)
  matrix_path <- NULL
  
  if (run_count) {
    log_info("Running CellRanger count...")
    matrix_path <- run_cellranger_count(config, output_dir, transcriptome_path)
    
    # Update data path in configuration for subsequent analysis
    config$data$data_dir <- matrix_path
    log_info(sprintf("Updated data path to CellRanger output: %s", matrix_path))
  } else {
    log_info("Skipping CellRanger count, using existing matrix files")
  }
  
  return(list(transcriptome = transcriptome_path, matrix = matrix_path, config = config))
}

# ============================================================================
# Data Loading Module
# ============================================================================

# Single sample data loading function
load_sc_data <- function(config, output_dir) {
  log_step("1", "Data loading started")
  
  data_dir <- get_nested_config(config, "data.data_dir", default = "./data")
  project_name <- get_nested_config(config, "project.project_name", default = "PanFamily_scRNA")
  min_cells <- get_nested_config(config, "quality_control.min_cells_per_feature", default = 3)
  min_features <- get_nested_config(config, "quality_control.min_features_per_cell", default = 200)
  
  # Read 10X data
  log_info(sprintf("Loading data from directory: %s", data_dir))
  raw_data <- Read10X(data.dir = data_dir)
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = raw_data,
    project = project_name,
    min.cells = min_cells,
    min.features = min_features
  )
  
  log_info(sprintf("Seurat object created successfully: %d genes, %d cells",
                   nrow(seurat_obj), ncol(seurat_obj)))
  
  # No need to save original object, reloading is fast
  
  return(seurat_obj)
}

# ============================================================================
# Quality Control Module
# ============================================================================

run_quality_control <- function(seurat_obj, config, output_dir) {
  log_step("2", "Quality control started")
  
  # Calculate QC metrics
  mt_pattern <- get_nested_config(config, "quality_control.plant_specific.mitochondrial_pattern", default = "^ATMG")
  cg_pattern <- get_nested_config(config, "quality_control.plant_specific.chloroplast_pattern", default = "^ATCG")
  
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mt_pattern)
  seurat_obj[["percent.cg"]] <- PercentageFeatureSet(seurat_obj, pattern = cg_pattern)
  
  # Visualization before QC
  plot_qc_before <- VlnPlot(seurat_obj, 
                           features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.cg"), 
                           ncol = 4)
  qc_before_file <- file.path(output_dir, "QC_before_filtering.png")
  save_plot(plot_qc_before, qc_before_file, width = 20, height = 10)
  
  # Generate data statistics summary
  stats_file <- generate_qc_summary(seurat_obj, output_dir)
  
  # Check if interaction is needed
  execution_mode <- get_nested_config(config, "execution.mode", default = "interactive")
  quality_control_interactive <- get_nested_config(config, "execution.interactive_settings.steps.quality_control", default = FALSE)
  
  log_info(sprintf("Execution mode: %s, QC interactive: %s", execution_mode, quality_control_interactive))
  
  if (execution_mode == "interactive" && quality_control_interactive) {
    log_info("QC parameters need user interactive adjustment")
    
    # Display statistics content
    cat("\n=== Interactive Prompt: QC Parameter Adjustment ===\n")
    cat("Pre-QC plot: ", qc_before_file, "\n", sep = "")
    cat("\n--- Data Statistics Summary ---\n")
    if (!is.null(stats_file) && file.exists(stats_file)) {
      stats_content <- readLines(stats_file)
      cat(paste(stats_content, collapse = "\n"), "\n")
    } else {
      cat("Statistics file generation failed\n")
    }
    
    cat("\n--- Current Filtering Parameters ---\n")
    cat("min_features_per_cell: ", get_nested_config(config, "quality_control.min_features_per_cell", default = 200), "\n", sep = "")
    cat("max_features_per_cell: ", get_nested_config(config, "quality_control.max_features_per_cell", default = 2500), "\n", sep = "")
    cat("max_mt_percent: ", get_nested_config(config, "quality_control.max_mt_percent", default = 5.0), "\n", sep = "")
    cat("max_count_per_cell: ", get_nested_config(config, "quality_control.max_count_per_cell", default = 20000), "\n", sep = "")
    
    # User choice with timeout (120 seconds)
    choice <- get_user_input_with_timeout("Do you need to modify QC parameters? (Y/N, will use current parameters if no response in 120 seconds): ", 120)
    
    if (toupper(choice) == "Y") {
      # Interactive parameter modification
      cat("\n=== Interactive Parameter Modification ===\n")
      cat("Enter new parameter values (press Enter directly to keep current value):\n\n")
      
      # Get new parameters (with timeout)
      input_min <- get_user_input_with_timeout(sprintf("min_features_per_cell [current:%d, will keep current value if no response in 120 seconds]: ", config$quality_control$min_features_per_cell), 120)
      new_min_features <- if (input_min == "TIMEOUT" || input_min == "") config$quality_control$min_features_per_cell else as.integer(input_min)
      
      input_max <- get_user_input_with_timeout(sprintf("max_features_per_cell [current:%d, will keep current value if no response in 120 seconds]: ", config$quality_control$max_features_per_cell), 120)
      new_max_features <- if (input_max == "TIMEOUT" || input_max == "") config$quality_control$max_features_per_cell else as.integer(input_max)
      
      input_mt <- get_user_input_with_timeout(sprintf("max_mt_percent [current:%.1f, will keep current value if no response in 120 seconds]: ", config$quality_control$max_mt_percent), 120)
      new_max_mt <- if (input_mt == "TIMEOUT" || input_mt == "") config$quality_control$max_mt_percent else as.numeric(input_mt)
      
      input_count <- get_user_input_with_timeout(sprintf("max_count_per_cell [current:%d, will keep current value if no response in 120 seconds]: ", config$quality_control$max_count_per_cell), 120)
      new_max_count <- if (input_count == "TIMEOUT" || input_count == "") config$quality_control$max_count_per_cell else as.integer(input_count)
      
      # Update configuration and save to yaml file
      config$quality_control$min_features_per_cell <- new_min_features
      config$quality_control$max_features_per_cell <- new_max_features  
      config$quality_control$max_mt_percent <- new_max_mt
      config$quality_control$max_count_per_cell <- new_max_count
      
      # Save updated configuration to yaml file
      yaml::write_yaml(config, opt$config)
      
      cat("\n=== Updated Parameters ===\n")
      cat("min_features_per_cell:", new_min_features, "\n")
      cat("max_features_per_cell:", new_max_features, "\n") 
      cat("max_mt_percent:", new_max_mt, "\n")
      cat("max_count_per_cell:", new_max_count, "\n")
      log_info("Parameters updated and saved to configuration file")
      
    } else {
      log_info("Continuing with current parameters")
    }
  }
  
  # Apply filtering criteria
  log_info("Applying QC filtering criteria")
  cells_before <- ncol(seurat_obj)
  
  log_info(sprintf("Using filtering parameters: min_features=%d, max_features=%d, max_mt=%.1f, max_count=%d",
                   config$quality_control$min_features_per_cell,
                   config$quality_control$max_features_per_cell,
                   config$quality_control$max_mt_percent,
                   config$quality_control$max_count_per_cell))
  
  seurat_obj <- subset(seurat_obj,
    subset = nFeature_RNA > config$quality_control$min_features_per_cell &
             nFeature_RNA < config$quality_control$max_features_per_cell &
             percent.mt < config$quality_control$max_mt_percent &
             nCount_RNA < config$quality_control$max_count_per_cell
  )
  
  cells_after <- ncol(seurat_obj)
  log_info(sprintf("Filtering completed: %d -> %d cells (%.1f%% retained)",
                   cells_before, cells_after, 100 * cells_after / cells_before))
  
  # Visualization after QC
  plot_qc_after <- VlnPlot(seurat_obj, 
                          features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.cg"), 
                          ncol = 4)
  save_plot(plot_qc_after, file.path(output_dir, "QC_after_filtering.png"), 
            width = 20, height = 10)
  
  # Filtering step is fast, no need to save intermediate object
  
  return(seurat_obj)
}

# ============================================================================
# Normalization and Dimension Reduction Module
# ============================================================================

run_normalization <- function(seurat_obj, config, output_dir) {
  log_step("3", "Data normalization started")
  
  # Seurat 5.0 compatibility handling
  log_info("Starting normalization workflow")
  
  # Normalize data
  seurat_obj <- NormalizeData(seurat_obj, 
                             normalization.method = "LogNormalize",
                             scale.factor = config$seurat$normalization$scale_factor)
  
  # Find highly variable genes (adjust based on gene count)
  total_genes <- nrow(seurat_obj)
  n_features_to_find <- min(config$seurat$variable_features$n_features, total_genes)
  
  if (total_genes < 50) {
    # For small gene counts (like gene family analysis), use all genes
    log_info(sprintf("Low gene count (%d genes), using all genes as variable features", total_genes))
    VariableFeatures(seurat_obj) <- rownames(seurat_obj)
  } else {
    # Normal case: find highly variable genes
  seurat_obj <- tryCatch({
    FindVariableFeatures(seurat_obj, 
                        selection.method = config$seurat$variable_features$method,
                          nfeatures = n_features_to_find,
                        layer = "data")
  }, error = function(e) {
    log_info("Using backup method to find highly variable genes")
    FindVariableFeatures(seurat_obj, 
                        selection.method = config$seurat$variable_features$method,
                          nfeatures = n_features_to_find)
  })
  }
  
  # Visualization of highly variable genes
  top10_hvg <- head(VariableFeatures(seurat_obj), 10)
  plot_hvg <- LabelPoints(plot = VariableFeaturePlot(seurat_obj), 
                         points = top10_hvg, repel = TRUE)
  save_plot(plot_hvg, file.path(output_dir, "variable_features.png"))
  
  # Scale data
  all.genes <- rownames(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj, features = all.genes)
  
  # Normalization step is fast, no need to save intermediate object
  
  return(seurat_obj)
}

run_dimension_reduction <- function(seurat_obj, config, output_dir) {
  log_step("4", "Dimension reduction analysis started")
  
  # Adjust PCA parameters based on gene count
  total_genes <- length(VariableFeatures(seurat_obj))
  n_pcs_compute <- min(config$seurat$pca$n_pcs_compute, total_genes - 1, ncol(seurat_obj) - 1)
  
  if (total_genes < 10) {
    log_info(sprintf("Too few genes (%d genes), skipping PCA dimension reduction, using raw expression matrix directly", total_genes))
    # For gene family analysis, create simplified dimension reduction result
    embed_matrix <- t(as.matrix(GetAssayData(seurat_obj, layer = "data")))
    seurat_obj[["pca"]] <- CreateDimReducObject(embeddings = embed_matrix[,1:min(2, ncol(embed_matrix))], 
                                               key = "PC_", assay = DefaultAssay(seurat_obj))
  } else {
    # PCA analysis
    log_info(sprintf("Computing %d principal components", n_pcs_compute))
    seurat_obj <- RunPCA(seurat_obj, 
                        features = VariableFeatures(seurat_obj),
                        npcs = n_pcs_compute)
  }
  
  # ElbowPlot (adjust ndims)
  max_dims <- min(50, n_pcs_compute)
  if (max_dims >= 3) {
    plot_elbow <- ElbowPlot(seurat_obj, ndims = max_dims)
  elbow_file <- file.path(output_dir, "PCA_ElbowPlot.png")
  save_plot(plot_elbow, elbow_file)
  } else {
    log_info("Too few principal components, skipping ElbowPlot")
    elbow_file <- "skipped_due_to_low_gene_count"
  }
  
  # PCA loadings
  plot_loadings <- VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca")
  save_plot(plot_loadings, file.path(output_dir, "PCA_loadings.png"))
  
  # Check if PCA interaction is needed
  execution_mode <- get_nested_config(config, "execution.mode", default = "interactive")
  pca_interactive <- get_nested_config(config, "execution.interactive_settings.steps.pca_selection", default = FALSE)
  
  if (execution_mode == "interactive" && pca_interactive) {
    log_info("PCA principal component count needs user interactive selection")
    
    # Display current configuration
    cat("\n=== Interactive Prompt: PCA Principal Component Count Selection ===\n")
    cat("ElbowPlot image: ", elbow_file, "\n", sep = "")
    cat("\n--- Current Configuration ---\n")
    cat("n_pcs_use: ", config$seurat$pca$n_pcs_use, "\n", sep = "")
    
    # User choice with timeout (120 seconds)
    choice <- get_user_input_with_timeout("Do you need to modify PCA principal component count? (Y/N, will use current parameters if no response in 120 seconds): ", 120)
    
    if (toupper(choice) == "Y") {
      cat("\n=== Principal Component Count Modification ===\n")
      input_pcs <- get_user_input_with_timeout(sprintf("Enter new principal component count [current:%d, will keep current value if no response in 120 seconds]: ", config$seurat$pca$n_pcs_use), 120)
      
      if (input_pcs != "TIMEOUT" && input_pcs != "") {
        new_n_pcs <- as.integer(input_pcs)
        config$seurat$pca$n_pcs_use <- new_n_pcs
        yaml::write_yaml(config, opt$config)
        cat("Principal component count updated to:", new_n_pcs, "\n")
        log_info("PCA parameters updated and saved to configuration file")
      } else {
        log_info("Keeping current PCA parameters unchanged")
      }
    } else {
      log_info("Continuing with current PCA parameters")
    }
  }
  
  # Clustering (adjust principal component count)
  n_pcs_use <- min(config$seurat$pca$n_pcs_use, n_pcs_compute)
  
  if (total_genes < 10) {
    # Too few genes, use all genes for clustering
    log_info("Too few genes, using raw expression data for clustering")
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:min(2, n_pcs_compute))
  } else {
    log_info("Using PCA space for clustering")
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:n_pcs_use)
  }
  
  seurat_obj <- FindClusters(seurat_obj, resolution = config$seurat$clustering$resolution)
  
  # UMAP and t-SNE (adjust dimensions)
  dims_for_umap <- 1:min(n_pcs_use, n_pcs_compute)
  
  if (length(dims_for_umap) >= 2) {
    log_info("Using PCA space for UMAP and t-SNE dimension reduction")
    seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = dims_for_umap)
    seurat_obj <- RunTSNE(seurat_obj, reduction = "pca", dims = dims_for_umap)
  } else {
    log_info("Insufficient principal components, skipping UMAP and t-SNE analysis")
    return(seurat_obj)
  }
  
  # Visualization - use beautified DimPlot (reference ggsci color scheme)
  # Create rich color scheme
  mycol <- c(pal_d3()(7), pal_aaas()(7), pal_uchicago()(7), pal_jama()(7))
  
  # UMAP plot - beautified colors, no grid background, square aspect ratio, corrected axes, interactive sizing
  plot_umap <- interactive_plot_sizing(
    plot_func = function() {
      DimPlot(seurat_obj, 
              reduction = "umap",  
              group.by = "seurat_clusters",
              label = TRUE,
              label.size = 4,
              pt.size = 0.5,
              cols = mycol) +
        theme_classic(base_size = 14) +
        theme(
          legend.position = "right",
          panel.border = element_rect(fill = NA, colour = "black", size = 0.8),
          aspect.ratio = 1
        ) +
        coord_fixed(ratio = 1) +
        labs(x = "", y = "") +
        ggtitle("UMAP (clusters)", subtitle = paste0("nCells = ", ncol(seurat_obj))) +
        theme(plot.margin = margin(0, 0, 0, 0)) +
        coord_flip() +
        scale_y_reverse()
    },
    plot_name = "UMAP_clusters",
    default_width = 8,
    default_height = 8,
    output_dir = output_dir,
    config = config
  )
  
  # t-SNE plot - beautified colors, no grid background, square aspect ratio, corrected axes, interactive sizing
  plot_tsne <- interactive_plot_sizing(
    plot_func = function() {
      DimPlot(seurat_obj,
              reduction = "tsne",  
              group.by = "seurat_clusters",
              label = TRUE,
              label.size = 4,
              pt.size = 0.5,
              cols = mycol) +
        theme_classic(base_size = 14) +
        theme(
          legend.position = "right",
          panel.border = element_rect(fill = NA, colour = "black", size = 0.8),
          aspect.ratio = 1
        ) +
        coord_fixed(ratio = 1) +
        labs(x = "", y = "") +
        ggtitle("t-SNE (clusters)", subtitle = paste0("nCells = ", ncol(seurat_obj))) +
        theme(plot.margin = margin(0, 0, 0, 0)) +
        coord_flip() +
        scale_y_reverse()
    },
    plot_name = "TSNE_clusters",
    default_width = 8,
    default_height = 8,
    output_dir = output_dir,
    config = config
  )
  
  # Save clustered result object - key checkpoint, contains complete analysis results
  saveRDS(seurat_obj, file.path(output_dir, "clustered_seurat_object.rds"))
  log_info("Clustered object saved, can continue analysis from here using continue_analysis.R")
  
  return(seurat_obj)
}

# ============================================================================
# Marker Gene Identification
# ============================================================================

find_marker_genes <- function(seurat_obj, config, output_dir) {
  log_step("5", "Marker gene identification started")
  
  # Find all marker genes
  markers_all <- FindAllMarkers(seurat_obj, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(markers_all, file.path(output_dir, "all_markers.csv"), row.names = FALSE)
  
  # Find upregulated marker genes
  markers_pos <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(markers_pos, file.path(output_dir, "positive_markers.csv"), row.names = FALSE)
  
  # Top markers heatmap
  if (nrow(markers_pos) > 0) {
    top_markers <- markers_pos %>% 
      group_by(cluster) %>% 
      top_n(n = 5, wt = avg_log2FC)
    
    if (nrow(top_markers) > 0) {
      plot_heatmap <- DoHeatmap(seurat_obj, features = top_markers$gene) + NoLegend()
      save_plot(plot_heatmap, file.path(output_dir, "top_markers_heatmap.png"), 
                width = 12, height = 10)
    }
  }
  
  return(list(all = markers_all, positive = markers_pos))
}

# ============================================================================
# Cell Annotation Module
# ============================================================================

run_cell_annotation <- function(seurat_obj, markers, config, output_dir) {
  log_step("6", "Cell type annotation started")
  
  annotation_method <- config$cell_annotation$method
  
  if (annotation_method == "scMayoMap") {
    log_info("Using scMayoMap for cell annotation")
    
    # Load annotation database
    db_file <- config$cell_annotation$scmayomap$database_file
    if (!file.exists(db_file)) {
      log_error(sprintf("Annotation database file does not exist: %s", db_file))
      return(seurat_obj)
    }
    
    # Prepare reference data
    demodata <- read.xlsx(db_file)
    unique_demodata <- demodata %>% distinct(gene, celltype, value, .keep_all = TRUE)
    db <- tidyr::spread(unique_demodata, key = c('celltype'), value = 'value')
    db[is.na(db)] <- 0
    
    # Run scMayoMap
    scMayoMap.obj <- scMayoMap(data = markers$positive, database = db)
    res <- scMayoMap.obj$res
    
    # Visualize scoring results
    plot_score <- scMayoMap.plot(scMayoMap.object = scMayoMap.obj)
    save_plot(plot_score, file.path(output_dir, "scmayomap_scores.png"))
    
    # Get cell type prediction results
    celltype <- scMayoMap.obj$markers %>% 
      group_by(cluster) %>% 
      slice_max(score, n = 1)
    
    result <- celltype %>%
      group_by(cluster) %>%
      filter(score == max(score)) %>%
      slice_head(n = 1) %>%
      ungroup()
    
    # Add celltype to meta.data
    for(i in 1:nrow(result)){
      seurat_obj@meta.data[which(seurat_obj$seurat_clusters == result$cluster[i]), "celltype"] <- result$celltype[i]
    }
    
    # Visualize annotation results - use beautified DimPlot (reference ggsci color scheme)
    # Create rich color scheme
    mycol_celltype <- c(pal_d3()(7), pal_aaas()(7), pal_uchicago()(7), pal_jama()(7))
    
    # Cell type UMAP plot - beautified colors, no grid background, square aspect ratio, corrected axes
    plot_celltype <- DimPlot(seurat_obj, 
                             reduction = "umap",  
                             group.by = "celltype",
                             label = TRUE,
                             label.size = 4,
                             pt.size = 0.5,
                             cols = mycol_celltype) +
                     theme_classic(base_size = 14) +
                     theme(
                       legend.position = "right",
                       panel.border = element_rect(fill = NA, colour = "black", size = 0.8),
                       aspect.ratio = 1
                     ) +
                     coord_fixed(ratio = 1) +
                     scale_x_reverse() +  # Reverse x-axis
                     scale_y_reverse() +  # Reverse y-axis
                     labs(x = "UMAP 1", y = "UMAP 2") +
                     ggtitle("UMAP (cell types)", subtitle = paste0("nCells = ", ncol(seurat_obj)))
    save_plot(plot_celltype, file.path(output_dir, "UMAP_celltypes.png"), 
              width = 8, height = 8)

    
    # Cell type proportion statistics
    celltype_counts <- table(seurat_obj$celltype)
    df <- as.data.frame(celltype_counts)
    names(df) <- c("celltype", "Freq")
    
    # Create labels (including percentages)
    labs <- paste0(df$celltype, " (", round(df$Freq / sum(df$Freq) * 100, 2), "%)")
    
    doughnut <- function(x, labels = names(x), outer.radius = 0.8, inner.radius = 0.4, col = NULL, border = "white") {
      pie(x, labels = labels, radius = outer.radius, col = col, border = border)
      symbols(0, 0, circles = inner.radius, inches = FALSE, add = TRUE, bg = "white")
    }
    
    tryCatch({
      library(RColorBrewer)
      colors <- brewer.pal(min(length(df$Freq), 11), "Set3")  # Set3 supports up to 12 colors
      if(length(df$Freq) > 11) {
        colors <- c(colors, rainbow(length(df$Freq) - 11))
      }
    }, error = function(e) {
      colors <- rainbow(length(df$Freq))
    })
    png(file.path(output_dir, "celltype_ratio.png"), width = 15, height = 10, 
        units = "in", res = 300)
    doughnut(df$Freq, labels = labs, col = colors, inner.radius = 0.4)
    dev.off()
  }
  
  # Save cell type annotation object - starting point for pseudotime analysis
  saveRDS(seurat_obj, file.path(output_dir, "annotated_seurat_object.rds"))
  log_info("Annotation object saved, contains cell type information, for pseudotime analysis")
  
  return(seurat_obj)
}


# Stacked violin plot function removed, expression levels too low for good visualization

# ============================================================================
# Gene Family Specificity Analysis
# ============================================================================

analyze_target_genes <- function(seurat_obj, config, output_dir) {
  log_step("7", "Gene family expression analysis started")
  
  # Read target gene list
  gene_list_file <- config$data$target_gene_list
  if (is.null(gene_list_file) || gene_list_file == "null" || !file.exists(gene_list_file)) {
    log_info("No gene family list provided, skipping specificity analysis")
    return(seurat_obj)
  }
  
  target_genes <- readLines(gene_list_file)
  target_genes <- target_genes[target_genes != ""]  # Remove empty lines
  
  # Check which genes exist in data
  available_genes <- target_genes[target_genes %in% rownames(seurat_obj)]
  missing_genes <- target_genes[!target_genes %in% rownames(seurat_obj)]
  
  log_info(sprintf("Total target genes: %d, available in data: %d, missing: %d",
                   length(target_genes), length(available_genes), length(missing_genes)))
  
  if (length(missing_genes) > 0) {
    writeLines(missing_genes, file.path(output_dir, "missing_genes.txt"))
  }
  
  if (length(available_genes) == 0) {
    log_error("No target genes found in data")
    return(seurat_obj)
  }
  
  # Feature plot - use configuration parameters
  ncol <- config$target_gene_analysis$plot_params$ncol
  feature_width <- config$target_gene_analysis$plot_params$feature_plot$width
  feature_height <- config$target_gene_analysis$plot_params$feature_plot$height
  
  log_info(sprintf("FeaturePlot: %d genes, %d per row, plot dimensions: %d x %d",
                   length(available_genes), ncol, feature_width, feature_height))
  
  # Interactive plot sizing
  plot_feature <- interactive_plot_sizing(
    plot_func = function() {
      FeaturePlot(seurat_obj, 
                  features = available_genes,
                  reduction = "umap",
                  pt.size = 0.8,
                  cols = c("lightgrey", "blue"),
                  label = FALSE,
                  order = TRUE,
                  ncol = ncol) &
        labs(x = "", y = "") &
        theme(plot.margin = margin(0, 0, 0, 0)) &
        coord_flip() &
        scale_y_reverse()
    },
    plot_name = "target_genes_feature",
    default_width = feature_width,
    default_height = feature_height,
    output_dir = output_dir,
    config = config
  )
  
  # Violin plot removed, expression levels too low for good visualization
  
  return(seurat_obj)
}

# ============================================================================
# Pseudotime Analysis Module
# ============================================================================

run_pseudotime_analysis <- function(seurat_obj, config, output_dir) {
  if (!config$pseudotime$enabled) {
    log_info("Pseudotime analysis disabled")
    return(NULL)
  }
  
  log_step("8", "Pseudotime analysis started")
  
  # Create Monocle2 object from Seurat object (reference family_sc.R)
  log_info("Creating Monocle2 object...")
  
  # Prepare expression matrix
  data <- as(as.matrix(seurat_obj[["RNA"]]$counts), 'sparseMatrix')
  
  # Prepare cell metadata
  pd <- new('AnnotatedDataFrame', data = seurat_obj@meta.data)
  
  # Prepare gene metadata
  fData <- data.frame(gene_short_name = row.names(data), 
                      row.names = row.names(data))
  fd <- new('AnnotatedDataFrame', data = fData)
  
  # Create CellDataSet object
  mycds <- newCellDataSet(data,
                          phenoData = pd,
                          featureData = fd,
                          expressionFamily = negbinomial.size())
  
  # Estimate dispersion of expression matrix
  log_info("Estimating size factors and dispersions...")
  mycds <- estimateSizeFactors(mycds)
  mycds <- estimateDispersions(mycds)
  
  # QC: filter gene expression levels
  log_info("Detecting gene expression...")
  mycds <- detectGenes(mycds, min_expr = config$pseudotime$monocle2_params$min_expr)
  
  # Selected genes expressed in at least 10 cells
  mycds_expressed_genes <- row.names(subset(fData(mycds),
                                          num_cells_expressed >= 10))
  
  # Differential gene analysis (reference family_sc.R)
  log_info("Performing differential gene analysis...")
  diff_test_res <- differentialGeneTest(mycds[mycds_expressed_genes,],
                                        fullModelFormulaStr = "~celltype",
                                        cores = 18)
  
  # Select significantly different genes as ordering genes
  ordering_genes <- row.names(subset(diff_test_res, qval < 0.05))
  mycds <- setOrderingFilter(mycds, ordering_genes)
  
  # Visualize ordering genes
  diff_ordering_genes <- plot_ordering_genes(mycds)
  save_plot(diff_ordering_genes, file.path(output_dir, "diff_ordering_genes.png"))
  
  # Dimension reduction analysis
  log_info("Dimension reduction analysis...")
  mycds <- reduceDimension(mycds,
                          max_components = config$pseudotime$monocle2_params$max_components,
                          method = config$pseudotime$monocle2_params$reduction_method)
  
  # Calculate cell pseudotime
  log_info("Calculating pseudotime...")
  mycds <- orderCells(mycds)
  saveRDS(mycds, file.path(output_dir, "monocle_cds.rds"))
  
  # Also save seurat_obj for branch analysis
  saveRDS(seurat_obj, file.path(output_dir, "seurat_for_branch.rds"))
  # Visualize results (reference multiple visualizations in family_sc.R)
  log_info("Generating visualization plots...")
  
  # 1. Trajectory plot colored by cell type, corrected axes, interactive sizing
  pse_celltype <- interactive_plot_sizing(
    plot_func = function() {
      plot_cell_trajectory(mycds, color_by = "celltype") +
        theme_classic(base_size = 14) +
        theme(
          legend.position = "right",
          panel.border = element_rect(fill = NA, colour = "black", size = 0.8),
          aspect.ratio = 1
        ) +
        coord_fixed(ratio = 1) +
        labs(x = "", y = "") +
        ggtitle("Pseudotime (cell types)", subtitle = paste0("nCells = ", ncol(mycds))) +
        theme(plot.margin = margin(0, 0, 0, 0)) +
        coord_flip() +
        scale_y_reverse()
    },
    plot_name = "pseudotime_celltype",
    default_width = 8,
    default_height = 6,
    output_dir = output_dir,
    config = config
  )
  
  # 2. Trajectory plot colored by state, corrected axes, interactive sizing
  pse_state <- interactive_plot_sizing(
    plot_func = function() {
      plot_cell_trajectory(mycds, color_by = "State") +
        theme_classic(base_size = 14) +
        theme(
          legend.position = "right",
          panel.border = element_rect(fill = NA, colour = "black", size = 0.8),
          aspect.ratio = 1
        ) +
        coord_fixed(ratio = 1) +
        labs(x = "", y = "") +
        ggtitle("Pseudotime (states)", subtitle = paste0("nCells = ", ncol(mycds))) +
        theme(plot.margin = margin(0, 0, 0, 0)) +
        coord_flip() +
        scale_y_reverse()
    },
    plot_name = "pseudotime_state",
    default_width = 8,
    default_height = 6,
    output_dir = output_dir,
    config = config
  )
  
  # 3. Trajectory plot colored by pseudotime, corrected axes, interactive sizing
  pse_pseudotime <- interactive_plot_sizing(
    plot_func = function() {
      plot_cell_trajectory(mycds, color_by = "Pseudotime") +
        theme_classic(base_size = 14) +
        theme(
          legend.position = "right",
          panel.border = element_rect(fill = NA, colour = "black", size = 0.8),
          aspect.ratio = 1
        ) +
        coord_fixed(ratio = 1) +
        labs(x = "", y = "") +
        ggtitle("Pseudotime (trajectory)", subtitle = paste0("nCells = ", ncol(mycds))) +
        theme(plot.margin = margin(0, 0, 0, 0)) +
        coord_flip() +
        scale_y_reverse()
    },
    plot_name = "pseudotime_trajectory",
    default_width = 8,
    default_height = 6,
    output_dir = output_dir,
    config = config
  )
  
  # 4. Combined plot, interactive sizing
  pse_combined <- interactive_plot_sizing(
    plot_func = function() {
      pse_celltype | pse_state | pse_pseudotime
    },
    plot_name = "pseudotime_combined",
    default_width = 15,
    default_height = 5,
    output_dir = output_dir,
    config = config
  )
  
  # 5. Faceted display by state, interactive sizing
  pse_state_facet <- interactive_plot_sizing(
    plot_func = function() {
      plot_cell_trajectory(mycds, color_by = "State") +
        facet_wrap(~State, nrow = 4) +
        labs(x = "", y = "") +
        ggtitle("Pseudotime by State", subtitle = paste0("nCells = ", ncol(mycds))) +
        theme(plot.margin = margin(0, 0, 0, 0)) +
        coord_flip() +
        scale_y_reverse()
    },
    plot_name = "pseudotime_state_facet",
    default_width = 12,
    default_height = 10,
    output_dir = output_dir,
    config = config
  )
  
  # 6. Faceted display by cell type, interactive sizing
  pse_celltype_facet <- interactive_plot_sizing(
    plot_func = function() {
      plot_cell_trajectory(mycds, color_by = "State") +
        facet_wrap(~celltype, nrow = 4) +
        labs(x = "", y = "") +
        ggtitle("Pseudotime by Cell Type", subtitle = paste0("nCells = ", ncol(mycds))) +
        theme(plot.margin = margin(0, 0, 0, 0)) +
        coord_flip() +
        scale_y_reverse()
    },
    plot_name = "pseudotime_celltype_facet",
    default_width = 12,
    default_height = 10,
    output_dir = output_dir,
    config = config
  )
  
  # Target gene expression analysis in pseudotime (before branch analysis)
  if (config$target_gene_analysis$enabled && 
      !is.null(config$target_gene_analysis$gene_list_file)) {
    
    log_info("Analyzing target gene expression in pseudotime (before branch analysis)...")
    
    # Read target genes
    target_genes <- readLines(config$target_gene_analysis$gene_list_file)
    target_genes <- trimws(target_genes)
    target_genes <- target_genes[target_genes != ""]
    
    # Find target genes that exist in data
    all_genes <- diff_test_res$gene_short_name
    existing_genes <- target_genes[target_genes %in% all_genes]
    
    if (length(existing_genes) > 0) {
      log_info(paste("Found", length(existing_genes), "target genes in pseudotime data"))
      
      # Jitter plot
      p1 <- plot_genes_jitter(mycds[existing_genes,],
                             grouping = "State",
                             min_expr = 0.1,
                             color_by = "celltype",
                             ncol = 3)
      save_plot(p1, file.path(output_dir, "target_genes_jitter_before_branch.png"))
      
      # Pseudotime plot
      p2 <- plot_genes_in_pseudotime(mycds[existing_genes,], 
                                    color_by = "celltype", 
                                    ncol = 3)
      save_plot(p2, file.path(output_dir, "target_genes_pseudotime_before_branch.png"))
      
      # Violin plot removed, expression levels too low for good visualization
      
      # Combined plot
      plot_combined <- p1 | p2
      save_plot(plot_combined, file.path(output_dir, "target_genes_combined_before_branch.png"), 
                width = 20, height = 6)
      
      # Heatmap display
      png(file.path(output_dir, "target_genes_heatmap_before_branch.png"), 
          width = 6, height = 5, units = "in", res = 300)
      plot_pseudotime_heatmap(mycds[existing_genes,],
                             num_clusters = 6,
                             cores = 1,
                             show_rownames = T)
      dev.off()
      
    } else {
      log_info("No target genes found in pseudotime data")
    }
  }
  
  # 4. Branch analysis
  if (config$pseudotime$branch_analysis$enabled) {
    log_info("Performing branch analysis...")
    
    # Get current configuration parameters for branch analysis
    current_branch_point <- get_nested_config(config, "pseudotime.monocle2_params.branch_point", default = 1)
    current_cores <- get_nested_config(config, "resources.cores", default = 4)
    
    # Check if interactive adjustment of branch analysis parameters is needed
    execution_mode <- get_nested_config(config, "execution.mode", default = "interactive")
    branch_interactive <- get_nested_config(config, "execution.interactive_settings.steps.branch_analysis", default = FALSE)
    
    if (execution_mode == "interactive" && branch_interactive) {
      log_info("Branch analysis parameters need user interactive adjustment")
      
      cat("\n=== Interactive Prompt: Branch Analysis Parameter Settings ===\n")
      cat("Current configuration parameters:\n")
      cat("branch_point: ", current_branch_point, "\n", sep = "")
      cat("cores: ", current_cores, "\n", sep = "")
      cat("progenitor_method: duplicate (fixed)\n")
      
      # User choice with timeout (120 seconds)
      choice <- get_user_input_with_timeout("Do you need to modify branch analysis parameters? (Y/N, will use current parameters if no response in 120 seconds): ", 120)
      
      if (toupper(choice) == "Y") {
        cat("\n=== Branch Analysis Parameter Modification ===\n")
        
        # Get new branch point parameter (with timeout)
        input_branch <- get_user_input_with_timeout(sprintf("Enter branch point [current:%d, will keep current value if no response in 120 seconds]: ", current_branch_point), 120)
        if (input_branch != "TIMEOUT" && input_branch != "") {
          new_branch_point <- as.integer(input_branch)
          if (!is.na(new_branch_point) && new_branch_point > 0) {
            current_branch_point <- new_branch_point
            config$pseudotime$monocle2_params$branch_point <- new_branch_point
          }
        }
        
        # Get new core count parameter (with timeout)
        input_cores <- get_user_input_with_timeout(sprintf("Enter number of cores to use [current:%d, will keep current value if no response in 120 seconds]: ", current_cores), 120)
        if (input_cores != "TIMEOUT" && input_cores != "") {
          new_cores <- as.integer(input_cores)
          if (!is.na(new_cores) && new_cores > 0 && new_cores <= 32) {  # Limit maximum core count
            current_cores <- new_cores
            config$resources$cores <- new_cores
          }
        }
        
        cat("\n=== Updated Parameters ===\n")
        cat("branch_point:", current_branch_point, "\n")
        cat("cores:", current_cores, "\n")
        log_info("Branch analysis parameters updated")
        
      } else {
        log_info("Continuing with current branch analysis parameters")
      }
    }
    
    # Run BEAM analysis
    log_info(sprintf("Running BEAM analysis with parameters: branch_point=%d, cores=%d", current_branch_point, current_cores))
    BEAM_res <- BEAM(mycds, 
                     branch_point = current_branch_point, 
                     progenitor_method = "duplicate",
                     cores = current_cores)
    BEAM_res <- BEAM_res[order(BEAM_res$qval),]
    BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
    
    # Save branch genes
    write.csv(BEAM_res, file.path(output_dir, "branch_genes.csv"), row.names = FALSE)
    
    # Save branch analysis results and working environment
    log_info("Saving branch analysis results and working environment...")
    saveRDS(BEAM_res, file.path(output_dir, "BEAM_results.rds"))
    saveRDS(mycds, file.path(output_dir, "monocle_cds_after_branch.rds"))
    if (exists("seurat_obj") && !is.null(seurat_obj)) {
      saveRDS(seurat_obj, file.path(output_dir, "seurat_after_branch.rds"))
    }
    
    # Visualize branch genes
    plot_genes_branched_heatmap(mycds[row.names(subset(BEAM_res, qval < 1e-4)),],
                                branch_point = current_branch_point,
                                num_clusters = 4,
                                cores = 1,
                                use_gene_short_name = T,
                                show_rownames = T)
    ggsave(file.path(output_dir, "branch_heatmap.png"), width = 10, height = 8)
    
    # Specific gene branch analysis
    if (config$target_gene_analysis$enabled && 
        !is.null(config$target_gene_analysis$gene_list_file)) {
      
      log_info("Starting specific gene branch analysis...")
      
      # Read target genes
      target_genes <- readLines(config$target_gene_analysis$gene_list_file)
      target_genes <- target_genes[target_genes != ""]  # Remove empty lines
      log_info(sprintf("Read %d target genes", length(target_genes)))
      
      # Debug information: show basic info about BEAM results and target genes
      log_info(sprintf("BEAM results contain %d genes", nrow(BEAM_res)))
      log_info(sprintf("BEAM result gene name examples: %s", paste(head(BEAM_res$gene_short_name, 3), collapse=", ")))
      log_info(sprintf("Target gene examples: %s", paste(head(target_genes, 3), collapse=", ")))
      
      # Find target genes in BEAM results
      target_beam <- BEAM_res[BEAM_res$gene_short_name %in% target_genes, ]
      
      if (nrow(target_beam) > 0) {
        log_info(sprintf("Found %d target genes in BEAM results", nrow(target_beam)))
        
        # Save target gene branch analysis results
        write.csv(target_beam, file.path(output_dir, "target_genes_branch_analysis.csv"), row.names = FALSE)
        
        # Select only significant target genes (qval < 0.05)
        significant_target_genes <- target_beam[target_beam$qval < 0.05, "gene_short_name"]
        
        if (length(significant_target_genes) > 0) {
          log_info(sprintf("Found %d significantly changed target genes (qval < 0.05)", length(significant_target_genes)))
          
          # Check if these genes exist in mycds data (key fix)
          all_genes_in_mycds <- rownames(mycds)
          existing_significant_genes <- significant_target_genes[significant_target_genes %in% all_genes_in_mycds]
          
          if (length(existing_significant_genes) > 0) {
            log_info(sprintf("Among them, %d genes exist in mycds data", length(existing_significant_genes)))
            
            # Create branch heatmap for significant target genes (reference pseudotime heatmap format)
            if (length(existing_significant_genes) <= 50) {  # Avoid too many genes causing unclear plots
              png(file.path(output_dir, "target_genes_branch_heatmap.png"), 
                  width = 12, height = 8, units = "in", res = 300)
              plot_pseudotime_heatmap(mycds[existing_significant_genes,],
                                     num_clusters = 6,
                                     cores = 1,
                                     show_rownames = T)
              dev.off()
            
              # If gene count is moderate, also create individual gene trajectory plots
              if (length(existing_significant_genes) <= 12) {
                p_branch_trajectory <- plot_genes_branched_pseudotime(mycds[existing_significant_genes, ],
                                                                      branch_point = current_branch_point,
                                                                      color_by = "celltype",
                                                                      ncol = 3)
                ggsave(file.path(output_dir, "target_genes_branch_trajectory.png"), 
                       plot = p_branch_trajectory, width = 15, height = 10)
              }
            } else {
              # If too many genes, only show top 50 most significant
              top_genes <- head(existing_significant_genes, 50)
              log_info(sprintf("Too many genes, only showing top 50 most significant genes"))
              
              png(file.path(output_dir, "target_genes_branch_heatmap_top50.png"), 
                  width = 12, height = 8, units = "in", res = 300)
              plot_pseudotime_heatmap(mycds[top_genes,],
                                     num_clusters = 6,
                                     cores = 1,
                                     show_rownames = T)
              dev.off()
            }
          } else {
            log_info("No significantly changed target genes exist in mycds data")
          }
          
            log_info("Specific gene branch analysis visualization completed")
        } else {
            log_info("No significantly changed target genes found (qval < 0.05)")
        }
      } else {
          log_info("No target genes found in BEAM results")
      }
    }
    
      log_info("Branch analysis completed")
  }
  
    # Save Monocle object and results
  saveRDS(mycds, file.path(output_dir, "monocle_cds_final.rds"))
  
    # Save pseudotime analysis results
  pse_analysis <- pData(mycds)
  write.csv(pse_analysis, file.path(output_dir, "pseudotime_analysis.csv"), row.names = FALSE)
  
    log_info("Monocle2 pseudotime analysis completed")
  
  return(mycds)
}
  

# ============================================================================
# Gene Family Specialized Analysis Mode (scLncR Strategy)
# ============================================================================

# Similar to scLncR's cut_cell_matrix function: extract only target genes to build count matrix
cut_target_genes_matrix <- function(data_matrix, target_genes, project_name) {
  # Find target genes in matrix
  existing_genes <- target_genes[target_genes %in% rownames(data_matrix)]
  if (length(existing_genes) == 0) {
    stop("No target genes found in data")
  }
  
  # Keep only target gene rows, keep all cell columns
  data_cut <- data_matrix[existing_genes, ]
  log_info(sprintf("Extracted %d target genes, keeping all %d cells",
                   nrow(data_cut), ncol(data_cut)))
  
  # Create Seurat object based on filtered matrix
  seu_obj <- CreateSeuratObject(counts = data_cut, project = project_name)
  return(seu_obj)
}

# Removed run_family_only_analysis function, using unified workflow

# ============================================================================
# Main Analysis Workflow
# ============================================================================

run_complete_analysis <- function(config, output_dir) {
  
  # PanFamily single cell RNA analysis workflow
  log_info("Starting PanFamily single cell RNA analysis workflow")
  
  # Check checkpoint, determine starting step
  checkpoint <- detect_checkpoint(output_dir)
  
  # Determine starting step based on checkpoint
  if (checkpoint$step == "branch_completed") {
    # Branch analysis completed, can proceed with subsequent analysis or end
    log_info("Branch analysis completed, continuing from analysis results")
    mycds <- readRDS(checkpoint$file)
    
    # Load BEAM results
    beam_file <- file.path(output_dir, "BEAM_results.rds")
    if (file.exists(beam_file)) {
      BEAM_res <- readRDS(beam_file)
      log_info("Successfully loaded BEAM analysis results")
    }
    
    # Try to load corresponding seurat_obj
    seurat_file <- file.path(output_dir, "seurat_after_branch.rds")
    if (file.exists(seurat_file)) {
      seurat_obj <- readRDS(seurat_file)
      log_info("Successfully loaded seurat_obj after branch analysis")
    } else {
      # Fallback: try to load from other seurat files
      seurat_file <- file.path(output_dir, "seurat_for_branch.rds")
      if (file.exists(seurat_file)) {
        seurat_obj <- readRDS(seurat_file)
        log_info("Successfully loaded seurat_obj from seurat_for_branch.rds")
      } else {
        annotated_file <- file.path(output_dir, "annotated_seurat_object.rds")
        if (file.exists(annotated_file)) {
          seurat_obj <- readRDS(annotated_file)
          log_info("Successfully loaded seurat_obj from annotated_seurat_object.rds")
        } else {
          seurat_obj <- NULL
          log_info("No seurat_obj file found, setting to NULL")
        }
      }
    }
    
    log_info("Branch analysis completed, can proceed with further analysis or visualization")
    
  } else if (checkpoint$step == "branch_analysis") {
    # Start from branch analysis
    log_info("Continuing from branch analysis step")
    mycds <- readRDS(checkpoint$file)
    
    # Try to load corresponding seurat_obj
    seurat_file <- file.path(output_dir, "seurat_for_branch.rds")
    if (file.exists(seurat_file)) {
      seurat_obj <- readRDS(seurat_file)
      log_info("Successfully loaded seurat_obj")
    } else {
      # Fallback: try to load from annotated_seurat_object.rds if no dedicated file
      annotated_file <- file.path(output_dir, "annotated_seurat_object.rds")
      if (file.exists(annotated_file)) {
        seurat_obj <- readRDS(annotated_file)
        log_info("Successfully loaded seurat_obj from annotated_seurat_object.rds")
      } else {
        seurat_obj <- NULL
        log_info("No seurat_obj file found, setting to NULL")
      }
    }
    
    # Only run branch analysis part
    if (get_nested_config(config, "pseudotime.branch_analysis.enabled", default = FALSE)) {
      log_info("Starting branch analysis...")
      
      # Get current configuration parameters for branch analysis
      current_branch_point <- get_nested_config(config, "pseudotime.monocle2_params.branch_point", default = 1)
      current_cores <- get_nested_config(config, "resources.cores", default = 4)
      
      # Check if interactive adjustment of branch analysis parameters is needed
      execution_mode <- get_nested_config(config, "execution.mode", default = "interactive")
      branch_interactive <- get_nested_config(config, "execution.interactive_settings.steps.branch_analysis", default = FALSE)
      
      if (execution_mode == "interactive" && branch_interactive) {
        log_info("Branch analysis parameters need user interactive adjustment")
        
        cat("\n=== Interactive Prompt: Branch Analysis Parameter Settings ===\n")
        cat("Current configuration parameters:\n")
        cat("branch_point: ", current_branch_point, "\n", sep = "")
        cat("cores: ", current_cores, "\n", sep = "")
        cat("progenitor_method: duplicate (fixed)\n")
        
        # User choice with timeout (120 seconds)
        choice <- get_user_input_with_timeout("Do you need to modify branch analysis parameters? (Y/N, will use current parameters if no response in 120 seconds): ", 120)
        
        if (toupper(choice) == "Y") {
          cat("\n=== Branch Analysis Parameter Modification ===\n")
          
          # Get new branch point parameter (with timeout)
          input_branch <- get_user_input_with_timeout(sprintf("Enter branch point [current:%d, will keep current value if no response in 120 seconds]: ", current_branch_point), 120)
          if (input_branch != "TIMEOUT" && input_branch != "") {
            new_branch_point <- as.integer(input_branch)
            if (!is.na(new_branch_point) && new_branch_point > 0) {
              current_branch_point <- new_branch_point
            }
          }
          
          # Get new core count parameter (with timeout)
          input_cores <- get_user_input_with_timeout(sprintf("Enter number of cores to use [current:%d, will keep current value if no response in 120 seconds]: ", current_cores), 120)
          if (input_cores != "TIMEOUT" && input_cores != "") {
            new_cores <- as.integer(input_cores)
            if (!is.na(new_cores) && new_cores > 0 && new_cores <= 32) {
              current_cores <- new_cores
            }
          }
          
          cat("\n=== Updated Parameters ===\n")
          cat("branch_point:", current_branch_point, "\n")
          cat("cores:", current_cores, "\n")
          log_info("Branch analysis parameters updated")
          
        } else {
          log_info("Continuing with current branch analysis parameters")
        }
      }
      
      # Run BEAM analysis
      log_info(sprintf("Running BEAM analysis with parameters: branch_point=%d, cores=%d", current_branch_point, current_cores))
      BEAM_res <- BEAM(mycds, 
                       branch_point = current_branch_point, 
                       progenitor_method = "duplicate",
                       cores = current_cores)
      BEAM_res <- BEAM_res[order(BEAM_res$qval),]
      BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
      
      # Save branch genes
      write.csv(BEAM_res, file.path(output_dir, "branch_genes.csv"), row.names = FALSE)
      
      # Save branch analysis results and working environment
      log_info("Saving branch analysis results and working environment...")
      saveRDS(BEAM_res, file.path(output_dir, "BEAM_results.rds"))
      saveRDS(mycds, file.path(output_dir, "monocle_cds_after_branch.rds"))
      if (exists("seurat_obj") && !is.null(seurat_obj)) {
        saveRDS(seurat_obj, file.path(output_dir, "seurat_after_branch.rds"))
      }
      
      # Visualize branch genes
      plot_genes_branched_heatmap(mycds[row.names(subset(BEAM_res, qval < 1e-4)),],
                                  branch_point = current_branch_point,
                                  num_clusters = 4,
                                  cores = 1,
                                  use_gene_short_name = T,
                                  show_rownames = T)
      ggsave(file.path(output_dir, "branch_heatmap.png"), width = 10, height = 8)
      
      # Specific gene branch analysis
      if (get_nested_config(config, "target_gene_analysis.enabled", default = FALSE) && 
          !is.null(get_nested_config(config, "target_gene_analysis.gene_list_file", default = NULL))) {
        
        log_info("Starting specific gene branch analysis...")
        
        # Read target genes
        target_genes <- readLines(get_nested_config(config, "target_gene_analysis.gene_list_file"))
        target_genes <- target_genes[target_genes != ""]  # Remove empty lines
        log_info(sprintf("Read %d target genes", length(target_genes)))
        
        # Debug information: show basic info about BEAM results and target genes
        log_info(sprintf("BEAM results contain %d genes", nrow(BEAM_res)))
        log_info(sprintf("BEAM result gene name examples: %s", paste(head(BEAM_res$gene_short_name, 3), collapse=", ")))
        log_info(sprintf("Target gene examples: %s", paste(head(target_genes, 3), collapse=", ")))
        
        # Find target genes in BEAM results
        target_beam <- BEAM_res[BEAM_res$gene_short_name %in% target_genes, ]
        
        if (nrow(target_beam) > 0) {
          log_info(sprintf("Found %d target genes in BEAM results", nrow(target_beam)))
          
          # Save target gene branch analysis results
          write.csv(target_beam, file.path(output_dir, "target_genes_branch_analysis.csv"), row.names = FALSE)
          
        # Select only significant target genes (qval < 0.05)
        significant_target_genes <- target_beam[target_beam$qval < 0.05, "gene_short_name"]
        
        if (length(significant_target_genes) > 0) {
          log_info(sprintf("Found %d significantly changed target genes (qval < 0.05)", length(significant_target_genes)))
          
          # Check if these genes exist in mycds data (key fix)
          all_genes_in_mycds <- rownames(mycds)
          existing_significant_genes <- significant_target_genes[significant_target_genes %in% all_genes_in_mycds]
          
          if (length(existing_significant_genes) > 0) {
            log_info(sprintf("Among them, %d genes exist in mycds data", length(existing_significant_genes)))
            
            # Create branch heatmap for significant target genes (reference pseudotime heatmap format)
            if (length(existing_significant_genes) <= 50) {  # Avoid too many genes causing unclear plots
              png(file.path(output_dir, "target_genes_branch_heatmap.png"), 
                  width = 12, height = 8, units = "in", res = 300)
              plot_pseudotime_heatmap(mycds[existing_significant_genes,],
                                     num_clusters = 6,
                                     cores = 1,
                                     show_rownames = T)
              dev.off()
              
              # If gene count is moderate, also create individual gene trajectory plots
              if (length(existing_significant_genes) <= 12) {
                p_branch_trajectory <- plot_genes_branched_pseudotime(mycds[existing_significant_genes, ],
                                                                      branch_point = current_branch_point,
                                                                      color_by = "celltype",
                                                                      ncol = 3)
                ggsave(file.path(output_dir, "target_genes_branch_trajectory.png"), 
                       plot = p_branch_trajectory, width = 15, height = 10)
              }
            } else {
              # If too many genes, only show top 50 most significant
              top_genes <- head(existing_significant_genes, 50)
              log_info(sprintf("Too many genes, only showing top 50 most significant genes"))
              
              png(file.path(output_dir, "target_genes_branch_heatmap_top50.png"), 
                  width = 12, height = 8, units = "in", res = 300)
              plot_pseudotime_heatmap(mycds[top_genes,],
                                     num_clusters = 6,
                                     cores = 1,
                                     show_rownames = T)
              dev.off()
            }
          } else {
            log_info("No significantly changed target genes exist in mycds data")
          }
            
            log_info("Specific gene branch analysis visualization completed")
          } else {
            log_info("No significantly changed target genes found (qval < 0.05)")
          }
        } else {
          log_info("No target genes found in BEAM results")
        }
      }
      
      log_info("Branch analysis completed")
    }
    
  } else if (checkpoint$step == "pseudotime") {
    # Start from pseudotime analysis
    log_info("Continuing from pseudotime analysis step")
    seurat_obj <- readRDS(checkpoint$file)
    
    # Step 8: Pseudotime analysis
    if (get_nested_config(config, "pseudotime.enabled", default = FALSE)) {
      mycds <- run_pseudotime_analysis(seurat_obj, config, output_dir)
    }
    
  } else if (checkpoint$step == "annotation") {
    # Continue after cell annotation
    log_info("Continuing from target gene analysis step")
    seurat_obj <- readRDS(checkpoint$file)
    
    # Step 7: Target gene family expression analysis
    if (get_nested_config(config, "target_gene_analysis.enabled", default = FALSE)) {
      seurat_obj <- analyze_target_genes(seurat_obj, config, output_dir)
    }
    
    # Step 8: Pseudotime analysis
    if (get_nested_config(config, "pseudotime.enabled", default = FALSE)) {
      mycds <- run_pseudotime_analysis(seurat_obj, config, output_dir)
    }
    
  } else {
    # Start complete analysis from beginning
    log_info("Starting complete analysis workflow from beginning")
    
    # Step 0: CellRanger workflow (optional)
    cellranger_result <- run_cellranger_pipeline(config, output_dir)
    if (!is.null(cellranger_result) && !is.null(cellranger_result$config)) {
      # If CellRanger was run, update configuration
      config <- cellranger_result$config
    }
    
    # Step 1: Data loading
    log_info("Single sample mode: loading single sample")
    seurat_obj <- load_sc_data(config, output_dir)
    
    # Step 2: Quality control
    seurat_obj <- run_quality_control(seurat_obj, config, output_dir)
    
    # Step 3: Normalization
    seurat_obj <- run_normalization(seurat_obj, config, output_dir)
    
    # Step 4: Dimension reduction and clustering
    seurat_obj <- run_dimension_reduction(seurat_obj, config, output_dir)
    
    # Step 5: Marker gene identification
    markers <- find_marker_genes(seurat_obj, config, output_dir)
    
    # Step 6: Cell annotation
    seurat_obj <- run_cell_annotation(seurat_obj, markers, config, output_dir)
    
    # Step 7: Target gene family expression analysis
    if (get_nested_config(config, "target_gene_analysis.enabled", default = FALSE)) {
      seurat_obj <- analyze_target_genes(seurat_obj, config, output_dir)
    }
    
    # Step 8: Pseudotime analysis
    if (get_nested_config(config, "pseudotime.enabled", default = FALSE)) {
      mycds <- run_pseudotime_analysis(seurat_obj, config, output_dir)
    }
  }
  
  log_info("PanFamily single cell analysis completed")
  
  # Safe return: ensure seurat_obj is defined
  if (exists("seurat_obj") && !is.null(seurat_obj)) {
    return(seurat_obj)
  } else {
    log_info("seurat_obj undefined or NULL, returning NULL")
    return(NULL)
  }
}

# ============================================================================
# Main Function
# ============================================================================

main <- function() {
  # Set working directory
  setwd(opt$output)
  
  # Load configuration
  config <- load_config(opt$config)
  
  # Ensure configuration parameters are complete, use defaults for missing values
  config <- ensure_config_complete(config)
  
  # Create output directory
  if (!dir.exists(opt$output)) {
    dir.create(opt$output, recursive = TRUE)
  }
  
  log_info("Starting PanFamily single cell RNA analysis")
  log_info(sprintf("Output directory: %s", opt$output))
  
  # Run unified analysis workflow
  result <- run_complete_analysis(config, opt$output)
  
  log_info("Analysis completed!")
}

# Error handling
tryCatch({
  main()
}, error = function(e) {
  log_error(sprintf("Analysis failed: %s", e$message))
  quit(status = 1)
})