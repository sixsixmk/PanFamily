#!/bin/bash

set -euo pipefail


log_info() {
    echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2
}

log_warn() {
    echo "[WARN] $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2
}

log_error() {
    echo "[ERROR] $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2
}

log_step() {
    echo "[STEP] $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2
}

log_sample() {
    echo "[SAMPLE] $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2
}


# Print usage information
print_usage() {
    cat << EOF
USAGE:
    $0 -c <CONFIG_FILE>

OPTIONS:
    -c <FILE>               YAML configuration file
    -h, --help              Show help information
    -v, --version           Show version information

YAML Configuration File Format:
    # Batch analysis configuration file
    batch_analysis:
      data_dir: "/path/to/fastq/directory"      # Sequencing data directory
      output_dir: "/path/to/output"             # Result output directory
      data_type: "PE"                           # PE (paired-end) or SE (single-end)
      single_script: "./singlesample_rnaseq.sh" # Single sample processing script path

    # Reference genome and annotation
    reference:
      genome_file: "/path/to/genome.fa"
      annotation_file: "/path/to/annotation.gtf"

    # Analysis parameters (passed to single sample script)
    analysis:
      threads: 20
      keep_intermediate: false

    # Data processing parameters
    data_processing:
      qc_before: true
      qc_after: true
      trimming: true
      trim_quality: 25
      trim_length: 35
      fastqc_options: ""

Output Files:
    - RNA.txt: Sample list (one sample name per line)
    - [SampleName]/: Independent result directory for each sample
      · featurecounts/: featureCounts quantification results
      · hisat2/: Alignment results
      · data_processing/: Data processing results
    - featureCounts/: Summarized featureCounts results
      · [SampleName].counts: Raw count files
      · [SampleName]_cut.counts: Extracted Geneid and count
    - gene_count_matrix.tsv: Merged gene expression matrix (Geneid × Samples)
    - batch_analysis.log: Batch analysis log
    - batch_analysis_report.txt: Batch analysis report

EXAMPLE:
    $0 -c batch_config.yaml

EOF
}

# Parse YAML configuration file
parse_yaml() {
    local yaml_file="$1"
    if [[ ! -f "$yaml_file" ]]; then
        log_error "YAML configuration file does not exist: $yaml_file"
        exit 1
    fi

    log_info "Parsing YAML configuration file: $yaml_file"

    # Use simple grep to extract parameters, avoid complex regular expressions
    DATA_DIR=$(grep -E "^\s*data_dir:" "$yaml_file" | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')
    OUTPUT_DIR=$(grep -E "^\s*output_dir:" "$yaml_file" | head -1 | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')
    DATA_TYPE=$(grep -E "^\s*data_type:" "$yaml_file" | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')
    GENOME_FILE=$(grep -E "^\s*genome_file:" "$yaml_file" | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')
    ANNOTATION_FILE=$(grep -E "^\s*annotation_file:" "$yaml_file" | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')
    THREADS=$(grep -E "^\s*threads:" "$yaml_file" | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')
    KEEP_INTERMEDIATE=$(grep -E "^\s*keep_intermediate:" "$yaml_file" | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')
    QC_BEFORE=$(grep -E "^\s*qc_before:" "$yaml_file" | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')
    QC_AFTER=$(grep -E "^\s*qc_after:" "$yaml_file" | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')
    TRIMMING=$(grep -E "^\s*trimming:" "$yaml_file" | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')
    TRIM_QUALITY=$(grep -E "^\s*trim_quality:" "$yaml_file" | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')
    TRIM_LENGTH=$(grep -E "^\s*trim_length:" "$yaml_file" | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')
    FASTQC_OPTIONS=$(grep -E "^\s*fastqc_options:" "$yaml_file" | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')

    # DEG analysis parameters
    DEG_ENABLED=$(grep -A 10 "deg_analysis:" "$yaml_file" | grep -E "^\s*enabled:" | head -1 | sed 's/#.*//' | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*$//' || true)
    DEG_GROUP_FILE=$(grep -A 10 "deg_analysis:" "$yaml_file" | grep -E "^\s*group_file:" | head -1 | sed 's/#.*//' | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*$//' || true)
    DEG_MATRIX_TYPE=$(grep -A 10 "deg_analysis:" "$yaml_file" | grep -E "^\s*matrix_type:" | head -1 | sed 's/#.*//' | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*$//' || true)
    DEG_PVALUE_CUTOFF=$(grep -A 10 "deg_analysis:" "$yaml_file" | grep -E "^\s*pvalue_cutoff:" | head -1 | sed 's/#.*//' | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*$//' || true)
    DEG_LOG2FC_CUTOFF=$(grep -A 10 "deg_analysis:" "$yaml_file" | grep -E "^\s*log2fc_cutoff:" | head -1 | sed 's/#.*//' | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*$//' || true)

    # Gene family parameters
    GENEFAMILY_LIST=$(grep -A 5 "gene_family:" "$yaml_file" | grep -E "^\s*genefamily_list:" | head -1 | sed 's/#.*//' | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*$//' || true)
    GENEFAMILY_PREFIX=$(grep -A 5 "gene_family:" "$yaml_file" | grep -E "^\s*genefamily_prefix:" | head -1 | sed 's/#.*//' | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*$//' || true)
    GENEFAMILY_EXTRACT_DEG=$(grep -A 5 "gene_family:" "$yaml_file" | grep -E "^\s*extract_from_deg:" | head -1 | sed 's/#.*//' | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*$//' || true)

    # Set default values
    DATA_DIR=${DATA_DIR:-""}
    OUTPUT_DIR=${OUTPUT_DIR:-"./results"}
    DATA_TYPE=${DATA_TYPE:-PE}
    THREADS=${THREADS:-8}
    KEEP_INTERMEDIATE=${KEEP_INTERMEDIATE:-false}
    QC_BEFORE=${QC_BEFORE:-true}
    QC_AFTER=${QC_AFTER:-true}
    TRIMMING=${TRIMMING:-true}
    TRIM_QUALITY=${TRIM_QUALITY:-25}
    TRIM_LENGTH=${TRIM_LENGTH:-35}
    FASTQC_OPTIONS=${FASTQC_OPTIONS:-""}

    # DEG analysis default values
    DEG_ENABLED=${DEG_ENABLED:-false}
    DEG_GROUP_FILE=${DEG_GROUP_FILE:-""}
    DEG_MATRIX_TYPE=${DEG_MATRIX_TYPE:-"count"}
    DEG_PVALUE_CUTOFF=${DEG_PVALUE_CUTOFF:-0.05}
    DEG_LOG2FC_CUTOFF=${DEG_LOG2FC_CUTOFF:-1}

    # Gene family default values
    GENEFAMILY_LIST=${GENEFAMILY_LIST:-""}
    GENEFAMILY_PREFIX=${GENEFAMILY_PREFIX:-""}
    GENEFAMILY_EXTRACT_DEG=${GENEFAMILY_EXTRACT_DEG:-false}

    log_info "YAML configuration parsing completed"
}

# Identify sequencing data files
identify_samples() {
    log_step "Identifying sequencing data samples"

    if [[ ! -d "$DATA_DIR" ]]; then
        log_error "Data directory does not exist: $DATA_DIR"
        exit 1
    fi

    cd "$DATA_DIR" || exit 1

    local samples=()

    if [[ "$DATA_TYPE" == "PE" ]]; then
        # Paired-end data: identify _1 and _2 pairs
        log_info "Identifying paired-end sequencing data..."

        for file in *_1.fastq.gz *_1.fq.gz *_1.fastq *_1.fq; do
            [[ ! -f "$file" ]] && continue

            # Extract sample name (remove _1 and suffix)
            local sample=$(basename "$file" | sed -E 's/_1\.(fastq|fq)(\.gz)?$//')

            # Check if corresponding _2 file exists
            local found_pair=false
            for ext in "fastq.gz" "fq.gz" "fastq" "fq"; do
                if [[ -f "${sample}_2.${ext}" ]]; then
                    found_pair=true
                    break
                fi
            done

            if [[ "$found_pair" == true ]]; then
                samples+=("$sample")
                log_info "  Found sample: $sample (paired-end)"
            else
                log_warn "  Sample $sample missing paired file, skipping"
            fi
        done

    else
        # Single-end data
        log_info "Identifying single-end sequencing data..."

        for file in *.fastq.gz *.fq.gz *.fastq *.fq; do
            [[ ! -f "$file" ]] && continue

            # Skip paired-end marked files
            [[ "$file" =~ _[12]\.(fastq|fq) ]] && continue

            # Extract sample name (remove suffix)
            local sample=$(basename "$file" | sed -E 's/\.(fastq|fq)(\.gz)?$//')

            samples+=("$sample")
            log_info "  Found sample: $sample (single-end)"
        done
    fi

    cd - > /dev/null

    if [[ ${#samples[@]} -eq 0 ]]; then
        log_error "No samples found in directory $DATA_DIR"
        log_error "Supported formats: .fastq.gz, .fq.gz, .fastq, .fq"
        exit 1
    fi

    log_info "Total ${#samples[@]} samples identified"

    # Generate RNA.txt
    local rna_txt="$OUTPUT_DIR/RNA.txt"
    > "$rna_txt"

    for sample in "${samples[@]}"; do
        echo "$sample" >> "$rna_txt"
    done

    log_info "Sample list saved to: $rna_txt"

    SAMPLES=("${samples[@]}")
}

# Generate temporary configuration file for single sample
generate_sample_config() {
    local sample="$1"
    local config_file="$OUTPUT_DIR/.tmp_${sample}_config.yaml"

    # Determine input file paths
    local read1=""
    local read2=""

    if [[ "$DATA_TYPE" == "PE" ]]; then
        # Find paired files
            for ext in "fastq.gz" "fq.gz" "fastq" "fq"; do
            if [[ -f "$DATA_DIR/${sample}_1.${ext}" ]]; then
                read1="$DATA_DIR/${sample}_1.${ext}"
                read2="$DATA_DIR/${sample}_2.${ext}"
                break
            fi
        done
    else
        # Find single-end file
        for ext in "fastq.gz" "fq.gz" "fastq" "fq"; do
            if [[ -f "$DATA_DIR/${sample}.${ext}" ]]; then
                read1="$DATA_DIR/${sample}.${ext}"
                break
            fi
        done
    fi

    if [[ -z "$read1" ]]; then
        log_error "Cannot find sequencing files for sample $sample"
        return 1
    fi

    # Generate configuration file
    cat > "$config_file" << EOF
# Auto-generated single sample configuration file
# Sample: $sample

analysis:
  sample_name: "$sample"
  threads: $THREADS
  keep_intermediate: $KEEP_INTERMEDIATE

input:
  read1: "$read1"
EOF

    if [[ "$DATA_TYPE" == "PE" ]]; then
        echo "  read2: \"$read2\"" >> "$config_file"
    fi

    cat >> "$config_file" << EOF
  genome_file: "$GENOME_FILE"
  annotation_file: "$ANNOTATION_FILE"

output:
  output_dir: "$OUTPUT_DIR/$sample"

data_processing:
  qc_before: $QC_BEFORE
  qc_after: $QC_AFTER
  trimming: $TRIMMING
  trim_quality: $TRIM_QUALITY
  trim_length: $TRIM_LENGTH
  fastqc_options: "$FASTQC_OPTIONS"
EOF

    echo "$config_file"
}

# Process single sample
process_sample() {
    local sample="$1"
    local sample_num="$2"
    local total="$3"

    log_sample "Start processing sample [$sample_num/$total]: $sample"

    # Generate temporary configuration file
    local config_file=$(generate_sample_config "$sample")

    if [[ -z "$config_file" ]]; then
        log_error "Configuration file generation failed for sample $sample, skipping"
        return 1
    fi

    # Call single sample processing script
    log_info "Calling single sample processing script: $SINGLE_SCRIPT"

    if bash "$SINGLE_SCRIPT" -c "$config_file"; then
        log_info "Sample $sample processing completed"
        rm -f "$config_file"
        return 0
    else
        log_error "Sample $sample processing failed"
        return 1
    fi
}

# Batch process all samples
batch_process() {
    log_step "Starting batch processing of ${#SAMPLES[@]} samples"

    local success_count=0
    local failed_count=0
    local failed_samples=()

    local total=${#SAMPLES[@]}
    local current=0

    for sample in "${SAMPLES[@]}"; do
        current=$((current + 1))

        echo ""
        echo "========================================================================"

        if process_sample "$sample" "$current" "$total"; then
            success_count=$((success_count + 1))
        else
            failed_count=$((failed_count + 1))
            failed_samples+=("$sample")
        fi

        echo "========================================================================"
        echo ""
    done

    # Output summary information
    log_step "Batch processing completed"
    log_info "Total samples: $total"
    log_info "Successful: $success_count"
    log_info "Failed: $failed_count"

    if [[ $failed_count -gt 0 ]]; then
        log_warn "Failed samples:"
        for sample in "${failed_samples[@]}"; do
            log_warn "  - $sample"
        done
    fi
}

# Collect and process featureCounts results
collect_featurecounts_results() {
    log_step "Collecting featureCounts results"

    local fc_dir="$OUTPUT_DIR/featureCounts"
    mkdir -p "$fc_dir"

    local collected_count=0
    local failed_collect=()

    for sample in "${SAMPLES[@]}"; do
        local sample_fc_file="$OUTPUT_DIR/$sample/featurecounts/${sample}.counts"

        if [[ -f "$sample_fc_file" ]]; then
            log_info "Collecting sample: $sample"
            cp "$sample_fc_file" "$fc_dir/"

            # Extract column 1 (Geneid) and column 7 (count), remove comment lines
            log_info "Processing sample: $sample (extracting Geneid and count columns)"
            cut -f 1,7 "$fc_dir/${sample}.counts" | grep -v "^#" > "$fc_dir/${sample}_cut.counts"

            ((collected_count++))
        else
            log_warn "featureCounts results not found for sample $sample: $sample_fc_file"
            failed_collect+=("$sample")
        fi
    done

    log_info "Successfully collected results from $collected_count samples"

    if [[ ${#failed_collect[@]} -gt 0 ]]; then
        log_warn "The following samples have no featureCounts results:"
        for sample in "${failed_collect[@]}"; do
            log_warn "  - $sample"
        done
    fi

    # Return number of collected samples
    echo "$collected_count"
}

# Merge featureCounts results into expression matrix
merge_count_matrix() {
    log_step "Merging gene expression matrix"

    # Extract strain name (from OUTPUT_DIR basename)
    local strain_name=$(basename "$OUTPUT_DIR")

    local fc_dir="$OUTPUT_DIR/featureCounts"
    local matrix_file="$OUTPUT_DIR/${strain_name}_count_matrix.tsv"

    # Check if processed files exist
    local cut_files=("$fc_dir"/*_cut.counts)
    if [[ ! -f "${cut_files[0]}" ]]; then
        log_error "No processed count files found"
        return 1
    fi

    # Run Python script to merge expression matrix
    log_info "Using Python script to merge expression matrix: $MERGE_SCRIPT"

    if command -v python3 &> /dev/null; then
        python3 "$MERGE_SCRIPT" "$fc_dir" "$matrix_file"

        if [[ -f "$matrix_file" ]]; then
            log_info "Expression matrix merging completed: $matrix_file"

            # Display matrix information
            local gene_count=$(tail -n +2 "$matrix_file" | wc -l)
            local sample_count=$(head -1 "$matrix_file" | awk '{print NF-1}')
            log_info "Matrix dimensions: $gene_count genes × $sample_count samples"

            return 0
        else
            log_error "Expression matrix generation failed"
            return 1
        fi
    else
        log_error "python3 not found, cannot merge expression matrix"
        log_error "Please install Python 3 and pandas: pip install pandas"
        return 1
    fi
}

# Calculate FPKM and TPM
calculate_fpkm_tpm() {
    log_step "Calculating FPKM and TPM"

    # Extract strain name (from OUTPUT_DIR basename)
    local strain_name=$(basename "$OUTPUT_DIR")

    local matrix_file="$OUTPUT_DIR/${strain_name}_count_matrix.tsv"

    # Check if count matrix exists
    if [[ ! -f "$matrix_file" ]]; then
        log_error "Count matrix does not exist: $matrix_file"
        return 1
    fi

    # Check GTF file
    if [[ -z "$ANNOTATION_FILE" ]] || [[ ! -f "$ANNOTATION_FILE" ]]; then
        log_warn "GTF annotation file does not exist, skipping FPKM/TPM calculation"
        return 0
    fi

    log_info "Using Python script to calculate FPKM/TPM: $FPKM_SCRIPT"
    log_info "Using GTF annotation file: $ANNOTATION_FILE"

    # Run Python script
    if command -v python3 &> /dev/null; then
        local output_prefix="$OUTPUT_DIR/${strain_name}"

        python3 "$FPKM_SCRIPT" \
            -c "$matrix_file" \
            -g "$ANNOTATION_FILE" \
            -o "$output_prefix"

        if [[ -f "${output_prefix}_fpkm_matrix.tsv" ]] && [[ -f "${output_prefix}_tpm_matrix.tsv" ]]; then
            log_info "FPKM/TPM calculation completed"
            log_info "  - FPKM: ${output_prefix}_fpkm_matrix.tsv"
            log_info "  - TPM: ${output_prefix}_tpm_matrix.tsv"
            return 0
        else
            log_warn "FPKM/TPM calculation may have failed, please check logs"
            return 1
        fi
    else
        log_warn "python3 not found, skipping FPKM/TPM calculation"
        return 0
    fi
}

# Differential expression gene analysis
run_deg_analysis() {
    # Check if DEG analysis is enabled
    if [[ "$DEG_ENABLED" != "true" ]] && [[ "$DEG_ENABLED" != "TRUE" ]] && [[ "$DEG_ENABLED" != "T" ]]; then
        # Not enabled, skip silently
        return 0
    fi

    log_step "Differential Expression Gene (DEG) Analysis"

    # Validate grouping file exists
    if [[ -z "$DEG_GROUP_FILE" ]] || [[ ! -f "$DEG_GROUP_FILE" ]]; then
        log_warn "DEG analysis is enabled, but grouping file does not exist or not specified: $DEG_GROUP_FILE"
        log_warn "Skipping differential analysis"
        return 0
    fi

    # Extract strain name (from OUTPUT_DIR basename)
    local strain_name=$(basename "$OUTPUT_DIR")

    # Determine which matrix file to use
    local input_matrix=""
    case "$DEG_MATRIX_TYPE" in
        count)
            input_matrix="$OUTPUT_DIR/${strain_name}_count_matrix.tsv"
            ;;
        fpkm)
            input_matrix="$OUTPUT_DIR/${strain_name}_fpkm_matrix.tsv"
            ;;
        tpm)
            input_matrix="$OUTPUT_DIR/${strain_name}_tpm_matrix.tsv"
            ;;
        *)
            log_error "Unsupported matrix type: $DEG_MATRIX_TYPE (supported: count, fpkm, tpm)"
            return 1
            ;;
    esac

    # Check if matrix file exists
    if [[ ! -f "$input_matrix" ]]; then
        log_warn "Matrix file does not exist: $input_matrix"
        log_warn "Skipping differential analysis"
        return 0
    fi

    log_info "Using matrix: $input_matrix ($DEG_MATRIX_TYPE)"
    log_info "Grouping file: $DEG_GROUP_FILE"
    log_info "P-value threshold: $DEG_PVALUE_CUTOFF"
    log_info "Log2FC threshold: $DEG_LOG2FC_CUTOFF"

    # Create DEG results directory
    local deg_dir="$OUTPUT_DIR/DEG_results"
    mkdir -p "$deg_dir"

    # Convert TSV to CSV (DEG.R requires CSV format)
    local csv_matrix="$deg_dir/input_matrix.csv"
    log_info "Converting matrix format TSV -> CSV"

    # Use awk to convert, replace tabs with commas
    awk 'BEGIN {FS="\t"; OFS=","} {$1=$1; print}' "$input_matrix" > "$csv_matrix"

    # Check if R is installed
    if ! command -v Rscript &> /dev/null; then
        log_error "Rscript not found, cannot run differential analysis"
        log_error "Please install R and DESeq2: install.packages('BiocManager'); BiocManager::install('DESeq2')"
        return 1
    fi

    # Auto-detect script directory
    local script_dir=$(cd "$(dirname "$0")" && pwd)
    local deg_script="$script_dir/DEG.R"

    if [[ ! -f "$deg_script" ]]; then
        log_error "DEG.R script does not exist: $deg_script"
        return 1
    fi

    log_info "Running DESeq2 differential analysis..."
    log_info "Using script: $deg_script"

    # Run DEG.R
    if Rscript "$deg_script" \
        "$csv_matrix" \
        "$DEG_GROUP_FILE" \
        "$deg_dir" \
        "$DEG_PVALUE_CUTOFF" \
        "$DEG_LOG2FC_CUTOFF" 2>&1 | while IFS= read -r line; do
            log_info "  $line"
        done; then

        log_info "Differential analysis completed"
        log_info "Results saved in: $deg_dir"

        # Count result files
        local deg_files=$(find "$deg_dir" -name "DEG_*.csv" -not -name "*input_matrix*" | wc -l)
        if [[ $deg_files -gt 0 ]]; then
            log_info "Generated $deg_files DEG result files"
        fi

        return 0
    else
        log_error "Differential analysis failed"
        return 1
    fi
}

# Extract gene family expression matrix
extract_genefamily_matrix() {
    # Check if gene family extraction is configured
    if [[ -z "$GENEFAMILY_LIST" ]] || [[ -z "$GENEFAMILY_PREFIX" ]]; then
        # Not configured, skip silently
        return 0
    fi

    log_step "Extracting gene family expression matrix: $GENEFAMILY_PREFIX"

    # Validate gene list file exists
    if [[ ! -f "$GENEFAMILY_LIST" ]]; then
        log_warn "Gene family list file does not exist: $GENEFAMILY_LIST, skipping extraction"
        return 0
    fi

    # Extract strain name (from OUTPUT_DIR basename)
    local strain_name=$(basename "$OUTPUT_DIR")

    # Create output directory
    local genefamily_dir="$OUTPUT_DIR/${GENEFAMILY_PREFIX}_matrix"
    mkdir -p "$genefamily_dir"

    log_info "Gene family list: $GENEFAMILY_LIST"
    log_info "Output directory: $genefamily_dir"

    # Count genes in the gene list
    local total_genes=$(grep -v "^#" "$GENEFAMILY_LIST" | grep -v "^$" | wc -l)
    log_info "Gene list contains $total_genes gene IDs"

    # Define matrix files to process
    declare -A matrix_files=(
        ["${strain_name}_count_matrix.tsv"]="${GENEFAMILY_PREFIX}_count_matrix.tsv"
        ["${strain_name}_fpkm_matrix.tsv"]="${GENEFAMILY_PREFIX}_fpkm_matrix.tsv"
        ["${strain_name}_tpm_matrix.tsv"]="${GENEFAMILY_PREFIX}_tpm_matrix.tsv"
    )

    local extracted_count=0

    # Process each matrix file for extraction
    for input_file in "${!matrix_files[@]}"; do
        local input_path="$OUTPUT_DIR/$input_file"
        local output_file="${matrix_files[$input_file]}"
        local output_path="$genefamily_dir/$output_file"

        if [[ ! -f "$input_path" ]]; then
            log_warn "Matrix file does not exist: $input_path, skipping"
            continue
        fi

        log_info "Extracting $input_file -> $output_file"

        # Use awk for extraction
        # 1. Read gene list into array
        # 2. Output header
        # 3. For each row, if first column is in gene list, output the row
        awk -v genelist="$GENEFAMILY_LIST" '
        BEGIN {
            # Read gene list
            while ((getline < genelist) > 0) {
                if ($0 !~ /^#/ && $0 !~ /^[[:space:]]*$/) {
                    genes[$1] = 1
                }
            }
            close(genelist)
        }
        NR == 1 {
            # Output header
            print $0
            next
        }
        {
            # Check if first column is in gene list
            if ($1 in genes) {
                print $0
                matched++
            }
        }
        END {
            if (matched > 0) {
                print "# Extracted " matched " genes" > "/dev/stderr"
            } else {
                print "# Warning: No genes matched" > "/dev/stderr"
            }
        }
        ' "$input_path" > "$output_path"

        # Check extraction results
        local extracted_genes=$(tail -n +2 "$output_path" | wc -l)
        if [[ $extracted_genes -gt 0 ]]; then
            log_info "  Successfully extracted $extracted_genes genes"
            ((extracted_count++))
        else
            log_warn "  No genes matched"
        fi
    done

    if [[ $extracted_count -gt 0 ]]; then
        log_info "Gene family matrix extraction completed"
        log_info "Results saved in: $genefamily_dir"
    else
        log_warn "No gene family matrices successfully extracted"
    fi

    # If configured to extract from DEG results
    if [[ "$GENEFAMILY_EXTRACT_DEG" == "true" ]] || [[ "$GENEFAMILY_EXTRACT_DEG" == "TRUE" ]] || [[ "$GENEFAMILY_EXTRACT_DEG" == "T" ]]; then
        local deg_dir="$OUTPUT_DIR/DEG_results"

        if [[ ! -d "$deg_dir" ]]; then
            log_warn "DEG results directory does not exist, skipping gene family member extraction from DEG"
        else
            log_info "Extracting gene family members from DEG results"

            # Find all DEG result files
            local deg_files=($(find "$deg_dir" -name "DEG_*.csv" -not -name "*input_matrix*" -type f))

            if [[ ${#deg_files[@]} -eq 0 ]]; then
                log_warn "No DEG result files found"
            else
                log_info "Found ${#deg_files[@]} DEG result files"

                # Extract from each DEG file
                for deg_file in "${deg_files[@]}"; do
                    local deg_basename=$(basename "$deg_file" .csv)
                    local output_file="$genefamily_dir/${GENEFAMILY_PREFIX}_${deg_basename}.csv"

                    log_info "  Processing: $deg_basename"

                    # Use awk to extract (CSV format)
                    awk -v genelist="$GENEFAMILY_LIST" '
                    BEGIN {
                        FS = ","
                        OFS = ","
                        # Read gene list
                        while ((getline < genelist) > 0) {
                            if ($0 !~ /^#/ && $0 !~ /^[[:space:]]*$/) {
                                genes[$1] = 1
                            }
                        }
                        close(genelist)
                    }
                    NR == 1 {
                        # Output header
                        print $0
                        next
                    }
                    {
                        # CSV file first column may have quotes
                        gene_id = $1
                        gsub(/^"|"$/, "", gene_id)  # Remove quotes
                        if (gene_id in genes) {
                            print $0
                            matched++
                        }
                    }
                    END {
                        if (matched > 0) {
                            print "    Extracted " matched " genes" > "/dev/stderr"
                        }
                    }
                    ' "$deg_file" > "$output_file"

                    local extracted=$(tail -n +2 "$output_file" | wc -l)
                    if [[ $extracted -gt 0 ]]; then
                        log_info "    Extracted $extracted genes"
                    fi
                done

                log_info "Gene family member extraction from DEG results completed"
            fi
        fi
    fi

    return 0
}

# Generate batch analysis report
generate_batch_report() {
    log_step "Generating batch analysis report"

    # Extract strain name (from OUTPUT_DIR basename)
    local strain_name=$(basename "$OUTPUT_DIR")

    local report_file="$OUTPUT_DIR/batch_analysis_report.txt"
    local matrix_file="$OUTPUT_DIR/${strain_name}_count_matrix.tsv"

    cat > "$report_file" << EOF
=======================================
    PanFamily RNAseq Batch Analysis Report
=======================================

Analysis Information:
- Data Type: $DATA_TYPE
- Quantification Method: featureCounts
- Analysis Time: $(date)

Configuration Parameters:
- Data Directory: $DATA_DIR
- Output Directory: $OUTPUT_DIR
- Reference Genome: $GENOME_FILE
- GTF Annotation: $ANNOTATION_FILE
- Threads: $THREADS

Sample List:
EOF

    for sample in "${SAMPLES[@]}"; do
        echo "  - $sample" >> "$report_file"
    done

    cat >> "$report_file" << EOF

Result Files:
- Sample list: $OUTPUT_DIR/RNA.txt
- Sample results: $OUTPUT_DIR/[SampleName]/
- featureCounts summary: $OUTPUT_DIR/featureCounts/
EOF

    if [[ -f "$matrix_file" ]]; then
        local gene_count=$(tail -n +2 "$matrix_file" | wc -l)
        local sample_count=$(head -1 "$matrix_file" | awk '{print NF-1}')
        cat >> "$report_file" << EOF
- Expression matrix: $matrix_file
  · Gene count: $gene_count
  · Sample count: $sample_count
EOF
    fi

    # Add FPKM/TPM information
    if [[ -f "$OUTPUT_DIR/${strain_name}_fpkm_matrix.tsv" ]]; then
        cat >> "$report_file" << EOF
- FPKM matrix: $OUTPUT_DIR/${strain_name}_fpkm_matrix.tsv
EOF
    fi

    if [[ -f "$OUTPUT_DIR/${strain_name}_tpm_matrix.tsv" ]]; then
        cat >> "$report_file" << EOF
- TPM matrix: $OUTPUT_DIR/${strain_name}_tpm_matrix.tsv
EOF
    fi

    # Add DEG analysis results information
    if [[ "$DEG_ENABLED" == "true" ]] || [[ "$DEG_ENABLED" == "TRUE" ]] || [[ "$DEG_ENABLED" == "T" ]]; then
        local deg_dir="$OUTPUT_DIR/DEG_results"
        if [[ -d "$deg_dir" ]]; then
            cat >> "$report_file" << EOF

Differential Expression Analysis Results:
- Matrix type used: $DEG_MATRIX_TYPE
- P-value threshold: $DEG_PVALUE_CUTOFF
- Log2FC threshold: $DEG_LOG2FC_CUTOFF
- Result directory: $deg_dir
EOF

            # Count DEG results
            if [[ -f "$deg_dir/DEG_up.csv" ]] && [[ -f "$deg_dir/DEG_down.csv" ]]; then
                local up_count=$(tail -n +2 "$deg_dir/DEG_up.csv" | wc -l)
                local down_count=$(tail -n +2 "$deg_dir/DEG_down.csv" | wc -l)
                local mid_count=0
                [[ -f "$deg_dir/DEG_intermediate.csv" ]] && mid_count=$(tail -n +2 "$deg_dir/DEG_intermediate.csv" | wc -l)

                cat >> "$report_file" << EOF
- Comparison: Treatment vs Control
  · Up-regulated genes: $up_count
  · Down-regulated genes: $down_count
  · Intermediate expression: $mid_count
EOF
            fi
        fi
    fi

    # Add gene family extraction information
    if [[ -n "$GENEFAMILY_PREFIX" ]] && [[ -d "$OUTPUT_DIR/${GENEFAMILY_PREFIX}_matrix" ]]; then
        cat >> "$report_file" << EOF

Gene Family Extraction Results ($GENEFAMILY_PREFIX):
EOF

        local genefamily_dir="$OUTPUT_DIR/${GENEFAMILY_PREFIX}_matrix"
        if [[ -f "$genefamily_dir/${GENEFAMILY_PREFIX}_count_matrix.tsv" ]]; then
            local gf_gene_count=$(tail -n +2 "$genefamily_dir/${GENEFAMILY_PREFIX}_count_matrix.tsv" | wc -l)
            cat >> "$report_file" << EOF
- Count matrix: $genefamily_dir/${GENEFAMILY_PREFIX}_count_matrix.tsv
  · Extracted gene count: $gf_gene_count
EOF
        fi

        if [[ -f "$genefamily_dir/${GENEFAMILY_PREFIX}_fpkm_matrix.tsv" ]]; then
            cat >> "$report_file" << EOF
- FPKM matrix: $genefamily_dir/${GENEFAMILY_PREFIX}_fpkm_matrix.tsv
EOF
        fi

        if [[ -f "$genefamily_dir/${GENEFAMILY_PREFIX}_tpm_matrix.tsv" ]]; then
            cat >> "$report_file" << EOF
- TPM matrix: $genefamily_dir/${GENEFAMILY_PREFIX}_tpm_matrix.tsv
EOF
        fi
    fi

    cat >> "$report_file" << EOF

=======================================
Batch Analysis Completed!
=======================================
EOF

    log_info "Batch analysis report generated: $report_file"
}

main() {
    log_info "Starting batch transcriptome analysis"

    # Create output directory
    mkdir -p "$OUTPUT_DIR"

    # Set log file
    LOG_FILE="$OUTPUT_DIR/batch_analysis.log"
    exec &> >(tee -a "$LOG_FILE")

    # Auto-detect script directory (use absolute path)
    local script_dir=$(cd "$(dirname "$0")" && pwd)
    log_info "Script directory: $script_dir"

    # Auto-detect single sample processing script singleRna.sh
    SINGLE_SCRIPT="$script_dir/singleRna.sh"
    if [[ ! -f "$SINGLE_SCRIPT" ]]; then
        log_error "Single sample processing script does not exist: $SINGLE_SCRIPT"
        log_error "Please ensure singleRna.sh and RNAseq.sh are in the same directory"
        exit 1
    fi
    log_info "Using single sample processing script: $SINGLE_SCRIPT"

    # Auto-detect merge_counts_fixed.py
    MERGE_SCRIPT="$script_dir/merge_counts_fixed.py"
    if [[ ! -f "$MERGE_SCRIPT" ]]; then
        log_error "Merge counts script does not exist: $MERGE_SCRIPT"
        log_error "Please ensure merge_counts_fixed.py and RNAseq.sh are in the same directory"
        exit 1
    fi
    log_info "Using merge counts script: $MERGE_SCRIPT"

    # Auto-detect calculate_fpkm_tpm.py
    FPKM_SCRIPT="$script_dir/calculate_fpkm_tpm.py"
    if [[ ! -f "$FPKM_SCRIPT" ]]; then
        log_error "FPKM/TPM calculation script does not exist: $FPKM_SCRIPT"
        log_error "Please ensure calculate_fpkm_tpm.py and RNAseq.sh are in the same directory"
        exit 1
    fi
    log_info "Using FPKM/TPM calculation script: $FPKM_SCRIPT"

    # Identify samples
    identify_samples

    # Batch process
    batch_process

    # Collect featureCounts results
    local collected=$(collect_featurecounts_results)

    # Merge expression matrix
    if [[ $collected -gt 0 ]]; then
        merge_count_matrix

        # Calculate FPKM and TPM
        calculate_fpkm_tpm

        # Differential expression analysis (if enabled)
        run_deg_analysis

        # Extract gene family matrix (if configured)
        extract_genefamily_matrix
    else
        log_warn "No featureCounts results collected, skipping merge steps"
    fi

    # Generate report
    generate_batch_report

    echo ""
    log_info "Batch analysis completed! Results saved in: $OUTPUT_DIR"
    log_info "View report: $OUTPUT_DIR/batch_analysis_report.txt"
    log_info "Sample list: $OUTPUT_DIR/RNA.txt"

    if [[ -f "$OUTPUT_DIR/gene_count_matrix.tsv" ]]; then
        log_info "Expression matrix: $OUTPUT_DIR/gene_count_matrix.tsv"
    fi

    if [[ -f "$OUTPUT_DIR/normalized_fpkm.tsv" ]]; then
        log_info "FPKM matrix: $OUTPUT_DIR/normalized_fpkm.tsv"
    fi

    if [[ -f "$OUTPUT_DIR/normalized_tpm.tsv" ]]; then
        log_info "TPM matrix: $OUTPUT_DIR/normalized_tpm.tsv"
    fi

    if [[ "$DEG_ENABLED" == "true" ]] || [[ "$DEG_ENABLED" == "TRUE" ]] || [[ "$DEG_ENABLED" == "T" ]]; then
        local deg_dir="$OUTPUT_DIR/DEG_results"
        if [[ -d "$deg_dir" ]]; then
            echo "" >&2
            log_info "Differential expression analysis results:"
            log_info "  Result directory: $deg_dir"
            if [[ -f "$deg_dir/DEG_all.csv" ]]; then
                log_info "  Analysis completed: Treatment vs Control"
            fi
        fi
    fi

    if [[ -n "$GENEFAMILY_PREFIX" ]] && [[ -d "$OUTPUT_DIR/${GENEFAMILY_PREFIX}_matrix" ]]; then
        echo "" >&2
        log_info "Gene family ($GENEFAMILY_PREFIX) extraction results:"
        local genefamily_dir="$OUTPUT_DIR/${GENEFAMILY_PREFIX}_matrix"
        if [[ -f "$genefamily_dir/${GENEFAMILY_PREFIX}_count_matrix.tsv" ]]; then
            log_info "  Count matrix: $genefamily_dir/${GENEFAMILY_PREFIX}_count_matrix.tsv"
        fi
        if [[ -f "$genefamily_dir/${GENEFAMILY_PREFIX}_fpkm_matrix.tsv" ]]; then
            log_info "  FPKM matrix: $genefamily_dir/${GENEFAMILY_PREFIX}_fpkm_matrix.tsv"
        fi
        if [[ -f "$genefamily_dir/${GENEFAMILY_PREFIX}_tpm_matrix.tsv" ]]; then
            log_info "  TPM matrix: $genefamily_dir/${GENEFAMILY_PREFIX}_tpm_matrix.tsv"
        fi
    fi
    echo ""
    
    # Explicitly return success
    return 0
}

# ============================================================================
# 参数解析
# ============================================================================

# 默认参数
CONFIG_FILE=""
DATA_DIR=""
OUTPUT_DIR=""
DATA_TYPE="PE"
GENOME_FILE=""
ANNOTATION_FILE=""
THREADS=8
KEEP_INTERMEDIATE="false"
QC_BEFORE="true"
QC_AFTER="true"
TRIMMING="true"
TRIM_QUALITY="25"
TRIM_LENGTH="35"
FASTQC_OPTIONS=""
DEG_ENABLED="false"
DEG_GROUP_FILE=""
DEG_MATRIX_TYPE="count"
DEG_PVALUE_CUTOFF="0.05"
DEG_LOG2FC_CUTOFF="1"
GENEFAMILY_LIST=""
GENEFAMILY_PREFIX=""
GENEFAMILY_EXTRACT_DEG="false"
while [[ $# -gt 0 ]]; do
    case $1 in
        -c)
            CONFIG_FILE="$2"
            shift 2
            ;;
        -h|--help)
            print_usage
            exit 0
            ;;
        -v|--version)
            echo "$PROGRAM_NAME"
            exit 0
            ;;
        *)
            log_error "Unknown parameter: $1"
            echo ""
            print_usage
            exit 1
            ;;
    esac
done

# Validate configuration file
if [[ -z "$CONFIG_FILE" ]]; then
    log_error "Missing required parameter: -c <configuration file>"
    echo ""
    print_usage
    exit 1
fi

if [[ ! -f "$CONFIG_FILE" ]]; then
    log_error "Configuration file does not exist: $CONFIG_FILE"
    exit 1
fi

# Parse configuration
parse_yaml "$CONFIG_FILE"

# Validate required parameters
if [[ -z "$DATA_DIR" ]] || [[ -z "$OUTPUT_DIR" ]]; then
    log_error "Configuration file missing required parameters: data_dir or output_dir"
    exit 1
fi

if [[ -z "$GENOME_FILE" ]] || [[ -z "$ANNOTATION_FILE" ]]; then
    log_error "Configuration file missing required parameters: genome_file or annotation_file"
    exit 1
fi

# Run main program
main