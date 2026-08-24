#!/bin/bash
# ============================================================================
# Pan-genome Gene Family Analysis Module
# ============================================================================

set -euo pipefail 

CONFIG_FILE="config.yaml"

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Pan-genome Gene Family Analysis Pipeline

OPTIONS:
    -c, --config FILE    Configuration file (default: config.yaml)
    -h, --help          Show this help message
    -v, --version       Show version information

DESCRIPTION:
    This pipeline performs the following steps:
    1. Create BLAST database from protein sequences
    2. Run BLASTP self-alignment
    3. Cluster genes into families based on identity
    4. Rename sequences with family IDs
    5. Extract representative sequences (longest per family)
    6. Generate final gene family list

EXAMPLE:
    $0 -c my_config.yaml

EOF
}


version() {
    echo "Pan-genome Pipeline v1.0.0"
}


log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $1"
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $1" >&2
}


parse_yaml_with_yq() {
    if ! command -v yq &> /dev/null; then
        log_error "yq is required but not installed. Please install yq first."
        log_error "Install: sudo wget -qO /usr/local/bin/yq https://github.com/mikefarah/yq/releases/latest/download/yq_linux_amd64 && sudo chmod +x /usr/local/bin/yq"
        exit 1
    fi
    

    INPUT_FASTA=$(yq -r '.INPUT_FASTA // ""' "$CONFIG_FILE")
    OUTPUT_DIR=$(yq -r '.OUTPUT_DIR // ""' "$CONFIG_FILE")
    BLAST_OUTFMT=$(yq -r '.BLAST_OUTFMT // "6"' "$CONFIG_FILE")
    BLAST_THREADS=$(yq -r '.BLAST_THREADS // "4"' "$CONFIG_FILE")
    IDENTITY_THRESHOLD=$(yq -r '.IDENTITY_THRESHOLD // "90"' "$CONFIG_FILE")
    GENE_PREFIX=$(yq -r '.GENE_PREFIX // "PF"' "$CONFIG_FILE")
    START_NUMBER=$(yq -r '.START_NUMBER // "1"' "$CONFIG_FILE")
    
    local prefix_lower=$(echo "$GENE_PREFIX" | tr '[:upper:]' '[:lower:]')
    INTERMEDIATE_DIR="$OUTPUT_DIR/intermediate"
    

    BLAST_DB_NAME="${prefix_lower}_db"
    BLAST_OUTPUT="${prefix_lower}_blast.csv"
    RENAMED_FASTA="${prefix_lower}_renamed.fasta"
    
  
    GROUPS_OUTPUT="${prefix_lower}_groups.csv"
    REPRESENTATIVE_FASTA="${GENE_PREFIX}_representatives.fasta"
    REPRESENTATIVE_PAIRS="${prefix_lower}_pairs.tsv"
    FINAL_LIST="${GENE_PREFIX}_families.list"

    local script_dir=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
    local scripts_dir="$script_dir"
    
    CLUSTER_SCRIPT="$scripts_dir/cluster_genes.R"
    RENAME_SCRIPT="$scripts_dir/rename_fa.py"
    LONGEST_SCRIPT="$scripts_dir/longgest_protein_ID.pl"
    
    log_info "Using scripts directory: $scripts_dir"
    log_info "Final result files (in $OUTPUT_DIR):"
    log_info "  - Gene groups: $GROUPS_OUTPUT"
    log_info "  - Representative FASTA: $REPRESENTATIVE_FASTA"
    log_info "  - Gene pairs: $REPRESENTATIVE_PAIRS"
    log_info "  - Final families list: $FINAL_LIST"
    log_info "Intermediate files (in $INTERMEDIATE_DIR):"
    log_info "  - BLAST database: $BLAST_DB_NAME"
    log_info "  - BLAST results: $BLAST_OUTPUT"
    log_info "  - Renamed FASTA: $RENAMED_FASTA"
    
    VERBOSE=$(yq -r '.VERBOSE // "false"' "$CONFIG_FILE")
    CLEANUP_INTERMEDIATE=$(yq -r '.CLEANUP_INTERMEDIATE // "true"' "$CONFIG_FILE")
}


load_config() {
    if [[ ! -f "$CONFIG_FILE" ]]; then
        log_error "Configuration file not found: $CONFIG_FILE"
        exit 1
    fi
    
    log_info "Loading configuration from: $CONFIG_FILE"
    parse_yaml_with_yq
}


check_dependencies() {
    log_info "Checking dependencies..."
    
    local deps=("makeblastdb" "blastp" "Rscript" "python3" "perl" "awk")
    local missing=()
    
    for dep in "${deps[@]}"; do
        if ! command -v "$dep" &> /dev/null; then
            missing+=("$dep")
        fi
    done
    
    if [[ ${#missing[@]} -gt 0 ]]; then
        log_error "Missing dependencies: ${missing[*]}"
        exit 1
    fi
    

    local scripts=("$CLUSTER_SCRIPT" "$RENAME_SCRIPT" "$LONGEST_SCRIPT")
    for script in "${scripts[@]}"; do
        if [[ ! -f "$script" ]]; then
            log_error "Script not found: $script"
            exit 1
        fi
    done
    
    log_info "All dependencies found"
}


check_input_files() {
    log_info "Checking input files..."
    
    if [[ ! -f "$INPUT_FASTA" ]]; then
        log_error "Input FASTA file not found: $INPUT_FASTA"
        exit 1
    fi
    

    if ! head -1 "$INPUT_FASTA" | grep -q "^>"; then
        log_error "Input file doesn't appear to be a valid FASTA file"
        exit 1
    fi
    
    local seq_count=$(grep -c "^>" "$INPUT_FASTA")
    log_info "Found $seq_count sequences in input file"
}


setup_directories() {
    log_info "Setting up directories..."
    

    mkdir -p "$OUTPUT_DIR"
    mkdir -p "$INTERMEDIATE_DIR"
    
    log_info "Created directories:"
    log_info "  - Main output: $OUTPUT_DIR"
    log_info "  - Intermediate files: $INTERMEDIATE_DIR"
    

    LOGFILE="$OUTPUT_DIR/pipeline.log"
    exec 1> >(tee -a "$LOGFILE")
    exec 2> >(tee -a "$LOGFILE" >&2)
}

# Step 1: Create BLAST database
create_blast_database() {
    log_info "Step 1: Creating BLAST database"
    
    local db_path="$INTERMEDIATE_DIR/$BLAST_DB_NAME"
    
    makeblastdb \
        -in "$INPUT_FASTA" \
        -dbtype prot \
        -out "$db_path" \
        -logfile "$OUTPUT_DIR/makeblastdb.log"
    
    log_info "BLAST database created: $db_path"
}

# Step 2: Run BLASTP self-alignment
run_blastp() {
    log_info "Step 2: Running BLASTP self-alignment"
    
    local blast_output="$INTERMEDIATE_DIR/$BLAST_OUTPUT"
    local db_path="$INTERMEDIATE_DIR/$BLAST_DB_NAME"
    
    blastp \
        -query "$INPUT_FASTA" \
        -db "$db_path" \
        -outfmt "$BLAST_OUTFMT" \
        -out "$blast_output" \
        -num_threads "$BLAST_THREADS"
    
    local hit_count=$(wc -l < "$blast_output")
    log_info "BLASTP completed with $hit_count hits"
}

# Step 3: Gene clustering
cluster_genes() {
    log_info "Step 3: Clustering genes into families"
    
    local blast_output="$INTERMEDIATE_DIR/$BLAST_OUTPUT"
    local groups_output="$OUTPUT_DIR/$GROUPS_OUTPUT"
    
    Rscript "$CLUSTER_SCRIPT" \
        --input "$blast_output" \
        --output "$groups_output" \
        --identity_threshold "$IDENTITY_THRESHOLD" \
        --prefix "$GENE_PREFIX" \
        --start_number "$START_NUMBER"
    
    local family_count=$(tail -n +2 "$groups_output" | cut -f2 | sort -u | wc -l)
    log_info "Generated $family_count gene families"
}

# Step 4: Sequence ID renaming
rename_sequences() {
    log_info "Step 4: Renaming sequence IDs"
    
    local groups_output="$OUTPUT_DIR/$GROUPS_OUTPUT"
    local renamed_output="$INTERMEDIATE_DIR/$RENAMED_FASTA"
    
    python3 "$RENAME_SCRIPT" \
        -m "$groups_output" \
        -i "$INPUT_FASTA" \
        -o "$renamed_output"
    
    local renamed_count=$(grep -c "^>" "$renamed_output")
    log_info "Renamed $renamed_count sequences"
}

# Step 5: Extract representative sequences
extract_representatives() {
    log_info "Step 5: Extracting representative sequences"
    
    local renamed_output="$INTERMEDIATE_DIR/$RENAMED_FASTA"
    local representative_output="$OUTPUT_DIR/$REPRESENTATIVE_FASTA"
    
    perl "$LONGEST_SCRIPT" \
        -i "$renamed_output" \
        -o "$representative_output"
    
    local rep_count=$(grep -c "^>" "$representative_output")
    log_info "Extracted $rep_count representative sequences"
}

# Step 6: Generate final gene family list
generate_final_list() {
    log_info "Step 6: Generating final gene family list"
    
    local representative_output="$OUTPUT_DIR/$REPRESENTATIVE_FASTA"
    local representative_pairs="$OUTPUT_DIR/$REPRESENTATIVE_PAIRS"
    local final_list="$OUTPUT_DIR/$FINAL_LIST"
    

    grep "^>" "$representative_output" | sed 's/>//' | \
    awk '{print $1"\t"$1}' > "$representative_pairs"

    shuf "$representative_pairs" | sort -t$'\t' -k1,1 -u > "$final_list"
    
    local list_count=$(wc -l < "$final_list")
    log_info "Generated final list with $list_count gene families"
}

# Clean up intermediate files
cleanup_intermediate_files() {
    if [[ "$CLEANUP_INTERMEDIATE" == "true" ]]; then
        log_info "Cleaning up intermediate files..."
        
        rm -rf "$INTERMEDIATE_DIR"
        rm -f "$OUTPUT_DIR"/*.log
        
        log_info "Cleanup completed"
    fi
}

# Generate analysis report
generate_report() {
    log_info "Generating analysis report"
    
    local report_file="$OUTPUT_DIR/analysis_report.txt"
    
    cat > "$report_file" << EOF
Pan-genome Gene Family Analysis Report
====================================

Analysis Date: $(date)
Input File: $INPUT_FASTA
Output Directory: $OUTPUT_DIR

Parameters:
- Identity Threshold: $IDENTITY_THRESHOLD%
- Gene Prefix: $GENE_PREFIX
- BLAST Threads: $BLAST_THREADS

Results:
- Input Sequences: $(grep -c "^>" "$INPUT_FASTA")
- Gene Families: $(tail -n +2 "$OUTPUT_DIR/$GROUPS_OUTPUT" | cut -f2 | sort -u | wc -l)
- Representative Sequences: $(grep -c "^>" "$OUTPUT_DIR/$REPRESENTATIVE_FASTA")
- Final Gene Family List: $(wc -l < "$OUTPUT_DIR/$FINAL_LIST") entries

Output Files:
1. $GROUPS_OUTPUT - All gene family ID mappings
2. $REPRESENTATIVE_FASTA - Representative protein sequences
3. $RENAMED_FASTA - All renamed protein sequences
4. $FINAL_LIST - Final gene family list for analysis

EOF

    log_info "Analysis report saved to: $report_file"
}

# Main function
main() {
    local start_time=$(date +%s)
    
    log_info "Starting Pan-genome Gene Family Analysis Pipeline"
    log_info "================================================"
    
    # Parse command-line arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -c|--config)
                CONFIG_FILE="$2"
                shift 2
                ;;
            -h|--help)
                usage
                exit 0
                ;;
            -v|--version)
                version
                exit 0
                ;;
            *)
                log_error "Unknown option: $1"
                usage
                exit 1
                ;;
        esac
    done
    
    # Run pipeline
    load_config
    check_dependencies
    check_input_files
    setup_directories
    
    create_blast_database
    run_blastp
    cluster_genes
    rename_sequences
    extract_representatives
    generate_final_list
    
    cleanup_intermediate_files
    generate_report
    
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    log_info "================================================"
    log_info "Pipeline completed successfully!"
    log_info "Total runtime: ${duration} seconds"
    log_info "Results available in: $OUTPUT_DIR"
    log_info "Key output files:"
    log_info "  - $OUTPUT_DIR/$FINAL_LIST (Final gene family list)"
    log_info "  - $OUTPUT_DIR/$GROUPS_OUTPUT (Gene family mappings)"
    log_info "  - $OUTPUT_DIR/$REPRESENTATIVE_FASTA (Representative sequences)"
    log_info "  - $OUTPUT_DIR/$RENAMED_FASTA (All renamed sequences)"
}

# Error handling
trap 'log_error "Pipeline failed at line $LINENO. Exit code: $?"' ERR

# Run main function
main "$@"
