#!/bin/bash

# ============================================================================
# MEME.sh - Single-genome Motif Analysis Module
# Motif discovery for protein/nucleotide sequences using MEME
# ============================================================================

set -euo pipefail

# Logging functions
log_info() {
    echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_error() {
    echo "[ERROR] $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2
}

log_warn() {
    echo "[WARN] $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

# Show help message
show_help() {
    cat << EOF
MEME.sh - Single-genome Motif Analysis Module

Usage:
    bash MEME.sh -c <config.yaml>
    bash MEME.sh --input <fasta> --output <dir> [options]

Configuration file parameters (YAML):
    input_fasta      - Input FASTA sequence file (required)
    output_dir       - Output directory (required)
    seq_type         - Sequence type: protein or dna (default: protein)
    nmotifs          - Number of motifs to discover (default: 10)
    minw             - Minimum motif width (default: 6)
    maxw             - Maximum motif width (default: 50)
    meme_options     - Additional MEME parameters (optional)
    threads          - Number of threads (default: 4)

Command-line options:
    -c, --config     Path to configuration file
    -i, --input      Input FASTA file
    -o, --output     Output directory
    -t, --type       Sequence type (protein/dna)
    -n, --nmotifs    Number of motifs
    --minw           Minimum motif width
    --maxw           Maximum motif width
    --threads        Number of threads
    -h, --help       Show help message

Examples:
    # Using a configuration file
    bash MEME.sh -c motif.yaml

    # Using command-line arguments (protein sequences)
    bash MEME.sh -i proteins.fa -o meme_output -t protein -n 10

    # Using command-line arguments (CDS sequences)
    bash MEME.sh -i cds.fa -o meme_output -t dna -n 5

EOF
    exit 0
}

# Parse YAML configuration file
parse_yaml() {
    local yaml_file="$1"
    
    if [[ ! -f "$yaml_file" ]]; then
        log_error "Configuration file not found: $yaml_file"
        exit 1
    fi
    
    # Parse YAML parameters
    INPUT_FASTA=$(grep -E "^input_fasta:" "$yaml_file" | sed 's/input_fasta:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "")
    OUTPUT_DIR=$(grep -E "^output_dir:" "$yaml_file" | sed 's/output_dir:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "")
    SEQ_TYPE=$(grep -E "^seq_type:" "$yaml_file" | sed 's/seq_type:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "protein")
    NMOTIFS=$(grep -E "^nmotifs:" "$yaml_file" | sed 's/nmotifs:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "10")
    MINW=$(grep -E "^minw:" "$yaml_file" | sed 's/minw:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "6")
    MAXW=$(grep -E "^maxw:" "$yaml_file" | sed 's/maxw:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "50")
    MEME_OPTIONS=$(grep -E "^meme_options:" "$yaml_file" | sed 's/meme_options:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "")
    THREADS=$(grep -E "^threads:" "$yaml_file" | sed 's/threads:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "4")
}

# Check dependencies
check_dependencies() {
    log_info "Checking dependencies..."
    
    if ! command -v meme &> /dev/null; then
        log_error "MEME software not found. Please install MEME Suite first."
        log_error "Installation method: conda install -c bioconda meme"
        exit 1
    fi
    
    log_info "MEME version: $(meme -version 2>&1 | head -1 || echo 'unknown')"
}

# Validate inputs
validate_inputs() {
    if [[ -z "$INPUT_FASTA" ]]; then
        log_error "Input FASTA file not specified"
        exit 1
    fi
    
    if [[ ! -f "$INPUT_FASTA" ]]; then
        log_error "Input file does not exist: $INPUT_FASTA"
        exit 1
    fi
    
    if [[ -z "$OUTPUT_DIR" ]]; then
        log_error "Output directory not specified"
        exit 1
    fi
    
    # Check sequence type
    if [[ "$SEQ_TYPE" != "protein" && "$SEQ_TYPE" != "dna" ]]; then
        log_warn "Unknown sequence type '$SEQ_TYPE', using default value 'protein'"
        SEQ_TYPE="protein"
    fi
    
    # Check number of sequences
    local seq_count=$(grep -c "^>" "$INPUT_FASTA" || echo "0")
    if [[ "$seq_count" -lt 2 ]]; then
        log_error "MEME requires at least 2 sequences, only $seq_count provided"
        exit 1
    fi
    log_info "Number of input sequences: $seq_count"
}

# Run MEME analysis
run_meme() {
    log_info "=========================================="
    log_info "Starting MEME Motif Analysis"
    log_info "=========================================="
    log_info "Input file: $INPUT_FASTA"
    log_info "Output directory: $OUTPUT_DIR"
    log_info "Sequence type: $SEQ_TYPE"
    log_info "Number of motifs: $NMOTIFS"
    log_info "Motif width: $MINW - $MAXW"
    log_info "Number of threads: $THREADS"
    
    # Create output directory
    mkdir -p "$OUTPUT_DIR"
    
    # Build MEME command
    local meme_cmd="meme $INPUT_FASTA -oc $OUTPUT_DIR -nmotifs $NMOTIFS -minw $MINW -maxw $MAXW"
    
    # Add sequence type parameter
    if [[ "$SEQ_TYPE" == "protein" ]]; then
        meme_cmd="$meme_cmd -protein"
    else
        meme_cmd="$meme_cmd -dna"
    fi
    
    # Add extra options
    if [[ -n "$MEME_OPTIONS" ]]; then
        meme_cmd="$meme_cmd $MEME_OPTIONS"
    fi
    
    log_info "Executing command: $meme_cmd"
    
    # Run MEME
    eval "$meme_cmd"
    
    local exit_code=$?
    if [[ $exit_code -ne 0 ]]; then
        log_error "MEME analysis failed with exit code: $exit_code"
        exit 1
    fi
    
    log_info "MEME analysis completed"
}

# Generate analysis report
generate_report() {
    log_info "Generating analysis report..."
    
    local report_file="${OUTPUT_DIR}/motif_analysis_report.txt"
    
    cat > "$report_file" << EOF
================================================================================
MEME Motif Analysis Report
================================================================================
Generated: $(date '+%Y-%m-%d %H:%M:%S')

Input Parameters:
-----------------
Input FASTA: $INPUT_FASTA
Sequence Type: $SEQ_TYPE
Number of Motifs: $NMOTIFS
Motif Width Range: $MINW - $MAXW

Output Files:
-------------
$(ls -la "$OUTPUT_DIR" 2>/dev/null || echo "No files found")

Key Output Files:
- meme.html    : Interactive HTML report (open in browser)
- meme.txt     : Text format results
- meme.xml     : XML format results (for downstream analysis)
- logo*.png    : Sequence logos for each motif

================================================================================
EOF
    
    log_info "Report saved to: $report_file"
}

# Main function
main() {
    # Default values
    INPUT_FASTA=""
    OUTPUT_DIR=""
    SEQ_TYPE="protein"
    NMOTIFS="10"
    MINW="6"
    MAXW="50"
    MEME_OPTIONS=""
    THREADS="4"
    CONFIG_FILE=""
    
    # Parse command-line arguments
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -h|--help)
                show_help
                ;;
            -c|--config)
                CONFIG_FILE="$2"
                shift 2
                ;;
            -i|--input)
                INPUT_FASTA="$2"
                shift 2
                ;;
            -o|--output)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            -t|--type)
                SEQ_TYPE="$2"
                shift 2
                ;;
            -n|--nmotifs)
                NMOTIFS="$2"
                shift 2
                ;;
            --minw)
                MINW="$2"
                shift 2
                ;;
            --maxw)
                MAXW="$2"
                shift 2
                ;;
            --threads)
                THREADS="$2"
                shift 2
                ;;
            *)
                log_error "Unknown parameter: $1"
                show_help
                ;;
        esac
    done
    
    # Parse configuration file if provided
    if [[ -n "$CONFIG_FILE" ]]; then
        parse_yaml "$CONFIG_FILE"
    fi
    
    # Check dependencies
    check_dependencies
    
    # Validate inputs
    validate_inputs
    
    # Run MEME
    run_meme
    
    # Generate report
    generate_report
    
    log_info "=========================================="
    log_info "Motif analysis completed!"
    log_info "Result directory: $OUTPUT_DIR"
    log_info "Please open meme.html in a browser to view the results"
    log_info "=========================================="
}

# Error handling
trap 'log_error "Script execution failed at line: $LINENO"' ERR

main "$@"
