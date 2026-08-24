#!/bin/bash

set -euo pipefail

log_info() {
    echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_error() {
    echo "[ERROR] $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2
}

log_warn() {
    echo "[WARN] $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"

show_help() {
    cat << EOF
Collinearity.sh - Collinearity Analysis Module (MCScanX)

Usage:
    bash Collinearity.sh -c <config.yaml>

Configuration file parameters (YAML):
    # Basic parameters
    project_name     - Project name/prefix (required, used for output naming)
    output_dir       - Output directory (required)
    
    # Input option 1: provide raw data, auto preparation
    gff_file         - GFF3 gene annotation file
    protein_file     - Protein FASTA file
    
    # Input option 2: provide prepared MCScanX format files
    mcscanx_gff      - MCScanX formatted .gff file
    mcscanx_blast    - MCScanX formatted .blast file
    
    # MCScanX parameters
    match_score      - Match score threshold (default: 50)
    match_size       - Minimum number of genes per collinearity block (default: 5)
    gap_penalty      - Gap penalty (default: -1)
    overlap_window   - Overlap window size (default: 5)
    e_value          - BLAST E-value threshold (default: 1e-10)
    num_threads      - Number of BLAST threads (default: 4)
    
    # Visualization parameters
    run_visualization - Whether to run visualization (default: true)
    plot_types       - Plot types: dot, dual, circle, bar (default: all)

Command-line options:
    -c, --config     Configuration file path (required)
    -h, --help       Show help information

Input file formats:
    1. GFF3 format (standard gene annotation):
       Chr1  .  gene  1000  2000  .  +  .  ID=Gene1;...
    
    2. MCScanX gff format:
       sp#Gene1  Chr1  1000  2000
       sp#Gene2  Chr1  3000  4000
    
    3. BLAST output format (-outfmt 6):
       Gene1  Gene2  95.00  200  10  0  1  200  1  200  1e-50  400

Output files:
    <project>.gff          - MCScanX formatted gene position file
    <project>.blast        - BLASTP alignment results
    <project>.collinearity - Collinearity block results
    <project>.html         - Collinearity statistics (if available)
    visualization/         - Visualization output directory
        dot_plot.png
        dual_synteny.png
        circle_plot.png
        bar_plot.png

Example:
    bash Collinearity.sh -c collinearity.yaml

EOF
    exit 0
}

# Parse YAML configuration file
parse_yaml() {
    local yaml_file="$1"
    
    if [[ ! -f "$yaml_file" ]]; then
        log_error "Configuration file does not exist: $yaml_file"
        exit 1
    fi
    
    log_info "Parsing configuration file: $yaml_file"
    
    # Parse basic parameters
    PROJECT_NAME=$(grep -E "^project_name:" "$yaml_file" | sed 's/project_name:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "")
    OUTPUT_DIR=$(grep -E "^output_dir:" "$yaml_file" | sed 's/output_dir:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "")
    
    # Input option 1: raw data
    GFF_FILE=$(grep -E "^gff_file:" "$yaml_file" | sed 's/gff_file:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "")
    PROTEIN_FILE=$(grep -E "^protein_file:" "$yaml_file" | sed 's/protein_file:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "")
    
    # Input option 2: prepared files
    MCSCANX_GFF=$(grep -E "^mcscanx_gff:" "$yaml_file" | sed 's/mcscanx_gff:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "")
    MCSCANX_BLAST=$(grep -E "^mcscanx_blast:" "$yaml_file" | sed 's/mcscanx_blast:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "")
    
    # MCScanX parameters
    MATCH_SCORE=$(grep -E "^match_score:" "$yaml_file" | sed 's/match_score:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "50")
    MATCH_SIZE=$(grep -E "^match_size:" "$yaml_file" | sed 's/match_size:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "5")
    GAP_PENALTY=$(grep -E "^gap_penalty:" "$yaml_file" | sed 's/gap_penalty:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "-1")
    OVERLAP_WINDOW=$(grep -E "^overlap_window:" "$yaml_file" | sed 's/overlap_window:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "5")
    E_VALUE=$(grep -E "^e_value:" "$yaml_file" | sed 's/e_value:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "1e-10")
    NUM_THREADS=$(grep -E "^num_threads:" "$yaml_file" | sed 's/num_threads:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "4")
    
    # Visualization parameters
    RUN_VISUALIZATION=$(grep -E "^run_visualization:" "$yaml_file" | sed 's/run_visualization:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "true")
    PLOT_TYPES=$(grep -E "^plot_types:" "$yaml_file" | sed 's/plot_types:[[:space:]]*//' | tr -d '"' | tr -d "'" | tr -d '[]' || echo "dot,dual,circle,bar")
}

# Check dependencies
check_dependencies() {
    log_info "Checking software dependencies..."
    
    # Check MCScanX
    if ! command -v MCScanX &> /dev/null; then
        log_error "MCScanX not found, please install it first"
        log_error "Installation method:"
        log_error "  git clone https://github.com/wyp1125/MCScanX.git"
        log_error "  cd MCScanX && make"
        log_error "  export PATH=\$PATH:/path/to/MCScanX"
        exit 1
    fi
    
    # Check BLAST+
    if ! command -v blastp &> /dev/null; then
        log_error "BLAST+ not found, please install: conda install -c bioconda blast"
        exit 1
    fi
    
    if ! command -v makeblastdb &> /dev/null; then
        log_error "makeblastdb not found, p_
