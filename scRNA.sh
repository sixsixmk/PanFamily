#!/bin/bash

# ============================================================================
# PanFamily scRNA Module - Single-cell RNA Analysis
# ============================================================================

set -euo pipefail

# Program information
PROGRAM_NAME="PanFamily scRNA"

# Logging functions
log_info() {
    echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_error() {
    echo "[ERROR] $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2
}

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
SC_R_SCRIPT="${SCRIPT_DIR}/SC.R"

# Print usage
print_usage() {
    cat << EOF
$PROGRAM_NAME - Single-cell RNA Analysis Module

USAGE:
    $0 -c <config.yaml>

REQUIRED ARGUMENTS:
    -c <FILE>        YAML configuration file

OPTIONS:
    -h, --help       Show this help message

DESCRIPTION:
    This module performs single-cell RNA analysis using the SC.R script.
    All analysis parameters should be specified in the YAML configuration file.

EXAMPLE:
    $0 -c Scripts/SC.yaml

EOF
}

# Check dependencies
check_dependencies() {
    log_info "Checking dependencies..."
    
    if ! command -v Rscript &> /dev/null; then
        log_error "Rscript not found. Please install R and ensure it's in your PATH."
        exit 1
    fi
    
    if [[ ! -f "$SC_R_SCRIPT" ]]; then
        log_error "SC.R script not found: $SC_R_SCRIPT"
        exit 1
    fi
    
    log_info "All dependencies are available."
}

# Parse arguments
CONFIG=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -c)
            CONFIG="$2"
            shift 2
            ;;
        -h|--help)
            print_usage
            exit 0
            ;;
        *)
            log_error "Unknown argument: $1"
            print_usage
            exit 1
            ;;
    esac
done

# Validate arguments
if [[ -z "$CONFIG" ]]; then
    log_error "Missing required argument: -c <config.yaml>"
    print_usage
    exit 1
fi

if [[ ! -f "$CONFIG" ]]; then
    log_error "Configuration file not found: $CONFIG"
    exit 1
fi

# Check dependencies
check_dependencies

# Run SC.R script
log_info "Starting scRNA analysis..."
log_info "Configuration file: $CONFIG"
log_info "Calling SC.R script: $SC_R_SCRIPT"

# Execute the R script with the configuration file
Rscript "$SC_R_SCRIPT" -c "$CONFIG"

EXIT_CODE=$?

if [[ $EXIT_CODE -eq 0 ]]; then
    log_info "scRNA analysis completed successfully!"
else
    log_error "scRNA analysis failed with exit code: $EXIT_CODE"
    exit $EXIT_CODE
fi

