#!/bin/bash

# ============================================================================
# PanFamily Gene Family Identification Module
# ============================================================================

set -euo pipefail

# Program information
PROGRAM_NAME="PanFamily GeneDetection"

# Logging functions
log_info() {
    echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_warn() {
    echo "[WARN] $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_error() {
    echo "[ERROR] $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2
}

log_step() {
    echo "[STEP] $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

print_banner() {
    echo "========================================================"
    echo "            Gene Family Identification Module"
    echo "========================================================"
}

# Print usage
print_usage() {
    cat << EOF
$PROGRAM_NAME - Gene Family Identification Module

USAGE:
    $0 -c <FILE> [OPTIONS]
    OR
    $0 --method <blastp|hmmer|both> -i <INPUT_FASTA> -o <OUTPUT_DIR> [--seed <SEED_FASTA>] [--hmm <HMM1,HMM2,...>] [--threads N] [--evalue-blast EV] [--evalue-hmm EV] [--coverage-threshold CV] [OPTIONS]

Required arguments:
    -c <FILE>               YAML configuration file (e.g., config.yaml)
    OR parameter mode:
    --method                Select method: blastp | hmmer | both
    -i, --input             Input protein FASTA file
    -o, --output            Output directory

Optional arguments:
    --keep-intermediate     Keep intermediate files
    --seed                  Seed protein file for BLASTP mode
    --hmm                   HMM file(s) for HMMER mode, comma-separated
    --threads               Number of threads (default: 12)
    --evalue-blast          E-value threshold for BLASTP (default: 1e-5)
    --evalue-hmm            E-value threshold for HMMER (default: 1e-5)
    --coverage-threshold    Coverage threshold for HMMER domain filtering (default: 0.9)
    -h, --help              Show this help message

EXAMPLES:
    $0 -c config.yaml
    $0 --method both -i proteins.fa -o out_dir --seed seeds.fa --hmm A.hmm,B.hmm --threads 16 --coverage-threshold 0.8
EOF
}

check_dependencies() {
    local required_tools=("yq" "seqkit")
    local missing_tools=()

    if ! command -v yq &> /dev/null; then
        missing_tools+=("yq")
    fi
    
    if [[ "$RUN_BLASTP" = true ]]; then
        required_tools+=("makeblastdb" "blastp")
    fi
    if [[ "$RUN_HMMER" = true ]]; then
        required_tools+=("hmmsearch")
    fi

    for tool in "${required_tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            missing_tools+=("$tool")
        fi
    done

    if [ ${#missing_tools[@]} -ne 0 ]; then
        log_error "Missing required tools: ${missing_tools[*]}"
        log_error "Please install these tools before running the analysis"
        exit 1
    fi

    log_info "Dependency check passed"
}

normalize_bool() {
    local v="$1"
    case "$v" in
        true|TRUE|True|t|T|1) echo true ;;
        false|FALSE|False|f|F|0|""|null) echo false ;;
        *) echo "$v" ;;
    esac
}

read_config() {
    log_info "Reading parameters from config file '$CONFIG_FILE'"
    
    if [[ ! -f "$CONFIG_FILE" ]]; then
        log_error "Config file does not exist: $CONFIG_FILE"
        exit 1
    fi
    
    RUN_BLASTP=$(normalize_bool "$(yq -r '.methods.blastp' "$CONFIG_FILE")")
    RUN_HMMER=$(normalize_bool "$(yq -r '.methods.hmmer' "$CONFIG_FILE")")
    
    INPUT_FILE=$(yq -r '.files.input_file' "$CONFIG_FILE")
    OUTPUT_DIR=$(yq -r '.files.output_dir' "$CONFIG_FILE")
    SEED_FILE=$(yq -r '.files.seed' "$CONFIG_FILE")
    HMM_FILES_STR=$(yq -r '.files.hmm' "$CONFIG_FILE")
    IFS=',' read -r -a HMM_FILES <<< "$HMM_FILES_STR"
    
    THREADS=$(yq -r '.parameters.threads' "$CONFIG_FILE")
    EVALUE_BLAST=$(yq -r '.parameters.evalue_blast' "$CONFIG_FILE")
    EVALUE_HMM=$(yq -r '.parameters.evalue_hmm' "$CONFIG_FILE")
    COVERAGE_THRESHOLD=$(yq -r '.parameters.coverage_threshold' "$CONFIG_FILE")
    
    if [[ ! -f "$INPUT_FILE" ]]; then
        log_error "Input file specified in config does not exist: $INPUT_FILE"
        exit 1
    fi
    
    if [[ "$RUN_BLASTP" = true ]] && [[ ! -f "$SEED_FILE" ]]; then
        log_error "Seed file specified in config does not exist: $SEED_FILE"
        exit 1
    fi
    if [[ "$RUN_HMMER" = true ]]; then
        for hmm_file in "${HMM_FILES[@]}"; do
            if [[ ! -f "$hmm_file" ]]; then
                log_error "HMM file specified in config does not exist: $hmm_file"
                exit 1
            fi
        done
    fi

    if [[ "$RUN_BLASTP" != true ]] && [[ "$RUN_HMMER" != true ]]; then
        log_error "No method specified in config (either blastp or hmmer must be true)"
        exit 1
    fi
}

read_args() {
    if [[ -z "${METHOD:-}" ]]; then
        log_error "Missing required parameter: --method"
        echo ""
        print_usage
        exit 1
    fi

    case "$METHOD" in
        blastp)
            RUN_BLASTP=true; RUN_HMMER=false ;;
        hmmer)
            RUN_BLASTP=false; RUN_HMMER=true ;;
        both)
            RUN_BLASTP=true; RUN_HMMER=true ;;
        *)
            log_error "--method only supports blastp|hmmer|both"
            exit 1 ;;
    esac

    if [[ -z "${INPUT_FILE:-}" || -z "${OUTPUT_DIR:-}" ]]; then
        log_error "Missing required parameters: -i/--input and -o/--output"
        echo ""
        print_usage
        exit 1
    fi
    THREADS=${THREADS:-12}
    EVALUE_BLAST=${EVALUE_BLAST:-1e-5}
    EVALUE_HMM=${EVALUE_HMM:-1e-5}
    COVERAGE_THRESHOLD=${COVERAGE_THRESHOLD:-0.9}

    if [[ ! -f "$INPUT_FILE" ]]; then
        log_error "Input file does not exist: $INPUT_FILE"
        exit 1
    fi
    if [[ "$RUN_BLASTP" = true ]]; then
        if [[ -z "${SEED_FILE:-}" || ! -f "$SEED_FILE" ]]; then
            log_error "BLASTP mode requires a --seed file"
            exit 1
        fi
    fi
    if [[ "$RUN_HMMER" = true ]]; then
        if [[ -z "${HMM_FILES_STR:-}" ]]; then
            log_error "HMMER mode requires --hmm file(s) (comma-separated)"
            exit 1
        fi
        IFS=',' read -r -a HMM_FILES <<< "$HMM_FILES_STR"
        for hmm_file in "${HMM_FILES[@]}"; do
            if [[ ! -f "$hmm_file" ]]; then
                log_error "HMM file does not exist: $hmm_file"
                exit 1
            fi
        done
    fi
}

create_directories() {
    log_info "Creating output directory structure"
    mkdir -p "$OUTPUT_DIR"/logs
    
    if [[ "$RUN_BLASTP" = true ]]; then
        mkdir -p "$OUTPUT_DIR"/blastp
    fi
    if [[ "$RUN_HMMER" = true ]]; then
        mkdir -p "$OUTPUT_DIR"/hmmer
    fi
    
    LOG_FILE="$OUTPUT_DIR/logs/identification_$(date +%Y%m%d_%H%M%S).log"
    log_info "Output directory: $OUTPUT_DIR"
    log_info "Log file: $LOG_FILE"
}

run_blastp() {
    log_step "Running BLASTP"
    local db_path="$OUTPUT_DIR/blastp/target_proteins_db"
    local blastp_output="$OUTPUT_DIR/blastp/blastp_result.txt"

    log_info "Building BLAST database from target proteins"
    makeblastdb -in "$INPUT_FILE" -dbtype prot -out "$db_path" 2>> "$LOG_FILE"

    log_info "Querying target protein database with seed proteins (E-value < $EVALUE_BLAST)"
    blastp -query "$SEED_FILE" \
           -db "$db_path" \
           -out "$blastp_output" \
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
           -evalue "$EVALUE_BLAST" \
           -num_threads "$THREADS" 2>> "$LOG_FILE"

    log_info "BLASTP search finished ✓"
}

run_hmmsearch() {
    log_step "Running HMMER search"
    local combined_domtbl="$OUTPUT_DIR/hmmer/hmmsearch_result.domtbl"
    : > "$combined_domtbl"

    log_info "Running hmmsearch for each HMM file and merging results (--noali --cut_tc)"
    for hmm_file in "${HMM_FILES[@]}"; do
        local base_name
        base_name=$(basename "$hmm_file")
        local per_out="$OUTPUT_DIR/hmmer/${base_name}.domtbl"

        log_info "hmmsearch: $base_name"
        hmmsearch --domtblout "$per_out" \
                  --cpu "$THREADS" \
                  --noali \
                  --cut_tc \
                  "$hmm_file" \
                  "$INPUT_FILE" >> "$LOG_FILE" 2>&1

        grep -v '^#' "$per_out" >> "$combined_domtbl"
    done

    log_info "HMMER search finished ✓"
}

integrate_results() {
    log_step "Integrating and filtering results"
    
    local blastp_candidates="$OUTPUT_DIR/blastp/blastp_result.txt"
    local hmmer_domtbl_combined="$OUTPUT_DIR/hmmer/hmmsearch_result.domtbl"
    local hmmer_filtered_ids="$OUTPUT_DIR/hmmer/hmmer_filtered_ids.txt"
    local candidate_genes="$OUTPUT_DIR/family_genes.txt"

    : > "$candidate_genes"

    if [[ "$RUN_HMMER" = true ]] && [[ -s "$hmmer_domtbl_combined" ]]; then
        log_info "Filtering HMMER results: e-value < $EVALUE_HMM and coverage >= $COVERAGE_THRESHOLD"
        awk -v e="$EVALUE_HMM" -v c="$COVERAGE_THRESHOLD" '($7<e && ($19-$18)/$6>=c){print $1}' "$hmmer_domtbl_combined" | sort -u > "$hmmer_filtered_ids"
    fi

    if [[ "$RUN_BLASTP" = true ]] && [[ "$RUN_HMMER" = true ]]; then
        log_info "Integrating BLASTP and HMMER results (intersection)"
        awk '{print $2}' "$blastp_candidates" | sort -u > "$OUTPUT_DIR/blastp_ids.txt"
        
        if [[ -s "$hmmer_filtered_ids" ]]; then
            sort -u "$hmmer_filtered_ids" > "$OUTPUT_DIR/hmmer_ids.txt"
        else
            : > "$OUTPUT_DIR/hmmer_ids.txt"
        fi
        
        if [[ -s "$OUTPUT_DIR/blastp_ids.txt" && -s "$OUTPUT_DIR/hmmer_ids.txt" ]]; then
            comm -12 "$OUTPUT_DIR/blastp_ids.txt" "$OUTPUT_DIR/hmmer_ids.txt" > "$candidate_genes"
            log_info "Intersection of BLASTP and HMMER results complete"
        else
            log_warn "BLASTP or HMMER results empty, unable to compute intersection"
            : > "$candidate_genes"
        fi
        
        rm "$OUTPUT_DIR/blastp_ids.txt" "$OUTPUT_DIR/hmmer_ids.txt"
    elif [[ "$RUN_BLASTP" = true ]]; then
        log_info "Extracting target gene IDs from BLASTP results (deduplicated)"
        awk '{print $2}' "$blastp_candidates" | sort -u > "$candidate_genes"
    elif [[ "$RUN_HMMER" = true ]]; then
        log_info "Extracting gene IDs from HMMER results (deduplicated)"
        if [[ -s "$hmmer_filtered_ids" ]]; then
            sort -u "$hmmer_filtered_ids" > "$candidate_genes"
        else
            : > "$candidate_genes"
        fi
    else
        log_error "No method specified in config, unable to integrate results"
        exit 1
    fi
    
    if [[ ! -s "$candidate_genes" ]]; then
        log_warn "No candidate genes found. Please check your input files and parameters."
        exit 0
    fi
    
    log_info "Number of candidate genes: $(wc -l < "$candidate_genes")"

    log_info "Extracting final protein sequences for candidate genes"
    seqkit grep -f "$candidate_genes" "$INPUT_FILE" > "$OUTPUT_DIR/family_protein.fa"

}

main() {
    print_banner
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            -c)
                CONFIG_FILE="$2"
                shift 2
                ;;
            --method)
                METHOD="$2"
                shift 2
                ;;
            -i|--input)
                INPUT_FILE="$2"
                shift 2
                ;;
            -o|--output)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            --seed)
                SEED_FILE="$2"
                shift 2
                ;;
            --hmm)
                HMM_FILES_STR="$2"
                shift 2
                ;;
            --threads)
                THREADS="$2"
                shift 2
                ;;
            --evalue-blast)
                EVALUE_BLAST="$2"
                shift 2
                ;;
            --evalue-hmm)
                EVALUE_HMM="$2"
                shift 2
                ;;
            --coverage-threshold)
                COVERAGE_THRESHOLD="$2"
                shift 2
                ;;
            --keep-intermediate)
                KEEP_INTERMEDIATE=true
                shift
                ;;
            -h|--help)
                print_usage
                exit 0
                ;;
            *)
                log_error "Unknown argument: $1"
                echo ""
                print_usage
                exit 1
                ;;
        esac
    done

    log_info "Starting $PROGRAM_NAME analysis"

    if [[ -n "${CONFIG_FILE:-}" ]]; then
        read_config
    else
        read_args
    fi
    log_info "Input file: $INPUT_FILE"
    check_dependencies
    create_directories
    
    if [[ "$RUN_BLASTP" = true ]]; then
        run_blastp
    fi
    if [[ "$RUN_HMMER" = true ]]; then
        run_hmmsearch
    fi
    
    integrate_results

    log_info "Results saved in: $OUTPUT_DIR"
    echo ""
}

KEEP_INTERMEDIATE=false
CONFIG_FILE=""
INPUT_FILE=""
OUTPUT_DIR=""
METHOD=""
SEED_FILE=""
HMM_FILES_STR=""
THREADS=""
EVALUE_BLAST=""
EVALUE_HMM=""
COVERAGE_THRESHOLD=""

main "$@"
