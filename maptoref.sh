#!/bin/bash
#######################################################################
# Script: blastp_mapping.sh
# Description: BLASTP-based gene sequence alignment and mapping table generation
#######################################################################

set -e
set -u
set -o pipefail
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' 
log_info() {
    echo -e "${GREEN}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}
usage() {
    cat << EOF
Usage: $0 -c <config_file>

Required arguments:
    -c    Configuration file path

Configuration file example (config.txt):
    REF_PROTEIN=/path/to/reference_proteins.fasta
    QUERY_PROTEIN=/path/to/gene_family_proteins.fasta
    OUTPUT_DIR=/path/to/output
    IDENTITY=100
    EVALUE=1e-5
    THREADS=8

Description:
    - REF_PROTEIN: Reference genome protein sequence file (FASTA)
    - QUERY_PROTEIN: Gene family protein sequence file (FASTA)
    - OUTPUT_DIR: Output directory
    - IDENTITY: Identity threshold (default 100, range 0–100)
    - EVALUE: E-value threshold (default 1e-5)
    - THREADS: Number of threads (default 8)

Output files:
    - blastp_result.txt: Raw BLASTP results
    - filtered_result.txt: Filtered results
    - gene_mapping.txt: Gene mapping table (query_id -> subject_id)

EOF
    exit 1
}

check_dependencies() {
    local tools=("makeblastdb" "blastp")
    for tool in "${tools[@]}"; do
        if ! command -v $tool &> /dev/null; then
            log_error "Required tool not found: $tool"
            log_error "Please ensure BLAST+ is installed and added to PATH"
            exit 1
        fi
    done
    log_info "Dependency check passed"
}

parse_config() {
    local config_file=$1
    
    if [ ! -f "$config_file" ]; then
        log_error "Configuration file not found: $config_file"
        exit 1
    fi
    
    log_info "Parsing configuration file: $config_file"
    while IFS='=' read -r key value; do
        [[ $key =~ ^#.*$ ]] && continue
        [[ -z $key ]] && continue
        key=$(echo "$key" | xargs)
        value=$(echo "$value" | xargs)
        
        case $key in
            REF_PROTEIN) REF_PROTEIN=$value ;;
            QUERY_PROTEIN) QUERY_PROTEIN=$value ;;
            OUTPUT_DIR) OUTPUT_DIR=$value ;;
            IDENTITY) IDENTITY=$value ;;
            EVALUE) EVALUE=$value ;;
            THREADS) THREADS=$value ;;
        esac
    done < "$config_file"
    IDENTITY=${IDENTITY:-100}
    EVALUE=${EVALUE:-1e-5}
    THREADS=${THREADS:-8}
}

validate_config() {
    log_info "Validating configuration parameters..."
    
    local errors=0
    
    # Check required files
    if [ -z "${REF_PROTEIN:-}" ] || [ ! -f "$REF_PROTEIN" ]; then
        log_error "Reference protein file not found: ${REF_PROTEIN:-not set}"
        ((errors++))
    fi
    
    if [ -z "${QUERY_PROTEIN:-}" ] || [ ! -f "$QUERY_PROTEIN" ]; then
        log_error "Query protein file not found: ${QUERY_PROTEIN:-not set}"
        ((errors++))
    fi
    
    if [ -z "${OUTPUT_DIR:-}" ]; then
        log_error "Output directory not set"
        ((errors++))
    fi
    
    # Check IDENTITY range
    if [ "$IDENTITY" -lt 0 ] || [ "$IDENTITY" -gt 100 ]; then
        log_error "IDENTITY must be between 0 and 100: $IDENTITY"
        ((errors++))
    fi
    
    if [ $errors -gt 0 ]; then
        log_error "Configuration validation failed with $errors error(s)"
        exit 1
    fi
    
    log_info "Configuration validation passed"
}
count_sequences() {
    local fasta_file=$1
    grep -c "^>" "$fasta_file" || echo "0"
}

# Build BLAST database
build_blast_db() {
    log_info "========== Step 1: Build BLAST database =========="
    
    local db_dir="$OUTPUT_DIR/blast_db"
    mkdir -p "$db_dir"
    
    local db_name="$db_dir/ref_proteins"
    
    log_info "Reference sequence file: $REF_PROTEIN"
    log_info "Database name: $db_name"
    
    # Count sequences
    local seq_count=$(count_sequences "$REF_PROTEIN")
    log_info "Number of reference sequences: $seq_count"
    
    # Build database
    makeblastdb \
        -in "$REF_PROTEIN" \
        -dbtype prot \
        -out "$db_name" \
        -parse_seqids \
        > "$OUTPUT_DIR/makeblastdb.log" 2>&1
    
    if [ $? -eq 0 ]; then
        log_info "BLAST database built successfully"
    else
        log_error "BLAST database build failed, see log: $OUTPUT_DIR/makeblastdb.log"
        exit 1
    fi
    
    echo "$db_name" > "$OUTPUT_DIR/db_path.txt"
}
run_blastp() {
    log_info "========== Step 2: Run BLASTP alignment =========="
    
    local db_name=$(cat "$OUTPUT_DIR/db_path.txt")
    local blast_output="$OUTPUT_DIR/blastp_result.txt"
    
    log_info "Query sequences: $QUERY_PROTEIN"
    log_info "Database: $db_name"
    log_info "E-value threshold: $EVALUE"
    log_info "Threads: $THREADS"
    local query_count=$(count_sequences "$QUERY_PROTEIN")
    log_info "Number of query sequences: $query_count"
    
    log_info "Starting alignment..."
    
    blastp \
        -query "$QUERY_PROTEIN" \
        -db "$db_name" \
        -out "$blast_output" \
        -evalue "$EVALUE" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
        -num_threads "$THREADS" \
        -max_target_seqs 5 \
        > "$OUTPUT_DIR/blastp.log" 2>&1
    
    if [ $? -eq 0 ]; then
        log_info "BLASTP alignment completed"
        local hit_count=$(wc -l < "$blast_output" || echo "0")
        log_info "Number of hits: $hit_count"
    else
        log_error "BLASTP failed, see log: $OUTPUT_DIR/blastp.log"
        exit 1
    fi
    
    local blast_with_header="$OUTPUT_DIR/blastp_result_with_header.txt"
    echo -e "query_id\tsubject_id\tidentity\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen" > "$blast_with_header"
    cat "$blast_output" >> "$blast_with_header"
    
    log_info "Result with header saved: $blast_with_header"
}

filter_results() {
    log_info "========== Step 3: Filter alignment results =========="
    
    local blast_output="$OUTPUT_DIR/blastp_result.txt"
    local filtered_output="$OUTPUT_DIR/filtered_result.txt"
    
    log_info "Identity threshold: ${IDENTITY}%"
    
    if [ ! -f "$blast_output" ] || [ ! -s "$blast_output" ]; then
        log_error "BLASTP result file is empty or missing"
        exit 1
    fi
    
    awk -v threshold="$IDENTITY" '$3 >= threshold' "$blast_output" > "$filtered_output"
    
    local filtered_count=$(wc -l < "$filtered_output" || echo "0")
    log_info "Filtered result count: $filtered_count"
    
    if [ "$filtered_count" -eq 0 ]; then
        log_warn "No alignments with identity >= ${IDENTITY}% found"
        touch "$OUTPUT_DIR/gene_mapping.txt"
        return
    fi
    
    local filtered_with_header="$OUTPUT_DIR/filtered_result_with_header.txt"
    echo -e "query_id\tsubject_id\tidentity\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen" > "$filtered_with_header"
    cat "$filtered_output" >> "$filtered_with_header"
    
    log_info "Filtered results saved: $filtered_output"
}
generate_mapping() {
    log_info "========== Step 4: Generate gene mapping table =========="
    
    local filtered_output="$OUTPUT_DIR/filtered_result.txt"
    local mapping_output="$OUTPUT_DIR/gene_mapping.txt"
    
    if [ ! -f "$filtered_output" ] || [ ! -s "$filtered_output" ]; then
        log_warn "Filtered results are empty, mapping table will not be generated"
        touch "$mapping_output"
        return
    fi
    
    awk '{print $1"\t"$2"\t"$12}' "$filtered_output" | \
    sort -k1,1 -k3,3rn | \
    awk '!seen[$1]++' | \
    cut -f1,2 > "$mapping_output"
    
    local mapping_count=$(wc -l < "$mapping_output" || echo "0")
    log_info "Number of mappings: $mapping_count"
    log_info "Mapping table saved: $mapping_output"
    
    if [ "$mapping_count" -gt 0 ]; then
        log_info "Mapping examples (first 10 lines):"
        head -n 10 "$mapping_output" | while read -r line; do
            echo "  $line"
        done
    fi
}

generate_report() {
    log_info "========== Step 5: Generate analysis report =========="
    
    local report_file="$OUTPUT_DIR/blastp_report.txt"
    
    local ref_count=$(count_sequences "$REF_PROTEIN")
    local query_count=$(count_sequences "$QUERY_PROTEIN")
    local total_hits=$(wc -l < "$OUTPUT_DIR/blastp_result.txt" || echo "0")
    local filtered_hits=$(wc -l < "$OUTPUT_DIR/filtered_result.txt" || echo "0")
    local mapping_count=$(wc -l < "$OUTPUT_DIR/gene_mapping.txt" || echo "0")
    
    cat > "$report_file" << EOF
================================
BLASTP Alignment Analysis Report
================================
Generation time: $(date '+%Y-%m-%d %H:%M:%S')

Configuration:
---------------------------------
Reference proteins: $REF_PROTEIN
Query proteins: $QUERY_PROTEIN
Identity threshold: ${IDENTITY}%
E-value threshold: $EVALUE
Threads: $THREADS

Input statistics:
---------------------------------
Reference sequence count: $ref_count
Query sequence count: $query_count

Alignment results:
---------------------------------
Total hits: $total_hits
Filtered hits: $filtered_hits (identity >= ${IDENTITY}%)
Gene mappings: $mapping_count

Output files:
---------------------------------
1. blastp_result.txt
2. blastp_result_with_header.txt
3. filtered_result.txt
4. filtered_result_with_header.txt
5. gene_mapping.txt

Database:
---------------------------------
BLAST database: $(cat "$OUTPUT_DIR/db_path.txt")

Notes:
---------------------------------
- Mapping format: query_id <TAB> subject_id
- For each query, the subject with the highest bitscore is selected
- Identity threshold is ${IDENTITY}%

================================
EOF

    log_info "Report saved: $report_file"
    cat "$report_file"
}

# Main workflow
main() {
    local config_file=""
    
    while getopts "c:h" opt; do
        case $opt in
            c) config_file=$OPTARG ;;
            h) usage ;;
            *) usage ;;
        esac
    done
    
    if [ -z "$config_file" ]; then
        log_error "Missing configuration file parameter"
        usage
    fi
    
    log_info "=========================================="
    log_info "BLASTP alignment and gene mapping workflow"
    log_info "=========================================="
    
    check_dependencies
    parse_config "$config_file"
    validate_config
    
    mkdir -p "$OUTPUT_DIR"
    log_info "Output directory: $OUTPUT_DIR"
    
    build_blast_db
    run_blastp
    filter_results
    generate_mapping
    generate_report
    
    log_info "=========================================="
    log_info "Analysis completed!"
    log_info "=========================================="
}

# Run main
main "$@"
