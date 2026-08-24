#!/bin/bash
set -e
set -u
set -o pipefail

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Logging functions
log_info() {
    echo -e "${GREEN}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

# Usage
usage() {
    cat << EOF
Usage: $0 -c <config_file>

Required argument:
    -c    Path to config file

Example config (config.txt):
    REF_GENOME=/path/to/reference.fa
    GENEPRED=/path/to/reference_refGene.txt
    VCF_FILE=/path/to/all_merge.vcf
    AVINPUT=/path/to/all_merge.avinput
    GENE_LIST=/path/to/gene_family.list
    OUTPUT_DIR=/path/to/output
    BUILD_VER=B73
    NEARGENE=2000
    ID_STRAIN=/path/to/id_strain.txt
    IDMAP=/path/to/Idmap.txt

Description:
    - REF_GENOME: Path to reference genome FASTA (required)
    - GENEPRED: Path to genePred annotation file (required)
    - VCF_FILE: Structural variant VCF (optional if AVINPUT is provided)
    - AVINPUT: ANNOVAR input file (required)
    - GENE_LIST: Gene family list file, 1 gene ID per line (required)
    - OUTPUT_DIR: Output directory (required)
    - BUILD_VER: Genome version name, e.g. B73 (required)
    - NEARGENE: Gene flank distance, default 2000 bp (optional)
    - ID_STRAIN: Gene ID to strain mapping file, tab-delimited (optional)
    - IDMAP: Gene ID renaming map file, tab-delimited (optional)

EOF
    exit 1
}

# Check required tools
check_dependencies() {
    local tools=("annotate_variation.pl" "retrieve_seq_from_fasta.pl")
    for tool in "${tools[@]}"; do
        if ! command -v $tool &> /dev/null; then
            log_error "Required tool not found: $tool"
            log_error "Please ensure ANNOVAR is correctly installed and added to PATH"
            exit 1
        fi
    done
    log_info "Dependency check passed"
}

# Parse config
parse_config() {
    local config_file=$1
    
    if [ ! -f "$config_file" ]; then
        log_error "Config file does not exist: $config_file"
        exit 1
    fi
    
    log_info "Parsing config: $config_file"
    
    # Read config
    while IFS='=' read -r key value; do
        # Skip comments and empty lines
        [[ $key =~ ^#.*$ ]] && continue
        [[ -z $key ]] && continue
        
        # Remove surrounding spaces
        key=$(echo "$key" | xargs)
        value=$(echo "$value" | xargs)
        
        case $key in
            REF_GENOME) REF_GENOME=$value ;;
            GENEPRED) GENEPRED=$value ;;
            VCF_FILE) VCF_FILE=$value ;;
            AVINPUT) AVINPUT=$value ;;
            GENE_LIST) GENE_LIST=$value ;;
            OUTPUT_DIR) OUTPUT_DIR=$value ;;
            BUILD_VER) BUILD_VER=$value ;;
            NEARGENE) NEARGENE=$value ;;
            ID_STRAIN) ID_STRAIN=$value ;;
            IDMAP) IDMAP=$value ;;
        esac
    done < "$config_file"
    
    # Default values
    NEARGENE=${NEARGENE:-2000}
}

# Validate config
validate_config() {
    log_info "Validating configuration..."
    
    local errors=0
    
    if [ -z "${REF_GENOME:-}" ] || [ ! -f "$REF_GENOME" ]; then
        log_error "Reference genome not found: ${REF_GENOME:-unset}"
        ((errors++))
    fi
    
    if [ -z "${GENEPRED:-}" ] || [ ! -f "$GENEPRED" ]; then
        log_error "GenePred file not found: ${GENEPRED:-unset}"
        ((errors++))
    fi
    
    if [ -z "${AVINPUT:-}" ] || [ ! -f "$AVINPUT" ]; then
        log_error "AVINPUT file not found: ${AVINPUT:-unset}"
        ((errors++))
    fi
    
    if [ -z "${GENE_LIST:-}" ] || [ ! -f "$GENE_LIST" ]; then
        log_error "Gene list not found: ${GENE_LIST:-unset}"
        ((errors++))
    fi
    
    if [ -z "${OUTPUT_DIR:-}" ]; then
        log_error "Output directory not set"
        ((errors++))
    fi
    
    if [ -z "${BUILD_VER:-}" ]; then
        log_error "Genome build name not set"
        ((errors++))
    fi
    
    # Optional parameters validation
    if [ -n "${ID_STRAIN:-}" ] && [ ! -f "$ID_STRAIN" ]; then
        log_error "ID_STRAIN file not found: $ID_STRAIN"
        ((errors++))
    fi
    
    if [ -n "${IDMAP:-}" ] && [ ! -f "$IDMAP" ]; then
        log_error "IDMAP file not found: $IDMAP"
        ((errors++))
    fi
    
    if [ $errors -gt 0 ]; then
        log_error "Configuration validation failed with $errors error(s)"
        exit 1
    fi
    
    log_info "Configuration validation passed"
}

# Build ANNOVAR database
build_annovar_db() {
    log_info "========== Step 1: Build ANNOVAR database =========="
    
    local db_dir="$OUTPUT_DIR/dbase"
    mkdir -p "$db_dir"
    
    log_info "Database directory: $db_dir"
    
    # Copy reference genome
    log_info "Preparing reference genome..."
    if [ ! -f "$db_dir/${BUILD_VER}.fa" ]; then
        cp "$REF_GENOME" "$db_dir/${BUILD_VER}.fa"
    fi
    
    # Copy GenePred
    log_info "Preparing GenePred file..."
    cp "$GENEPRED" "$db_dir/${BUILD_VER}_refGene.txt"
    
    # Generate refGene.fa
    log_info "Generating gene sequence file ${BUILD_VER}_refGene.fa ..."
    cd "$db_dir"
    retrieve_seq_from_fasta.pl \
        --format refGene \
        --seqfile "${BUILD_VER}.fa" \
        "${BUILD_VER}_refGene.txt" \
        --out "${BUILD_VER}_refGene.fa"
    cd - > /dev/null
    
    log_info "Database build completed"
    echo "${BUILD_VER}_refGene.txt" > "$OUTPUT_DIR/db_files.txt"
    echo "${BUILD_VER}_refGene.fa" >> "$OUTPUT_DIR/db_files.txt"
}

# Run ANNOVAR
run_annovar() {
    log_info "========== Step 2: Run ANNOVAR annotation =========="
    
    local db_dir="$OUTPUT_DIR/dbase"
    local output_prefix="$OUTPUT_DIR/sv_annotation"
    
    log_info "Input AVINPUT: $AVINPUT"
    log_info "Database dir: $db_dir"
    log_info "Output prefix: $output_prefix"
    log_info "Near gene distance: ${NEARGENE}bp"
    
    annotate_variation.pl \
        -geneanno \
        --neargene $NEARGENE \
        -buildver $BUILD_VER \
        -dbtype refGene \
        -outfile "$output_prefix" \
        -exonsort \
        "$AVINPUT" \
        "$db_dir"
    
    log_info "ANNOVAR annotation completed"
    log_info "Output file: ${output_prefix}.variant_function"
}

# Extract target genes
extract_target_genes() {
    log_info "========== Step 3: Extract structural variants for target genes =========="
    
    local annotation_file="$OUTPUT_DIR/sv_annotation.variant_function"
    local output_file="$OUTPUT_DIR/target_genes_sv.txt"
    
    if [ ! -f "$annotation_file" ]; then
        log_error "Annotation file not found: $annotation_file"
        exit 1
    fi
    
    log_info "Gene list: $GENE_LIST"
    log_info "Annotation file: $annotation_file"
    
    local gene_count=$(wc -l < "$GENE_LIST")
    log_info "Number of genes to extract: $gene_count"
    
    grep -F -f "$GENE_LIST" "$annotation_file" > "$output_file" || true
    
    local match_count=$(wc -l < "$output_file")
    log_info "Number of matched SV records: $match_count"
    
    if [ $match_count -eq 0 ]; then
        log_warn "No matched structural variant records found"
    else
        log_info "Result saved to: $output_file"
    fi
}

# Generate report
generate_report() {
    log_info "========== Generate analysis report =========="
    
    local report_file="$OUTPUT_DIR/analysis_report.txt"
    
    cat > "$report_file" << EOF
================================
Structural Variant Annotation Report
================================
Generated at: $(date '+%Y-%m-%d %H:%M:%S')

Configuration:
---------------------------------
Reference genome: $REF_GENOME
GenePred file: $GENEPRED
AVINPUT file: $AVINPUT
Gene list: $GENE_LIST
Genome build: $BUILD_VER
Near gene distance: ${NEARGENE}bp

Input statistics:
---------------------------------
Number of SVs: $(wc -l < "$AVINPUT")
Number of target genes: $(wc -l < "$GENE_LIST")

Output files:
---------------------------------
Annotation result: $OUTPUT_DIR/sv_annotation.variant_function
Target gene SV: $OUTPUT_DIR/target_genes_sv.txt
Matched SV count: $(wc -l < "$OUTPUT_DIR/target_genes_sv.txt")

Database files:
---------------------------------
$(cat "$OUTPUT_DIR/db_files.txt")

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
        log_error "Missing config file argument"
        usage
    fi
    
    log_info "=========================================="
    log_info "Structural Variant ANNOVAR Annotation Workflow"
    log_info "=========================================="
    
    check_dependencies
    parse_config "$config_file"
    validate_config
    
    mkdir -p "$OUTPUT_DIR"
    log_info "Output dir: $OUTPUT_DIR"
    
    build_annovar_db
    run_annovar
    extract_target_genes
    generate_report
    
    log_info "=========================================="
    log_info "Analysis completed!"
    log_info "=========================================="
}

main "$@"
