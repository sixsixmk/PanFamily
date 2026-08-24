#!/bin/bash
set -e 
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

log_info() {
    echo -e "${GREEN}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

# Usage
show_usage() {
cat << EOF
Usage: $0 -c <config.yaml>

Multi-strain Ka/Ks analysis pipeline

Options:
  -c <config.yaml>    Configuration file (YAML format)
  -h                  Show this help message

Example:
  $0 -c pankaks.yaml
EOF
}

CONFIG_FILE=""

while getopts "c:h" opt; do
    case $opt in
        c)
            CONFIG_FILE="$OPTARG"
            ;;
        h)
            show_usage
            exit 0
            ;;
        \?)
            log_error "Invalid option: -$OPTARG"
            show_usage
            exit 1
            ;;
    esac
done

# Check if config file is provided
if [ -z "$CONFIG_FILE" ]; then
    log_error "No configuration file provided!"
    show_usage
    exit 1
fi

if [ ! -f "$CONFIG_FILE" ]; then
    log_error "Configuration file '$CONFIG_FILE' not found!"
    exit 1
fi

log_info "Reading configuration from: $CONFIG_FILE"

################################################################################
# Parse YAML configuration file
################################################################################

parse_yaml() {
    local prefix=$2
    local s='[[:space:]]*'
    local w='[a-zA-Z0-9_]*'
    local fs=$(echo @|tr @ '\034')
    sed -ne "s|^\($s\):|\1|" \
         -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
         -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p" $1 |
    awk -F$fs '{
        indent = length($1)/2;
        vname[indent] = $2;
        for (i in vname) {if (i > indent) {delete vname[i]}}
        if (length($3) > 0) {
            vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
            printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
        }
    }'
}

# Parse YAML and export variables
eval $(parse_yaml "$CONFIG_FILE" "config_")

################################################################################
# Validate configuration
################################################################################

log_info "Validating configuration..."

# Required parameters
required_params=(
    "config_gene_family_pep"
    "config_gene_family_cds"
    "config_strains_pep_dir"
    "config_strains_cds_dir"
    "config_strain_list"
    "config_output_dir"
)

for param in "${required_params[@]}"; do
    if [ -z "${!param}" ]; then
        log_error "Required parameter '${param#config_}' is missing in config file!"
        exit 1
    fi
done

# Check if directories and files exist
for file in "$config_gene_family_pep" "$config_gene_family_cds"; do
    if [ ! -f "$file" ]; then
        log_error "File not found: $file"
        exit 1
    fi
done

for dir in "$config_strains_pep_dir" "$config_strains_cds_dir"; do
    if [ ! -d "$dir" ]; then
        log_error "Directory not found: $dir"
        exit 1
    fi
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ -z "$config_processor_file" ]; then
    config_processor_file="${SCRIPT_DIR}/proc"
fi

if [ -z "$config_kaks_script" ]; then
    config_kaks_script="${SCRIPT_DIR}/KAKS.sh"
fi

if [ ! -f "$config_processor_file" ]; then
    log_error "Processor file not found: $config_processor_file"
    log_error "Please create a 'proc' file in the same directory as this script"
    exit 1
fi

if [ ! -f "$config_kaks_script" ]; then
    log_error "KAKS.sh script not found: $config_kaks_script"
    log_error "Expected location: ${SCRIPT_DIR}/../KAKS.sh"
    exit 1
fi

config_gene_family_pep=$(realpath "$config_gene_family_pep")
config_gene_family_cds=$(realpath "$config_gene_family_cds")
config_strains_pep_dir=$(realpath "$config_strains_pep_dir")
config_strains_cds_dir=$(realpath "$config_strains_cds_dir")
config_processor_file=$(realpath "$config_processor_file")
config_kaks_script=$(realpath "$config_kaks_script")
config_output_dir=$(realpath -m "$config_output_dir")

log_info "Using processor file: $config_processor_file"
log_info "Using KAKS.sh script: $config_kaks_script"

EVALUE=${config_blastp_evalue:-1e-5}
THREADS=${config_blastp_threads:-18}

log_info "Configuration validated successfully"

################################################################################
# Create output directories
################################################################################

log_info "Creating output directory structure..."

mkdir -p "$config_output_dir"
mkdir -p "$config_output_dir/all_strains"
mkdir -p "$config_output_dir/merged_matrix"

################################################################################
# Process strains
################################################################################

if [ -f "$config_strain_list" ]; then
    # Read from file (one strain per line)
    log_info "Reading strain list from file: $config_strain_list"
    mapfile -t STRAINS < <(grep -v '^#' "$config_strain_list" | grep -v '^$' | sed 's/[[:space:]]*$//')
else
    log_info "Parsing comma-separated strain list"
    IFS=',' read -ra STRAINS <<< "$config_strain_list"
fi

STRAIN_COUNT=${#STRAINS[@]}

log_info "Processing $STRAIN_COUNT strains: ${STRAINS[*]}"

FAILED_STRAINS=()
SUCCESSFUL_STRAINS=()

for strain in "${STRAINS[@]}"; do
    strain=$(echo "$strain" | xargs)  
    
    log_info "=========================================="
    log_info "Processing strain: $strain"
    log_info "=========================================="
    
    # Check if strain files exist
    STRAIN_PEP="${config_strains_pep_dir}/${strain}.pep"
    STRAIN_CDS="${config_strains_cds_dir}/${strain}.cds"
    
    if [ ! -f "$STRAIN_PEP" ]; then
        log_error "Strain PEP file not found: $STRAIN_PEP"
        FAILED_STRAINS+=("$strain")
        continue
    fi
    
    if [ ! -f "$STRAIN_CDS" ]; then
        log_error "Strain CDS file not found: $STRAIN_CDS"
        FAILED_STRAINS+=("$strain")
        continue
    fi
    
    # Create strain output directory
    STRAIN_OUTPUT="${config_output_dir}/${strain}"
    mkdir -p "$STRAIN_OUTPUT"
    
    # Generate strain-specific config file
    STRAIN_CONFIG="${STRAIN_OUTPUT}/config_${strain}.yaml"
    
    cat > "$STRAIN_CONFIG" << EOF
# Auto-generated config for strain: $strain
gene_family_pep: "$config_gene_family_pep"
gene_family_cds: "$config_gene_family_cds"
genome_pep: "$STRAIN_PEP"
genome_cds: "$STRAIN_CDS"
output_dir: "$STRAIN_OUTPUT"
blastp:
  evalue: $EVALUE
  threads: $THREADS
processor_file: "$config_processor_file"
EOF
    
    log_info "Running KAKS.sh for strain: $strain"
    
    # Run KAKS.sh
    if bash "$config_kaks_script" -c "$STRAIN_CONFIG"; then
        log_info "KAKS.sh completed successfully for strain: $strain"
        SUCCESSFUL_STRAINS+=("$strain")
        
        # Copy and rename results to all_strains directory
        if [ -f "${STRAIN_OUTPUT}/results/all_kaks_results.txt" ]; then
            cp "${STRAIN_OUTPUT}/results/all_kaks_results.txt" \
               "${config_output_dir}/all_strains/${strain}_all_kaks.txt"
        fi
        
        if [ -f "${STRAIN_OUTPUT}/results/kaks_simplified.txt" ]; then
            cp "${STRAIN_OUTPUT}/results/kaks_simplified.txt" \
               "${config_output_dir}/all_strains/${strain}_kaks_simplified.txt"
        fi
        
        # Copy visualization CSV
        if [ -f "${STRAIN_OUTPUT}/results/kaks.csv" ]; then
            cp "${STRAIN_OUTPUT}/results/kaks.csv" \
               "${config_output_dir}/all_strains/${strain}_kaks.csv"
        fi
        
    else
        log_error "KAKS.sh failed for strain: $strain"
        FAILED_STRAINS+=("$strain")
    fi
done

################################################################################
# Merge results (long format with strain column)
################################################################################

log_info "=========================================="
log_info "Merging results..."
log_info "=========================================="

if [ ${#SUCCESSFUL_STRAINS[@]} -eq 0 ]; then
    log_error "No strains were processed successfully!"
    exit 1
fi

# Merge simplified Ka/Ks results (long format: Gene, Ka_Ks, Strain)
MERGED_KAKS_FILE="${config_output_dir}/merged_matrix/merged_kaks_simplified.txt"

# Create header
echo -e "Gene\tKa_Ks\tStrain" > "$MERGED_KAKS_FILE"

# Append all strain results with strain name in third column
for strain in "${SUCCESSFUL_STRAINS[@]}"; do
    STRAIN_SIMPLIFIED="${config_output_dir}/all_strains/${strain}_kaks_simplified.txt"
    if [ -f "$STRAIN_SIMPLIFIED" ]; then
        # Skip header, add strain column
        tail -n +2 "$STRAIN_SIMPLIFIED" | awk -v strain="$strain" -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, strain}' >> "$MERGED_KAKS_FILE"
    fi
done

MERGED_COUNT=$(($(wc -l < "$MERGED_KAKS_FILE") - 1))
log_info "Merged Ka/Ks results: $MERGED_COUNT entries"
log_info "Output file: $MERGED_KAKS_FILE"

################################################################################
# Merge visualization CSV files (kaks.csv)
################################################################################

log_info "Merging visualization CSV files..."

MERGED_KAKS_CSV="${config_output_dir}/merged_matrix/kaks.csv"

# Create header
echo "Gene,Ka,Ks,KaKs" > "$MERGED_KAKS_CSV"

# Append all strain kaks.csv (skip headers)
for strain in "${SUCCESSFUL_STRAINS[@]}"; do
    STRAIN_KAKS_CSV="${config_output_dir}/all_strains/${strain}_kaks.csv"
    if [ -f "$STRAIN_KAKS_CSV" ]; then
        tail -n +2 "$STRAIN_KAKS_CSV" >> "$MERGED_KAKS_CSV"
    fi
done

MERGED_KAKS_COUNT=$(tail -n +2 "$MERGED_KAKS_CSV" | wc -l)
log_info "Merged visualization CSV created: $MERGED_KAKS_CSV ($MERGED_KAKS_COUNT entries)"

################################################################################
# Generate summary
################################################################################

log_info "Generating summary..."

SUMMARY_FILE="${config_output_dir}/merged_matrix/summary.txt"

cat > "$SUMMARY_FILE" << EOF
================================================================================
PanKAKS Multi-Strain Analysis Summary
================================================================================
Analysis Date: $(date)
Configuration File: $CONFIG_FILE

Input Files:
  - Gene Family PEP: $config_gene_family_pep
  - Gene Family CDS: $config_gene_family_cds
  - Strains PEP Directory: $config_strains_pep_dir
  - Strains CDS Directory: $config_strains_cds_dir

Strains Processed:
  - Total strains: $STRAIN_COUNT
  - Successful: ${#SUCCESSFUL_STRAINS[@]}
  - Failed: ${#FAILED_STRAINS[@]}

Successful Strains:
$(for s in "${SUCCESSFUL_STRAINS[@]}"; do echo "  - $s"; done)

EOF

if [ ${#FAILED_STRAINS[@]} -gt 0 ]; then
    cat >> "$SUMMARY_FILE" << EOF
Failed Strains:
$(for s in "${FAILED_STRAINS[@]}"; do echo "  - $s"; done)

EOF
fi

cat >> "$SUMMARY_FILE" << EOF
Output Files:
  - Merged Ka/Ks file: $MERGED_KAKS_FILE
  - Visualization CSV: $MERGED_KAKS_CSV
  - Individual results: ${config_output_dir}/all_strains/
  - Strain directories: ${config_output_dir}/[strain_name]/

================================================================================
EOF

log_info "Summary saved to: $SUMMARY_FILE"

################################################################################
# Completion
################################################################################

echo ""
log_info "=========================================="
log_info "PanKAKS Analysis Completed!"
log_info "=========================================="
echo ""
log_info "Results directory: $config_output_dir"
log_info "Merged Ka/Ks file: $MERGED_KAKS_FILE"
log_info "Visualization CSV: $MERGED_KAKS_CSV"
log_info "Summary file: $SUMMARY_FILE"
echo ""
log_info "Successful strains: ${#SUCCESSFUL_STRAINS[@]}/${STRAIN_COUNT}"
echo ""
cat "$SUMMARY_FILE"

exit 0
