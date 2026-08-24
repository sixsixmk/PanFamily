
#!/bin/bash

################################################################################
# KAKS Analysis Pipeline
# Description: Automated Ka/Ks analysis for gene family natural selection pressure
# Author: Bioinformatics Pipeline
# Usage: ./KAKS.sh -c config.yaml
################################################################################

set -e  # Exit on error

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored messages
log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

# Function to show usage
show_usage() {
    echo "Usage: $0 -c <config.yaml>"
    echo ""
    echo "Options:"
    echo "  -c    Configuration file (YAML format)"
    echo "  -h    Show this help message"
    echo ""
    echo "Example:"
    echo "  $0 -c config.yaml"
}

# Parse command line arguments
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

# Check if config file exists
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
# Validate required parameters
################################################################################

log_info "Validating configuration parameters..."

required_params=(
    "config_gene_family_pep"
    "config_gene_family_cds"
    "config_genome_pep"
    "config_genome_cds"
    "config_output_dir"
    "config_processor_file"
)

for param in "${required_params[@]}"; do
    if [ -z "${!param}" ]; then
        log_error "Required parameter '${param#config_}' is missing in config file!"
        exit 1
    fi
done

# Check if input files exist and convert to absolute paths
for file in "$config_gene_family_pep" "$config_gene_family_cds" "$config_genome_pep" "$config_genome_cds"; do
    if [ ! -f "$file" ]; then
        log_error "Input file not found: $file"
        exit 1
    fi
done

# Convert relative paths to absolute paths
config_gene_family_pep=$(realpath "$config_gene_family_pep")
config_gene_family_cds=$(realpath "$config_gene_family_cds")
config_genome_pep=$(realpath "$config_genome_pep")
config_genome_cds=$(realpath "$config_genome_cds")
config_output_dir=$(realpath -m "$config_output_dir")

################################################################################
# Set default parameters
################################################################################

EVALUE=${config_blastp_evalue:-1e-5}
THREADS=${config_blastp_threads:-18}
OUTFMT=${config_blastp_outfmt:-6}
MUSCLE=${config_muscle_path:-muscle}
PARAAT=${config_paraat_path:-ParaAT.pl}

OUTPUT_DIR="$config_output_dir"
WORK_DIR="${OUTPUT_DIR}/workdir"
RESULTS_DIR="${OUTPUT_DIR}/results"

log_info "Output directory: $OUTPUT_DIR"
log_info "BLASTP E-value: $EVALUE"
log_info "Number of threads: $THREADS"

################################################################################
# Create output directories
################################################################################

log_info "Creating output directories..."
mkdir -p "$WORK_DIR"
mkdir -p "$RESULTS_DIR"

cd "$WORK_DIR"

################################################################################
# Build BLAST database
################################################################################

log_info "Building BLAST database..."

GENOME_DB="${WORK_DIR}/ref"
makeblastdb -in "$config_genome_pep" -dbtype prot -out "$GENOME_DB" > makeblastdb.log 2>&1

if [ $? -eq 0 ]; then
    log_info "BLAST database created successfully"
else
    log_error "Failed to create BLAST database"
    exit 1
fi

################################################################################
# Run BLASTP
################################################################################

log_info "Running BLASTP search..."

BLASTP_OUTPUT="${WORK_DIR}/genefamily_vs_genome.blastp"
blastp -query "$config_gene_family_pep" \
       -db "$GENOME_DB" \
       -outfmt "$OUTFMT" \
       -evalue "$EVALUE" \
       -num_threads "$THREADS" \
       -out "$BLASTP_OUTPUT"

if [ $? -eq 0 ]; then
    log_info "BLASTP search completed successfully"
    BLAST_HITS=$(wc -l < "$BLASTP_OUTPUT")
    log_info "Total BLAST hits: $BLAST_HITS"
else
    log_error "BLASTP search failed"
    exit 1
fi

################################################################################
# Filter best hits
################################################################################

log_info "Filtering best hits (lowest E-value per query)..."

BEST_HITS="${WORK_DIR}/best_hits.txt"
awk '{
    if(!seen[$1] || $11 < best[$1]) {
        seen[$1]=1; 
        best[$1]=$11; 
        line[$1]=$0
    }
} 
END {
    for(q in line) print line[q]
}' "$BLASTP_OUTPUT" > "$BEST_HITS"

BEST_HITS_COUNT=$(wc -l < "$BEST_HITS")
log_info "Best hits selected: $BEST_HITS_COUNT"

if [ $BEST_HITS_COUNT -eq 0 ]; then
    log_error "No best hits found! Check your BLASTP parameters."
    exit 1
fi

################################################################################
# Generate gene pairs file
################################################################################

log_info "Generating gene pairs file..."

GENE_PAIRS="${WORK_DIR}/gene_pairs.txt"
awk '{print $1"\t"$2}' "$BEST_HITS" > "$GENE_PAIRS"

log_info "Gene pairs file created: $GENE_PAIRS"

################################################################################
# Combine and clean sequence files
################################################################################

log_info "Preparing sequence files..."

# Combine CDS files
ALL_GENES_CDS_TMP="${WORK_DIR}/all_genes_tmp.cds"
cat "$config_gene_family_cds" "$config_genome_cds" > "$ALL_GENES_CDS_TMP"

# Combine PEP files
ALL_GENES_PEP_TMP="${WORK_DIR}/all_genes_tmp.pep"
cat "$config_gene_family_pep" "$config_genome_pep" > "$ALL_GENES_PEP_TMP"

# Clean sequence headers (remove annotations, keep only ID)
ALL_GENES_CDS="${WORK_DIR}/all_genes.cds"
ALL_GENES_PEP="${WORK_DIR}/all_genes.pep"

awk '/^>/ {print $1; next} {print}' "$ALL_GENES_CDS_TMP" > "$ALL_GENES_CDS"
awk '/^>/ {print $1; next} {print}' "$ALL_GENES_PEP_TMP" > "$ALL_GENES_PEP"

log_info "Sequence files prepared"

################################################################################
# Prepare processor file
################################################################################

log_info "Preparing processor file..."

# Check if processor file exists
if [ ! -f "$config_processor_file" ]; then
    log_error "Processor file not found: $config_processor_file"
    exit 1
fi

# Copy processor file to work directory for ParaAT to use
cp "$config_processor_file" "${WORK_DIR}/proc"
# Use relative path for ParaAT (since we're already in WORK_DIR)
PROC_FILE="proc"
log_info "Using processor file: $config_processor_file (copied to workdir/proc)"

################################################################################
# Run ParaAT for alignment and Ka/Ks calculation
################################################################################

log_info "Running ParaAT for sequence alignment and Ka/Ks calculation..."
log_info "This may take several minutes depending on the number of gene pairs..."

PARAAT_OUTPUT="${WORK_DIR}/paraat_output"

$PARAAT -h "$GENE_PAIRS" \
        -n "$ALL_GENES_CDS" \
        -a "$ALL_GENES_PEP" \
        -p "$PROC_FILE" \
        -m muscle \
        -f axt \
        -k \
        -o "$PARAAT_OUTPUT"

if [ $? -eq 0 ]; then
    log_info "ParaAT analysis completed successfully"
else
    log_error "ParaAT analysis failed"
    exit 1
fi

################################################################################
# Collect and merge Ka/Ks results
################################################################################

log_info "Collecting Ka/Ks results..."

KAKS_RESULTS="${RESULTS_DIR}/all_kaks_results.txt"

# Find all .kaks files and merge them
KAKS_FILES=$(find "$PARAAT_OUTPUT" -name "*.kaks" 2>/dev/null)

if [ -z "$KAKS_FILES" ]; then
    log_error "No Ka/Ks result files found!"
    exit 1
fi

# Get header from first file
FIRST_FILE=$(echo "$KAKS_FILES" | head -1)
head -1 "$FIRST_FILE" > "$KAKS_RESULTS"

# Append all results (skip headers)
find "$PARAAT_OUTPUT" -name "*.kaks" -exec tail -n +2 {} \; >> "$KAKS_RESULTS"

RESULT_COUNT=$(tail -n +2 "$KAKS_RESULTS" | wc -l)
log_info "Total Ka/Ks results: $RESULT_COUNT"

################################################################################
# Create simplified Ka/Ks results
################################################################################

log_info "Creating simplified Ka/Ks results..."

SIMPLIFIED_RESULTS="${RESULTS_DIR}/kaks_simplified.txt"

# Extract TPS gene and Ka/Ks value
awk -F'\t' 'BEGIN {OFS="\t"; print "TPS_Gene", "Ka_Ks"}
NR>1 {
    split($1, genes, "-")
    tps_gene = genes[1]
    ka_ks = $5
    print tps_gene, ka_ks
}' "$KAKS_RESULTS" > "$SIMPLIFIED_RESULTS"

SIMPLIFIED_COUNT=$(tail -n +2 "$SIMPLIFIED_RESULTS" | wc -l)
log_info "Simplified results created: $SIMPLIFIED_COUNT genes"
log_info "Simplified results saved to: $SIMPLIFIED_RESULTS"

################################################################################
# Generate summary statistics
################################################################################

log_info "Generating summary statistics..."

SUMMARY_FILE="${RESULTS_DIR}/summary.txt"

cat > "$SUMMARY_FILE" << EOF
================================================================================
Ka/Ks Analysis Summary
================================================================================
Analysis Date: $(date)
Configuration File: $CONFIG_FILE

Input Files:
  - Gene Family PEP: $config_gene_family_pep
  - Gene Family CDS: $config_gene_family_cds
  - Genome PEP: $config_genome_pep
  - Genome CDS: $config_genome_cds

BLASTP Parameters:
  - E-value: $EVALUE
  - Threads: $THREADS
  - Output format: $OUTFMT

Results:
  - Total BLAST hits: $BLAST_HITS
  - Best hits (gene pairs): $BEST_HITS_COUNT
  - Ka/Ks calculations: $RESULT_COUNT

Output Files:
  - Ka/Ks results: $KAKS_RESULTS
  - Visualization CSV: $KAKS_VIS_CSV
  - Gene pairs: $GENE_PAIRS
  - BLAST results: $BLASTP_OUTPUT
  - ParaAT output: $PARAAT_OUTPUT

Selection Pressure Summary:
EOF

# Calculate selection pressure statistics
awk -F'\t' 'NR>1 {
    ka_ks = $5
    if (ka_ks != "NA" && ka_ks > 0) {
        if (ka_ks < 1) purifying++
        else if (ka_ks > 1) positive++
        else neutral++
        total++
        sum += ka_ks
    }
}
END {
    if (total > 0) {
        printf "  - Purifying selection (Ka/Ks < 1): %d (%.2f%%)\n", purifying, purifying/total*100
        printf "  - Neutral selection (Ka/Ks = 1): %d (%.2f%%)\n", neutral, neutral/total*100
        printf "  - Positive selection (Ka/Ks > 1): %d (%.2f%%)\n", positive, positive/total*100
        printf "  - Average Ka/Ks: %.4f\n", sum/total
    }
}' "$KAKS_RESULTS" >> "$SUMMARY_FILE"

echo "================================================================================" >> "$SUMMARY_FILE"

log_info "Summary saved to: $SUMMARY_FILE"

################################################################################
# Copy important files to results directory
################################################################################

log_info "Organizing results..."

cp "$GENE_PAIRS" "$RESULTS_DIR/"
cp "$BLASTP_OUTPUT" "$RESULTS_DIR/"
cp "$BEST_HITS" "$RESULTS_DIR/"

# Create a simplified results table
SIMPLE_RESULTS="${RESULTS_DIR}/kaks_simple.txt"
awk -F'\t' 'NR==1 {print "Gene_Pair\tMethod\tKa\tKs\tKa_Ks\tP_Value"} 
            NR>1 {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' "$KAKS_RESULTS" > "$SIMPLE_RESULTS"

log_info "Simplified results saved to: $SIMPLE_RESULTS"

################################################################################
# Generate visualization-ready CSV file (kaks.csv)
################################################################################

log_info "Generating visualization-ready kaks.csv..."

KAKS_VIS_CSV="${RESULTS_DIR}/kaks.csv"

# Extract Gene (from Gene_Pair), Ka, Ks, KaKs for visualization module
awk -F'\t' 'BEGIN {OFS=","; print "Gene,Ka,Ks,KaKs"}
NR>1 {
    # Extract gene name from Gene_Pair (format: GeneA-GeneB)
    split($1, genes, "-")
    gene = genes[1]
    ka = $3
    ks = $4
    kaks = $5
    # Skip invalid values
    if (ka != "NA" && ks != "NA" && kaks != "NA" && ks > 0) {
        print gene, ka, ks, kaks
    }
}' "$KAKS_RESULTS" > "$KAKS_VIS_CSV"

KAKS_VIS_COUNT=$(tail -n +2 "$KAKS_VIS_CSV" | wc -l)
log_info "Visualization CSV created: $KAKS_VIS_CSV ($KAKS_VIS_COUNT entries)"

################################################################################
# Completion
################################################################################

echo ""
log_info "=========================================="
log_info "Ka/Ks Analysis Pipeline Completed!"
log_info "=========================================="
echo ""
log_info "Results directory: $RESULTS_DIR"
log_info "Main results file: $KAKS_RESULTS"
log_info "Simplified results: $SIMPLIFIED_RESULTS"
log_info "Visualization CSV: $KAKS_VIS_CSV"
log_info "Summary file: $SUMMARY_FILE"
echo ""
log_info "Please check the summary file for analysis statistics."
echo ""

# Display summary
cat "$SUMMARY_FILE"

exit 0

