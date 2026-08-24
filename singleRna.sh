#!/bin/bash



set -euo pipefail  

PROGRAM_NAME="PanFamily RNAseq"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log_info() {
    echo -e "${GREEN}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2
}

log_step() {
    echo -e "${BLUE}[STEP]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

print_banner() {
    echo "============================================================================"
    echo "                           $PROGRAM_NAME "
    echo "============================================================================"
}

print_usage() {
    cat << EOF
USAGE:
    $0 [OPTIONS]

Configuration Methods (choose one):
    Method 1: YAML configuration file (recommended)
    -c <FILE>               YAML configuration file (contains all parameters)

    Method 2: Command line parameters
    -i <FILE> [<FILE>]      Input sequencing files (fastq/fastq.gz)
                              - Single-end (SE): provide only one file
                              - Paired-end (PE): provide two files
    -s <NAME>               Sample name
    -o <DIR>                Output directory

Reference Genome Related (required for command line mode):
    -g <FILE>               Reference genome FASTA file
    -a <FILE>               GTF annotation file

Data Processing Parameters (optional):
    --qc-before             Run quality control before trimming (default: enabled)
    --qc-after              Run quality control after trimming (default: enabled)
    --no-trimming           Skip trimming step
    --no-qc                 Skip all quality control steps

Alignment/Workflow Parameters (optional):
    -p <INT>                Number of threads (default: 8)
    --keep-intermediate     Keep intermediate files
    -h, --help              Show help information
    -v, --version           Show version information

OUTPUT STRUCTURE:
    output_dir/
    ├── data_processing/   Data processing results
    │   ├── qc_before/     Quality control before trimming
    │   ├── trimmed/       Trimmed reads
    │   └── qc_after/      Quality control after trimming
    ├── hisat2/            Alignment results (BAM files)
    ├── featurecounts/     featureCounts quantification results
    └── logs/              Log files

YAML Configuration File Example:
    # Complete YAML configuration (recommended method)
    analysis:
      sample_name: "sample_01"
      threads: 8
      keep_intermediate: false

    input:
      read1: "/path/to/sample_R1.fastq.gz"
      read2: "/path/to/sample_R2.fastq.gz"
      genome_file: "/path/to/genome.fa"
      annotation_file: "/path/to/annotation.gtf"

    output:
      output_dir: "/path/to/output_directory"

    data_processing:
      qc_before: true
      qc_after: true
      trimming: true
      trim_quality: 25
      trim_length: 35
      fastqc_options: ""

EXAMPLES:
    # Use complete YAML configuration file (recommended)
    $0 -c config.yaml

    # YAML + command line parameters
    $0 -c config.yaml -p 16

    # Traditional command line mode (still supported)
    $0 -i sample_R1.fq.gz sample_R2.fq.gz -s sample_pe \\
       -g /path/to/genome.fa -a annotation.gtf -o results_pe \\
       --qc-before --qc-after

EOF
}

# Complete function for parsing YAML configuration file
parse_yaml() {
    local yaml_file="$1"
    if [[ ! -f "$yaml_file" ]]; then
        log_error "YAML configuration file does not exist: $yaml_file"
        exit 1
    fi

    log_info "Parsing YAML configuration file: $yaml_file"

    # Read YAML configuration (simple parsing supporting nested structures)
    while IFS=': ' read -r key value; do
        # Remove spaces and comments, and handle quotes
        key=$(echo "$key" | sed 's/^[[:space:]]*//' | sed 's/[[:space:]]*$//')
        value=$(echo "$value" | sed 's/^[[:space:]]*//' | sed 's/[[:space:]]*$//' | sed 's/#.*//' | sed 's/^"//;s/"$//')

        # Skip empty values and comment lines
        [[ -z "$value" || "$key" =~ ^# ]] && continue

        case "$key" in
            # Basic analysis parameters
            "sample_name")
                SAMPLE_NAME="$value"
                ;;
            "threads")
                THREADS="$value"
                ;;
            "keep_intermediate")
                KEEP_INTERMEDIATE=$(echo "$value" | tr '[:upper:]' '[:lower:]')
                ;;
            # Input files
            "read1")
                READ1="$value"
                ;;
            "read2")
                READ2="$value"
                ;;
            "genome_file")
                GENOME_FILE="$value"
                ;;
            "annotation_file")
                ANNOTATION_FILE="$value"
                ;;
            # Output directory
            "output_dir")
                OUTPUT_DIR="$value"
                ;;
            # Data processing parameters
            "qc_before")
                QC_BEFORE=$(echo "$value" | tr '[:upper:]' '[:lower:]')
                ;;
            "qc_after")
                QC_AFTER=$(echo "$value" | tr '[:upper:]' '[:lower:]')
                ;;
            "trimming")
                TRIMMING=$(echo "$value" | tr '[:upper:]' '[:lower:]')
                ;;
            "trim_quality")
                TRIM_QUALITY="$value"
                ;;
            "trim_length")
                TRIM_LENGTH="$value"
                ;;
            "fastqc_options")
                FASTQC_OPTIONS="$value"
                ;;
        esac
    done < <(grep -E "^\s*(sample_name|threads|keep_intermediate|read1|read2|genome_file|annotation_file|output_dir|qc_before|qc_after|trimming|trim_quality|trim_length|fastqc_options):" "$yaml_file")

    log_info "YAML configuration parsing completed"
}

# Check required software tools
check_dependencies() {
    local required_tools=("samtools")
    local missing_tools=()

    # Add tools based on data processing options
    if [[ "$QC_BEFORE" == "true" ]] || [[ "$QC_AFTER" == "true" ]]; then
        required_tools+=("fastqc")
    fi

    if [[ "$TRIMMING" == "true" ]]; then
        required_tools+=("trim_galore")
    fi

    # Tools required for featureCounts
    required_tools+=("hisat2" "featureCounts")

    for tool in "${required_tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            missing_tools+=("$tool")
        fi
    done

    if [ ${#missing_tools[@]} -ne 0 ]; then
        log_error "Missing required tools: ${missing_tools[*]}"
        log_error "Please install these tools before running analysis"
        exit 1
    fi

    log_info "Dependency tools check passed ✓"
}

# Create output directory structure
create_directories() {
    log_info "Creating output directory structure"

    mkdir -p "$OUTPUT_DIR"/{data_processing,hisat2,logs,featurecounts}
    mkdir -p "$OUTPUT_DIR/data_processing"/{qc_before,trimmed,qc_after}

    # Set log file
    LOG_FILE="$OUTPUT_DIR/logs/${SAMPLE_NAME}_$(date +%Y%m%d_%H%M%S).log"

    log_info "Output directory: $OUTPUT_DIR"
    log_info "Log file: $LOG_FILE"
}

# Quality control before trimming
run_qc_before() {
    log_info "Running quality control before trimming"

    local qc_dir="$OUTPUT_DIR/data_processing/qc_before"

    if [ "$DATA_TYPE" = "PE" ]; then
        fastqc $FASTQC_OPTIONS -t "$THREADS" -o "$qc_dir" "$READ1" "$READ2" 2>> "$LOG_FILE"
    else
        fastqc $FASTQC_OPTIONS -t "$THREADS" -o "$qc_dir" "$READ1" 2>> "$LOG_FILE"
    fi

    log_info "Quality control before trimming completed ✓"
}

# Quality control after trimming
run_qc_after() {
    log_info "Running quality control after trimming"

    local qc_dir="$OUTPUT_DIR/data_processing/qc_after"

    if [ "$DATA_TYPE" = "PE" ]; then
        fastqc $FASTQC_OPTIONS -t "$THREADS" -o "$qc_dir" "$CLEAN_READ1" "$CLEAN_READ2" 2>> "$LOG_FILE"
    else
        fastqc $FASTQC_OPTIONS -t "$THREADS" -o "$qc_dir" "$CLEAN_READ1" 2>> "$LOG_FILE"
    fi

    log_info "Quality control after trimming completed ✓"
}

# Sequence trimming
run_trimming() {
    log_info "Using Trim Galore for quality trimming"

    local trim_dir="$OUTPUT_DIR/data_processing/trimmed"

    if [ "$DATA_TYPE" = "PE" ]; then
        trim_galore \
            -j "$THREADS" \
            -q "$TRIM_QUALITY" \
            --phred33 \
            --length "$TRIM_LENGTH" \
            --paired \
            -o "$trim_dir" \
            "$READ1" "$READ2" 2>> "$LOG_FILE"

        local base1=$(basename "$READ1" .fastq.gz)
        local base2=$(basename "$READ2" .fastq.gz)
        base1=$(basename "$base1" .fq.gz)
        base2=$(basename "$base2" .fq.gz)

        CLEAN_READ1="$trim_dir/${base1}_val_1.fq.gz"
        CLEAN_READ2="$trim_dir/${base2}_val_2.fq.gz"
    else
        trim_galore \
            -j "$THREADS" \
            -q "$TRIM_QUALITY" \
            --phred33 \
            --length "$TRIM_LENGTH" \
            -o "$trim_dir" \
            "$READ1" 2>> "$LOG_FILE"

        local base=$(basename "$READ1" .fastq.gz)
        base=$(basename "$base" .fq.gz)

        CLEAN_READ1="$trim_dir/${base}_trimmed.fq.gz"
        CLEAN_READ2=""
    fi

    log_info "Sequence trimming completed ✓"
}

# Main data processing workflow
run_data_processing() {
    log_step "Data processing workflow"

    # Initialize cleaned reads as original reads
    CLEAN_READ1="$READ1"
    CLEAN_READ2="$READ2"

    # 1. Quality control before trimming
    if [[ "$QC_BEFORE" == "true" ]]; then
        run_qc_before
    fi

    # 2. Sequence trimming
    if [[ "$TRIMMING" == "true" ]]; then
        run_trimming
    fi

    # 3. Quality control after trimming
    if [[ "$QC_AFTER" == "true" ]]; then
        run_qc_after
    fi

    log_info "Data processing workflow completed ✓"
}

# Prepare HISAT2 index (build or use existing)
prepare_hisat2_index() {
    log_step "Preparing HISAT2 index"

    # Check if genome file provided with -g exists and use it to build index
    if [[ ! -n "$GENOME_FILE" ]]; then
        log_error "Missing genome FASTA file required for HISAT2 index (-g)"
        exit 1
    fi

    local index_dir="$OUTPUT_DIR/hisat2"
    local index_base="$index_dir/genome_index"

    # Check if index already exists
    if [[ ! -f "${index_base}.1.ht2" ]]; then
        log_info "HISAT2 index does not exist, starting to build using genome file: $GENOME_FILE"
        hisat2-build -p "$THREADS" "$GENOME_FILE" "$index_base" 2>> "$LOG_FILE"
    else
        log_info "HISAT2 index already exists, skipping build"
    fi

    HISAT2_INDEX="$index_base"
}

# Sequence alignment (for alignment-based methods)
run_alignment() {
    log_step "Sequence alignment"

    prepare_hisat2_index

    log_info "Using HISAT2 for sequence alignment"
    local sam_file="$OUTPUT_DIR/hisat2/${SAMPLE_NAME}.sam"
    local sorted_bam="$OUTPUT_DIR/hisat2/${SAMPLE_NAME}.sorted.bam"

    if [ "$DATA_TYPE" = "PE" ]; then
        hisat2 -p "$THREADS" \
               -x "$HISAT2_INDEX" \
               -1 "$CLEAN_READ1" \
               -2 "$CLEAN_READ2" \
               -S "$sam_file" 2>> "$LOG_FILE"
    else
        hisat2 -p "$THREADS" \
               -x "$HISAT2_INDEX" \
               -U "$CLEAN_READ1" \
               -S "$sam_file" 2>> "$LOG_FILE"
    fi

    # Convert SAM to BAM and sort
    log_info "Converting and sorting BAM file"
    if samtools view -bS -@ "$THREADS" "$sam_file" | \
       samtools sort -@ "$THREADS" -o "$sorted_bam" - 2>> "$LOG_FILE"; then
        log_info "BAM conversion and sorting completed successfully"
        
        # Clean up SAM file immediately after successful conversion
        if [ "$KEEP_INTERMEDIATE" = false ]; then
            log_info "Removing SAM file to save disk space"
            rm -f "$sam_file"
        fi
    else
        log_error "Failed to convert SAM to BAM"
        # Clean up SAM file even on failure to save disk space
        rm -f "$sam_file"
        return 1
    fi

    # Create index
    if ! samtools index "$sorted_bam" 2>> "$LOG_FILE"; then
        log_error "Failed to create BAM index"
        return 1
    fi

    # Alignment statistics
    log_info "Generating alignment statistics"
    samtools flagstat "$sorted_bam" > "$OUTPUT_DIR/hisat2/${SAMPLE_NAME}.flagstat" 2>> "$LOG_FILE"

    SORTED_BAM="$sorted_bam"
}

# featureCounts quantification
run_featurecounts() {
    log_info "Using featureCounts for gene quantification"

    if [ "$DATA_TYPE" = "PE" ]; then
        featureCounts \
            -T "$THREADS" \
            -p \
            -t exon \
            -g gene_id \
            -a "$ANNOTATION_FILE" \
            -o "$OUTPUT_DIR/featurecounts/${SAMPLE_NAME}.counts" \
            "$SORTED_BAM" 2>> "$LOG_FILE"
    else
        featureCounts \
            -T "$THREADS" \
            -t exon \
            -g gene_id \
            -a "$ANNOTATION_FILE" \
            -o "$OUTPUT_DIR/featurecounts/${SAMPLE_NAME}.counts" \
            "$SORTED_BAM" 2>> "$LOG_FILE"
    fi

    log_info "featureCounts quantification completed ✓"
}




# Run quantification analysis (fixed to use featureCounts)
run_quantification() {
    log_step "Gene quantification (featureCounts)"

    run_alignment
    run_featurecounts
}

# Generate analysis report
generate_report() {
    log_step "Generating analysis report"

    local report_file="$OUTPUT_DIR/${SAMPLE_NAME}_analysis_report.txt"

    cat > "$report_file" << EOF
=======================================
PanFamily RNAseq Analysis Report
=======================================

Analysis Information:
- Sample Name: $SAMPLE_NAME
- Data Type: $DATA_TYPE
- Quantification Method: featureCounts
- Analysis Time: $(date)

Input Files:
- Read1: $READ1
- Read2: $READ2
- Genome: $GENOME_FILE
- Annotation File: $ANNOTATION_FILE

Data Processing Parameters:
- Quality Control Before Trimming: $QC_BEFORE
- Sequence Trimming: $TRIMMING
- Trimming Quality Threshold: $TRIM_QUALITY
- Minimum Length Threshold: $TRIM_LENGTH
- Quality Control After Trimming: $QC_AFTER

Analysis Parameters:
- Threads: $THREADS
- Keep Intermediate Files: $KEEP_INTERMEDIATE

Result Files:
- Output Directory: $OUTPUT_DIR
- Data Processing Results: $OUTPUT_DIR/data_processing/
- Quantification Results: $OUTPUT_DIR/featurecounts/
- Log File: $LOG_FILE

=======================================
Analysis completed! Please check result files in corresponding directories
=======================================
EOF

    log_info "Analysis report generated: $report_file"
}

# Main program workflow
main() {
    print_banner

    log_info "Starting $PROGRAM_NAME analysis"
    log_info "Sample: $SAMPLE_NAME | Method: featureCounts | Threads: $THREADS | Data Type: $DATA_TYPE"

    create_directories
    check_dependencies
    run_data_processing
    run_quantification
    generate_report

    echo ""
    log_info "🎉 Analysis completed! Results saved in: $OUTPUT_DIR"
    log_info "📊 View report: $OUTPUT_DIR/${SAMPLE_NAME}_analysis_report.txt"
    echo ""
}

# ============================================================================
# Parameter parsing
# ============================================================================

# Default parameter settings
THREADS=8
KEEP_INTERMEDIATE=false
SAMPLE_NAME=""
DATA_TYPE=""
YAML_CONFIG=""

# Default data processing parameters
QC_BEFORE="true"
QC_AFTER="true"
TRIMMING="true"
TRIM_QUALITY="25"
TRIM_LENGTH="35"
FASTQC_OPTIONS=""

# Required parameters
GENOME_FILE=""
ANNOTATION_FILE=""
READ1=""
READ2=""
OUTPUT_DIR=""

# Parse command line parameters
while [[ $# -gt 0 ]]; do
    case $1 in
        -g)
            GENOME_FILE="$2"
            shift 2
            ;;
        -a)
            ANNOTATION_FILE="$2"
            shift 2
            ;;
        -i)
            READ1="$2"
            shift 2
            if [[ $# -gt 0 ]] && [[ ! "$1" =~ ^- ]]; then
                READ2="$1"
                shift
            fi
            ;;
        -o)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -p)
            THREADS="$2"
            shift 2
            ;;
        -s)
            SAMPLE_NAME="$2"
            shift 2
            ;;
        -c)
            YAML_CONFIG="$2"
            shift 2
            ;;
        --qc-before)
            QC_BEFORE="true"
            shift
            ;;
        --qc-after)
            QC_AFTER="true"
            shift
            ;;
        --no-trimming)
            TRIMMING="false"
            shift
            ;;
        --no-qc)
            QC_BEFORE="false"
            QC_AFTER="false"
            shift
            ;;
        --keep-intermediate)
            KEEP_INTERMEDIATE=true
            shift
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

# If YAML configuration file is provided, parse it first
if [[ -n "$YAML_CONFIG" ]]; then
    parse_yaml "$YAML_CONFIG"
fi

# Parameter validation
if [[ -z "$READ1" ]] || [[ -z "$OUTPUT_DIR" ]]; then
    log_error "Missing required parameters: input file(-i) and output directory(-o), or complete YAML configuration file(-c)"
    echo ""
    print_usage
    exit 1
fi

# Determine data type
if [[ -n "$READ2" ]]; then
    DATA_TYPE="PE"
else
    DATA_TYPE="SE"
fi

# Check if files exist
for file in "$READ1" "$READ2" "$GENOME_FILE" "$ANNOTATION_FILE" "$YAML_CONFIG"; do
    if [[ -n "$file" ]] && [[ ! -f "$file" ]]; then
        log_error "File does not exist: $file"
        exit 1
    fi
done

# Auto-extract sample name
if [[ -z "$SAMPLE_NAME" ]]; then
    SAMPLE_NAME=$(basename "$READ1" | sed 's/_R[12].*//g' | sed 's/_1.*//g' | sed 's/\..*//g')
    log_info "Auto-extracted sample name: $SAMPLE_NAME"
fi

# Check parameters required for featureCounts
if [[ -z "$GENOME_FILE" ]]; then
    log_error "Missing required genome FASTA file (-g)"
    print_usage
    exit 1
fi

if [[ -z "$ANNOTATION_FILE" ]]; then
    log_error "Missing required GTF annotation file (-a)"
    print_usage
    exit 1
fi

# Run main program
main