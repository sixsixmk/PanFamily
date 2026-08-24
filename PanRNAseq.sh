#!/bin/bash

set -euo pipefail


PROGRAM_NAME="PanRNAseq"
VERSION="1.0.0"


log_info() {
    echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2
}

log_warn() {
    echo "[WARN] $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2
}

log_error() {
    echo "[ERROR] $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2
}

log_step() {
    echo "[STEP] $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2
}

log_strain() {
    echo "[STRAIN] $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2
}

# Print usage information
print_usage() {
    cat << EOF
USAGE:
    $0 -c <CONFIG_FILE>

OPTIONS:
    -c <FILE>               YAML configuration file
    -h, --help              Show help information
    -v, --version           Show version information

YAML Configuration File Format:
    # Pan-transcriptome batch analysis configuration file
    pan_analysis:
      strain_list: "/path/to/strain_list.txt"      # Strain list file
      transcriptome_dir: "/path/to/transcriptomes" # Transcriptome data directory
      genome_dir: "/path/to/genomes"               # Genome directory
      annotation_dir: "/path/to/annotations"       # Annotation file directory
      output_dir: "/path/to/output"                # Result output directory
      data_type: "PE"                              # PE (paired-end) or SE (single-end)
    
    # Analysis parameters
    analysis:
      threads: 20
      keep_intermediate: false
    
    # Data processing parameters
    data_processing:
      qc_before: true
      qc_after: true
      trimming: true
      trim_quality: 25
      trim_length: 35
      fastqc_options: ""

Strain List File Format (strain_list.txt):
    StrainName\tSampleName
    
    One sample per line, example:
    StrainA    SampleA1
    StrainA    SampleA2
    StrainA    SampleA3
    StrainB    SampleB1
    StrainB    SampleB2

File Naming Rules:
    - Transcriptome data: transcriptome_dir/SampleName_1.fastq.gz, SampleName_2.fastq.gz (PE mode)
    - Genome: genome_dir/StrainName.fa or StrainName.fasta
    - Annotation file: annotation_dir/StrainName.gtf or StrainName.gff3

Output Files:
    - [StrainName]/: Independent analysis directory for each strain
      · linked_data/: Soft-linked transcriptome data
      · [SampleName]/: Analysis results for each sample
      · featureCounts/: Summarized featureCounts results
      · gene_count_matrix.tsv: Gene expression matrix within strain
      · normalized_fpkm.tsv: FPKM normalized matrix
      · normalized_tpm.tsv: TPM normalized matrix
      · strain_analysis_report.txt: Strain analysis report
    - pan_summary/: Pan-transcriptome summary results
      · all_strains_summary.txt: Summary information for all strains
      · strain_comparison.txt: Comparison between strains
    - pan_analysis.log: Pan-transcriptome analysis log

EXAMPLE:
    $0 -c pan_config.yaml

EOF
}
parse_yaml() {
    local yaml_file="$1"
    if [[ ! -f "$yaml_file" ]]; then
        log_error "YAML configuration file does not exist: $yaml_file"
        exit 1
    fi
    
    log_info "Parsing YAML configuration file: $yaml_file"

    # Extract configuration parameters
    STRAIN_LIST=$(grep -E "^\s*strain_list:" "$yaml_file" | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')
    TRANSCRIPTOME_DIR=$(grep -E "^\s*transcriptome_dir:" "$yaml_file" | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')
    GENOME_DIR=$(grep -E "^\s*genome_dir:" "$yaml_file" | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')
    ANNOTATION_DIR=$(grep -E "^\s*annotation_dir:" "$yaml_file" | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')
    OUTPUT_DIR=$(grep -E "^\s*output_dir:" "$yaml_file" | head -1 | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')
    DATA_TYPE=$(grep -E "^\s*data_type:" "$yaml_file" | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')
    THREADS=$(grep -E "^\s*threads:" "$yaml_file" | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')
    KEEP_INTERMEDIATE=$(grep -E "^\s*keep_intermediate:" "$yaml_file" | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')
    QC_BEFORE=$(grep -E "^\s*qc_before:" "$yaml_file" | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')
    QC_AFTER=$(grep -E "^\s*qc_after:" "$yaml_file" | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')
    TRIMMING=$(grep -E "^\s*trimming:" "$yaml_file" | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')
    TRIM_QUALITY=$(grep -E "^\s*trim_quality:" "$yaml_file" | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')
    TRIM_LENGTH=$(grep -E "^\s*trim_length:" "$yaml_file" | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')
    FASTQC_OPTIONS=$(grep -E "^\s*fastqc_options:" "$yaml_file" | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')
    
    # Sample grouping parameters (shared by DEG and SV_pre)
    SAMPLE_GROUP_FILE=$(grep -A 3 "sample_grouping:" "$yaml_file" | grep -E "^\s*group_file:" | head -1 | sed 's/#.*//' | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*$//' || true)
    
    # DEG analysis parameters
    DEG_ENABLED=$(grep -A 10 "deg_analysis:" "$yaml_file" | grep -E "^\s*enabled:" | head -1 | sed 's/#.*//' | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*$//' || true)
    DEG_MATRIX_TYPE=$(grep -A 10 "deg_analysis:" "$yaml_file" | grep -E "^\s*matrix_type:" | head -1 | sed 's/#.*//' | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*$//' || true)
    DEG_PVALUE_CUTOFF=$(grep -A 10 "deg_analysis:" "$yaml_file" | grep -E "^\s*pvalue_cutoff:" | head -1 | sed 's/#.*//' | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*$//' || true)
    DEG_LOG2FC_CUTOFF=$(grep -A 10 "deg_analysis:" "$yaml_file" | grep -E "^\s*log2fc_cutoff:" | head -1 | sed 's/#.*//' | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*$//' || true)
    
    # Gene family parameters
    GENEFAMILY_LIST=$(grep -E "^\s*genefamily_list:" "$yaml_file" | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')
    GENEFAMILY_PREFIX=$(grep -E "^\s*genefamily_prefix:" "$yaml_file" | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')
    GENEFAMILY_EXTRACT_DEG=$(grep -E "^\s*extract_from_deg:" "$yaml_file" | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*#.*$//' | sed 's/[[:space:]]*$//')
    
    # Gene cluster parameters (original ID mapping)
    GENE_CLUSTER_ENABLED=$(grep -A 15 "gene_cluster:" "$yaml_file" | grep -E "^\s*enabled:" | head -1 | sed 's/#.*//' | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*$//' || true)
    GENE_CLUSTER_FILE=$(grep -A 15 "gene_cluster:" "$yaml_file" | grep -E "^\s*cluster_file:" | head -1 | sed 's/#.*//' | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*$//' || true)
    GENE_CLUSTER_OUTPUT_DIR=$(grep -A 15 "gene_cluster:" "$yaml_file" | grep -E "^\s*output_dir:" | head -1 | sed 's/#.*//' | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*$//' || true)
    
    # SV_pre parameters
    SV_PRE_ENABLED=$(grep -A 5 "sv_pre:" "$yaml_file" | grep -E "^\s*enabled:" | head -1 | sed 's/#.*//' | sed 's/.*:[[:space:]]*//' | tr -d '"' | sed 's/[[:space:]]*$//' || true)
    
    # Set default values
    STRAIN_LIST=${STRAIN_LIST:-""}
    TRANSCRIPTOME_DIR=${TRANSCRIPTOME_DIR:-""}
    GENOME_DIR=${GENOME_DIR:-""}
    ANNOTATION_DIR=${ANNOTATION_DIR:-""}
    OUTPUT_DIR=${OUTPUT_DIR:-"./pan_results"}
    DATA_TYPE=${DATA_TYPE:-PE}
    THREADS=${THREADS:-8}
    KEEP_INTERMEDIATE=${KEEP_INTERMEDIATE:-false}
    QC_BEFORE=${QC_BEFORE:-true}
    QC_AFTER=${QC_AFTER:-true}
    TRIMMING=${TRIMMING:-true}
    TRIM_QUALITY=${TRIM_QUALITY:-25}
    TRIM_LENGTH=${TRIM_LENGTH:-35}
    FASTQC_OPTIONS=${FASTQC_OPTIONS:-""}
    
    # Sample grouping default values
    SAMPLE_GROUP_FILE=${SAMPLE_GROUP_FILE:-""}
    
    # DEG analysis default values
    DEG_ENABLED=${DEG_ENABLED:-false}
    DEG_MATRIX_TYPE=${DEG_MATRIX_TYPE:-"count"}
    DEG_PVALUE_CUTOFF=${DEG_PVALUE_CUTOFF:-0.05}
    DEG_LOG2FC_CUTOFF=${DEG_LOG2FC_CUTOFF:-1}
    
    # Gene family default values
    GENEFAMILY_LIST=${GENEFAMILY_LIST:-""}
    GENEFAMILY_PREFIX=${GENEFAMILY_PREFIX:-""}
    GENEFAMILY_EXTRACT_DEG=${GENEFAMILY_EXTRACT_DEG:-false}
    
    # Gene cluster default values
    GENE_CLUSTER_ENABLED=${GENE_CLUSTER_ENABLED:-false}
    GENE_CLUSTER_FILE=${GENE_CLUSTER_FILE:-""}
    GENE_CLUSTER_OUTPUT_DIR=${GENE_CLUSTER_OUTPUT_DIR:-"gene_cluster_results"}
    
    # SV_pre default values
    SV_PRE_ENABLED=${SV_PRE_ENABLED:-false}

    log_info "YAML configuration parsing completed"
}

# Validate configuration parameters
validate_config() {
    log_step "Validating configuration parameters"
    
    local has_error=false
    
    # Validate strain list file
    if [[ -z "$STRAIN_LIST" ]] || [[ ! -f "$STRAIN_LIST" ]]; then
        log_error "Strain list file does not exist or not specified: $STRAIN_LIST"
        has_error=true
    fi
    
    # Validate transcriptome directory
    if [[ -z "$TRANSCRIPTOME_DIR" ]] || [[ ! -d "$TRANSCRIPTOME_DIR" ]]; then
        log_error "Transcriptome data directory does not exist or not specified: $TRANSCRIPTOME_DIR"
        has_error=true
    fi
    
    # Validate genome directory
    if [[ -z "$GENOME_DIR" ]] || [[ ! -d "$GENOME_DIR" ]]; then
        log_error "Genome directory does not exist or not specified: $GENOME_DIR"
        has_error=true
    fi
    
    # Validate annotation file directory
    if [[ -z "$ANNOTATION_DIR" ]] || [[ ! -d "$ANNOTATION_DIR" ]]; then
        log_error "Annotation file directory does not exist or not specified: $ANNOTATION_DIR"
        has_error=true
    fi
    
    if [[ "$has_error" == true ]]; then
        log_error "Configuration validation failed, please check configuration file"
        exit 1
    fi
    
    log_info "Configuration parameters validated successfully"
}

# Parse strain list
parse_strain_list() {
    log_step "Parsing strain list"
    
    declare -gA STRAIN_SAMPLES  # Associative array: strain name -> sample list (space-separated)
    declare -ga STRAINS         # Strain name array (maintain order, deduplicated)
    declare -A strain_seen      # For deduplication
    
    local line_num=0
    while IFS=$'\t' read -r strain sample || [[ -n "$strain" ]]; do
        line_num=$((line_num + 1))
        
        # Skip empty lines and comment lines
        [[ -z "$strain" ]] && continue
        [[ "$strain" =~ ^#.* ]] && continue
        
        # Remove leading and trailing whitespace
        strain=$(echo "$strain" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        sample=$(echo "$sample" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        
        if [[ -z "$strain" ]] || [[ -z "$sample" ]]; then
            log_warn "Line $line_num format incorrect, skipping"
            continue
        fi
        
        # If new strain, add to STRAINS array
        if [[ -z "${strain_seen[$strain]:-}" ]]; then
            STRAINS+=("$strain")
            strain_seen["$strain"]=1
            STRAIN_SAMPLES["$strain"]=""
        fi
        
        # Add sample to the strain's sample list (space-separated)
        if [[ -z "${STRAIN_SAMPLES[$strain]:-}" ]]; then
            STRAIN_SAMPLES["$strain"]="$sample"
        else
            STRAIN_SAMPLES["$strain"]="${STRAIN_SAMPLES[$strain]} $sample"
        fi
        
    done < "$STRAIN_LIST"
    
    if [[ ${#STRAINS[@]} -eq 0 ]]; then
        log_error "No strains parsed from strain list"
        exit 1
    fi
    
    # Output parsing results
    log_info "Total ${#STRAINS[@]} strains parsed"
    for strain in "${STRAINS[@]}"; do
        local sample_count=$(echo "${STRAIN_SAMPLES[$strain]}" | wc -w)
        log_info "  Strain: $strain, Sample count: $sample_count"
    done
}

# Find genome file for strain
find_genome_file() {
    local strain="$1"
    
    for ext in "fa" "fasta" "fna" "fa.gz" "fasta.gz" "fna.gz"; do
        local genome_file="$GENOME_DIR/${strain}.${ext}"
        if [[ -f "$genome_file" ]]; then
            echo "$genome_file"
            return 0
        fi
    done
    
    # Return empty when not found, no error
    return 1
}

# Find annotation file for strain
find_annotation_file() {
    local strain="$1"
    
    for ext in "gtf" "gff" "gff3" "gtf.gz" "gff.gz" "gff3.gz"; do
        local annot_file="$ANNOTATION_DIR/${strain}.${ext}"
        if [[ -f "$annot_file" ]]; then
            echo "$annot_file"
            return 0
        fi
    done
    
    # Return empty when not found, no error
    return 1
}

# Create soft links for strain (transcriptome data, genome, annotation files)
create_strain_links() {
    local strain="$1"
    local samples="$2"
    local strain_dir="$OUTPUT_DIR/$strain"
    
    log_info "Creating data soft links for strain $strain"
    
    mkdir -p "$strain_dir"
    
    local linked_count=0
    local failed_samples=()
    
    # 1. Link transcriptome data
    # Convert sample list to array (space-separated)
    read -ra sample_array <<< "$samples"
    
    for sample in "${sample_array[@]}"; do
        # Remove whitespace
        sample=$(echo "$sample" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        
        if [[ "$DATA_TYPE" == "PE" ]]; then
            # Paired-end data
            local found=false
            for ext in "fastq.gz" "fq.gz" "fastq" "fq"; do
                local read1="$TRANSCRIPTOME_DIR/${sample}_1.${ext}"
                local read2="$TRANSCRIPTOME_DIR/${sample}_2.${ext}"
                
                if [[ -f "$read1" ]] && [[ -f "$read2" ]]; then
                    ln -sf "$(realpath "$read1")" "$strain_dir/${sample}_1.${ext}"
                    ln -sf "$(realpath "$read2")" "$strain_dir/${sample}_2.${ext}"
                    log_info "  Linked transcriptome sample: $sample (PE)"
                    ((linked_count++))
                    found=true
                    break
                fi
            done
            
            if [[ "$found" == false ]]; then
                log_warn "  Paired files not found for sample $sample, skipping"
                failed_samples+=("$sample")
            fi
        else
            # Single-end data
            local found=false
            for ext in "fastq.gz" "fq.gz" "fastq" "fq"; do
                local read1="$TRANSCRIPTOME_DIR/${sample}.${ext}"
                
                if [[ -f "$read1" ]]; then
                    ln -sf "$(realpath "$read1")" "$strain_dir/${sample}.${ext}"
                    log_info "  Linked transcriptome sample: $sample (SE)"
                    ((linked_count++))
                    found=true
                    break
                fi
            done
            
            if [[ "$found" == false ]]; then
                log_warn "  Data file not found for sample $sample, skipping"
                failed_samples+=("$sample")
            fi
        fi
    done
    
    # 2. Link genome file
    local genome_file=$(find_genome_file "$strain")
    if [[ -n "$genome_file" ]]; then
        local genome_basename=$(basename "$genome_file")
        ln -sf "$(realpath "$genome_file")" "$strain_dir/$genome_basename"
        log_info "  Linked genome: $genome_basename"
    else
        log_warn "  Genome file not found for strain $strain"
    fi
    
    # 3. Link annotation file
    local annotation_file=$(find_annotation_file "$strain")
    if [[ -n "$annotation_file" ]]; then
        local annotation_basename=$(basename "$annotation_file")
        ln -sf "$(realpath "$annotation_file")" "$strain_dir/$annotation_basename"
        log_info "  Linked annotation file: $annotation_basename"
    else
        log_warn "  Annotation file not found for strain $strain"
    fi
    
    log_info "Strain $strain: Successfully linked $linked_count transcriptome samples"
    
    if [[ ${#failed_samples[@]} -gt 0 ]]; then
        log_warn "The following samples have no data files (skipped):"
        for sample in "${failed_samples[@]}"; do
            log_warn "  - $sample"
        done
    fi
    
    echo "$linked_count"
}

# Extract samples for specified strain from 3-column DEG grouping file and generate 2-column format file
# Parameters: $1 - strain name, $2 - 3-column grouping file path, $3 - output 2-column format file path
# Return: 0 - success, 1 - failure
generate_strain_deg_group() {
    local strain="$1"
    local input_file="$2"
    local output_file="$3"
    
    # Check if input file exists
    if [[ ! -f "$input_file" ]]; then
        return 1
    fi
    
    # Create output file directory
    local output_dir=$(dirname "$output_file")
    mkdir -p "$output_dir"
    
    # Extract samples for this strain
    local control_samples=()
    local treatment_samples=()
    local line_num=0
    
    while IFS=$'\t' read -r strain_col sample_col group_col || [[ -n "$strain_col" ]]; do
        line_num=$((line_num + 1))
        
        # Skip header
        if [[ $line_num -eq 1 ]]; then
            if [[ "$strain_col" == "Strain" ]] || [[ "$strain_col" =~ ^[[:space:]]*Strain ]]; then
                continue
            fi
        fi
        
        # Skip empty lines
        [[ -z "$strain_col" ]] && continue
        [[ -z "$sample_col" ]] && continue
        
        # Remove leading and trailing whitespace
        strain_col=$(echo "$strain_col" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        sample_col=$(echo "$sample_col" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        group_col=$(echo "$group_col" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        
        # Only process samples for specified strain
        if [[ "$strain_col" == "$strain" ]]; then
            if [[ "$group_col" == "Control" ]]; then
                control_samples+=("$sample_col")
            elif [[ "$group_col" == "Treatment" ]]; then
                treatment_samples+=("$sample_col")
            fi
        fi
    done < "$input_file"
    
    # Check if samples were found
    if [[ ${#control_samples[@]} -eq 0 ]] && [[ ${#treatment_samples[@]} -eq 0 ]]; then
        return 1
    fi
    
    # Generate 2-column format grouping file
    # Write header
    echo -e "Control\tTreatment" > "$output_file"
    
    # Find maximum number of lines
    local max_lines=${#control_samples[@]}
    if [[ ${#treatment_samples[@]} -gt $max_lines ]]; then
        max_lines=${#treatment_samples[@]}
    fi
    
    # Write samples line by line
    for ((i=0; i<max_lines; i++)); do
        local ctrl_sample=""
        local treat_sample=""
        
        if [[ $i -lt ${#control_samples[@]} ]]; then
            ctrl_sample="${control_samples[$i]}"
        fi
        
        if [[ $i -lt ${#treatment_samples[@]} ]]; then
            treat_sample="${treatment_samples[$i]}"
        fi
        
        echo -e "${ctrl_sample}\t${treat_sample}" >> "$output_file"
    done
    
    return 0
}

# Generate configuration file for single strain
generate_strain_config() {
    local strain="$1"
    local strain_dir="$OUTPUT_DIR/$strain"
    local config_file="$strain_dir/strain_config.yaml"
    
    # Find genome and annotation files (should be in strain directory after linking)
    local genome_file=""
    local annotation_file=""
    
    # First try to find linked files in strain directory
    for ext in "fa" "fasta" "fna" "fa.gz" "fasta.gz" "fna.gz"; do
        if [[ -f "$strain_dir/${strain}.${ext}" ]]; then
            genome_file="$strain_dir/${strain}.${ext}"
            break
        fi
    done
    
    for ext in "gtf" "gff" "gff3" "gtf.gz" "gff.gz" "gff3.gz"; do
        if [[ -f "$strain_dir/${strain}.${ext}" ]]; then
            annotation_file="$strain_dir/${strain}.${ext}"
            break
        fi
    done
    
    # If not in strain directory, try original directory
    if [[ -z "$genome_file" ]]; then
        genome_file=$(find_genome_file "$strain")
    fi
    
    if [[ -z "$annotation_file" ]]; then
        annotation_file=$(find_annotation_file "$strain")
    fi
    
    # Check if required files were found
    if [[ -z "$genome_file" ]]; then
        log_error "Genome file not found for strain $strain"
        return 1
    fi
    
    if [[ -z "$annotation_file" ]]; then
        log_error "Annotation file not found for strain $strain"
        return 1
    fi
    
    log_info "Generating configuration file for strain $strain"
    log_info "  Genome: $genome_file"
    log_info "  Annotation: $annotation_file"
    
    # Process DEG grouping file (generate 2-column format from 3-column format)
    local strain_deg_group_file=""
    if [[ -n "$SAMPLE_GROUP_FILE" ]] && [[ -f "$SAMPLE_GROUP_FILE" ]]; then
        # Extract samples for this strain from 3-column format and generate 2-column format
        local temp_deg_file="$strain_dir/.deg_group_${strain}.txt"
        
        if generate_strain_deg_group "$strain" "$SAMPLE_GROUP_FILE" "$temp_deg_file"; then
            strain_deg_group_file="$temp_deg_file"
            
            if [[ "$DEG_ENABLED" == "true" ]] || [[ "$DEG_ENABLED" == "TRUE" ]]; then
                log_info "  DEG grouping file: $strain_deg_group_file"
            fi
        else
            if [[ "$DEG_ENABLED" == "true" ]] || [[ "$DEG_ENABLED" == "TRUE" ]]; then
                log_warn "No samples found for strain $strain in grouping file, DEG analysis for this strain will be skipped"
            fi
        fi
    elif [[ "$DEG_ENABLED" == "true" ]] || [[ "$DEG_ENABLED" == "TRUE" ]]; then
        if [[ -n "$SAMPLE_GROUP_FILE" ]]; then
            log_warn "Sample grouping file does not exist: $SAMPLE_GROUP_FILE"
        fi
    fi
    
    # Generate configuration file
    cat > "$config_file" << EOF
# Transcriptome analysis configuration file for strain $strain
# Auto-generated at: $(date)

batch_analysis:
  data_dir: "$strain_dir"
  output_dir: "$strain_dir"
  data_type: "$DATA_TYPE"

reference:
  genome_file: "$genome_file"
  annotation_file: "$annotation_file"

analysis:
  threads: $THREADS
  keep_intermediate: $KEEP_INTERMEDIATE

data_processing:
  qc_before: $QC_BEFORE
  qc_after: $QC_AFTER
  trimming: $TRIMMING
  trim_quality: $TRIM_QUALITY
  trim_length: $TRIM_LENGTH
  fastqc_options: "$FASTQC_OPTIONS"

# Differential expression analysis (optional)
deg_analysis:
  enabled: $DEG_ENABLED
  group_file: "$strain_deg_group_file"
  matrix_type: "$DEG_MATRIX_TYPE"
  pvalue_cutoff: $DEG_PVALUE_CUTOFF
  log2fc_cutoff: $DEG_LOG2FC_CUTOFF

# Gene family extraction (optional)
gene_family:
  genefamily_list: "$GENEFAMILY_LIST"
  genefamily_prefix: "$GENEFAMILY_PREFIX"
  extract_from_deg: $GENEFAMILY_EXTRACT_DEG
EOF
    
    echo "$config_file"
}

# Process single strain
process_strain() {
    local strain="$1"
    local strain_num="$2"
    local total="$3"
    
    log_strain "Start processing strain [$strain_num/$total]: $strain"
    
    local strain_dir="$OUTPUT_DIR/$strain"
    mkdir -p "$strain_dir"
    
    # Get sample list
    local samples="${STRAIN_SAMPLES[$strain]}"
    
    # Create soft links
    local linked_count=$(create_strain_links "$strain" "$samples")
    
    if [[ $linked_count -eq 0 ]]; then
        log_error "Strain $strain has no valid sample data, skipping"
        return 1
    fi
    
    # Generate configuration file
    local config_file=$(generate_strain_config "$strain")
    
    if [[ -z "$config_file" ]]; then
        log_error "Configuration file generation failed for strain $strain, skipping"
        return 1
    fi
    
    # Call RNAseq.sh for analysis
    log_info "Calling RNAseq.sh to analyze strain: $strain"
    
    if bash "$RNASEQ_SCRIPT" -c "$config_file" >&2; then
        log_info "Strain $strain analysis completed"
        
        # Generate strain report
        generate_strain_report "$strain" "$strain_dir"
        
        return 0
    else
        log_error "Strain $strain analysis failed"
        return 1
    fi
}

# Collect DEG results from all strains to unified directory
collect_deg_results() {
    log_step "Collecting DEG results from all strains"
    
    local deg_summary_dir="$OUTPUT_DIR/DEG_results"
    mkdir -p "$deg_summary_dir"
    
    local collected_count=0
    
    for strain in "${STRAINS[@]}"; do
        local strain_dir="$OUTPUT_DIR/$strain"
        local strain_deg_dir="$strain_dir/DEG_results"
        
        # Check if strain's DEG results directory exists
        if [[ ! -d "$strain_deg_dir" ]]; then
            continue
        fi
        
        # Check if there are DEG files
        local deg_files=$(find "$strain_deg_dir" -name "DEG_*.csv" 2>/dev/null | wc -l)
        if [[ $deg_files -eq 0 ]]; then
            continue
        fi
        
        # Create strain subdirectory
        local strain_deg_summary="$deg_summary_dir/$strain"
        mkdir -p "$strain_deg_summary"
        
        # Copy DEG result files
        cp "$strain_deg_dir"/DEG_*.csv "$strain_deg_summary/" 2>/dev/null || true
        
        local copied_files=$(ls "$strain_deg_summary"/DEG_*.csv 2>/dev/null | wc -l)
        if [[ $copied_files -gt 0 ]]; then
            log_info "Strain $strain: Collected $copied_files DEG files"
            collected_count=$((collected_count + 1))
        fi
    done
    
    if [[ $collected_count -gt 0 ]]; then
        log_info "✓ Successfully collected DEG results from $collected_count strains"
        log_info "  Saved to: $deg_summary_dir"
    else
        log_warn "No DEG results found from any strain"
    fi
    
    return 0
}

# Gene clustering and matrix integration
cluster_and_integrate_matrices() {
    log_step "Gene clustering and matrix integration"
    
    # Check if enabled
    if [[ "$GENE_CLUSTER_ENABLED" != "true" ]] && [[ "$GENE_CLUSTER_ENABLED" != "TRUE" ]]; then
        log_info "Gene clustering function not enabled, skipping"
        return 0
    fi
    
    # Check gene cluster file
    if [[ -z "$GENE_CLUSTER_FILE" ]] || [[ ! -f "$GENE_CLUSTER_FILE" ]]; then
        log_warn "Gene cluster file does not exist or not specified: $GENE_CLUSTER_FILE"
        log_warn "Skipping gene clustering and matrix integration"
        return 0
    fi
    
    log_info "Gene cluster file: $GENE_CLUSTER_FILE"
    
    # Create output directory
    local cluster_dir="$OUTPUT_DIR/$GENE_CLUSTER_OUTPUT_DIR"
    mkdir -p "$cluster_dir"/{count,fpkm,tpm,deg}
    
    log_info "Output directory: $cluster_dir"
    
    # Check Python and rename script
    if ! command -v python3 &> /dev/null; then
        log_error "python3 not found, cannot execute gene clustering"
        return 1
    fi
    
    local script_dir=$(cd "$(dirname "$0")" && pwd)
    local rename_script="$script_dir/rename_gene_ids.py"
    
    if [[ ! -f "$rename_script" ]]; then
        log_error "Clustering script does not exist: $rename_script"
        return 1
    fi
    
    # Rename matrix files for each strain
    local clustered_count=0
    
    for strain in "${STRAINS[@]}"; do
        local strain_dir="$OUTPUT_DIR/$strain"
        
        log_info "Processing strain: $strain"
        
        # Rename count matrix
        if [[ -f "$strain_dir/${strain}_count_matrix.tsv" ]]; then
            python3 "$rename_script" \
                -i "$strain_dir/${strain}_count_matrix.tsv" \
                -o "$cluster_dir/count/${strain}_count_matrix.tsv" \
                -m "$GENE_CLUSTER_FILE" \
                -s "$strain" \
                -t matrix 2>&1 | while IFS= read -r line; do
                    log_info "  $line"
                done
        fi
        
        # Rename FPKM matrix
        if [[ -f "$strain_dir/${strain}_fpkm_matrix.tsv" ]]; then
            python3 "$rename_script" \
                -i "$strain_dir/${strain}_fpkm_matrix.tsv" \
                -o "$cluster_dir/fpkm/${strain}_fpkm_matrix.tsv" \
                -m "$GENE_CLUSTER_FILE" \
                -s "$strain" \
                -t matrix 2>&1 | while IFS= read -r line; do
                    log_info "  $line"
                done
        fi
        
        # Rename TPM matrix
        if [[ -f "$strain_dir/${strain}_tpm_matrix.tsv" ]]; then
            python3 "$rename_script" \
                -i "$strain_dir/${strain}_tpm_matrix.tsv" \
                -o "$cluster_dir/tpm/${strain}_tpm_matrix.tsv" \
                -m "$GENE_CLUSTER_FILE" \
                -s "$strain" \
                -t matrix 2>&1 | while IFS= read -r line; do
                    log_info "  $line"
                done
        fi
        
        # Rename DEG results
        local strain_deg_dir="$strain_dir/DEG_results"
        if [[ -d "$strain_deg_dir" ]]; then
            mkdir -p "$cluster_dir/deg/$strain"
            
            for deg_file in "$strain_deg_dir"/DEG_*.csv; do
                [[ ! -f "$deg_file" ]] && continue
                
                local deg_basename=$(basename "$deg_file")
                python3 "$rename_script" \
                    -i "$deg_file" \
                    -o "$cluster_dir/deg/$strain/$deg_basename" \
                    -m "$GENE_CLUSTER_FILE" \
                    -t deg 2>&1 | while IFS= read -r line; do
                        log_info "  $line"
                    done
            done
        fi
        
        clustered_count=$((clustered_count + 1))
    done
    
    log_info "✓ Completed gene clustering for $clustered_count strains"
    
    # Integrate clustered matrices
    local merge_script="$script_dir/merge_renamed_matrices.py"
    
    if [[ ! -f "$merge_script" ]]; then
        log_warn "Matrix integration script does not exist: $merge_script"
        log_warn "Skipping matrix integration step"
        return 0
    fi
    
    log_info "Integrating clustered matrices..."
    
    # Integrate count matrix
    if [[ $(ls "$cluster_dir/count"/*.tsv 2>/dev/null | wc -l) -gt 0 ]]; then
        log_info "Integrating count matrix..."
        python3 "$merge_script" \
            "$cluster_dir/count" \
            "$cluster_dir/pan_count_matrix.tsv" \
            "count" 2>&1 | while IFS= read -r line; do
                log_info "  $line"
            done
    fi
    
    # Integrate FPKM matrix
    if [[ $(ls "$cluster_dir/fpkm"/*.tsv 2>/dev/null | wc -l) -gt 0 ]]; then
        log_info "Integrating FPKM matrix..."
        python3 "$merge_script" \
            "$cluster_dir/fpkm" \
            "$cluster_dir/pan_fpkm_matrix.tsv" \
            "fpkm" 2>&1 | while IFS= read -r line; do
                log_info "  $line"
            done
    fi
    
    # Integrate TPM matrix
    if [[ $(ls "$cluster_dir/tpm"/*.tsv 2>/dev/null | wc -l) -gt 0 ]]; then
        log_info "Integrating TPM matrix..."
        python3 "$merge_script" \
            "$cluster_dir/tpm" \
            "$cluster_dir/pan_tpm_matrix.tsv" \
            "tpm" 2>&1 | while IFS= read -r line; do
                log_info "  $line"
            done
    fi
    
    log_info "✓ Gene clustering and matrix integration completed"
    
    # Extract gene family from renamed matrices
    if [[ -n "$GENEFAMILY_PREFIX" ]]; then
        log_info "Extracting gene family ($GENEFAMILY_PREFIX) from renamed matrices..."
        
        local genefamily_dir="$cluster_dir/${GENEFAMILY_PREFIX}_genefamily"
        mkdir -p "$genefamily_dir"/{count,fpkm,tpm}
        
        # Extract from each strain's renamed matrix
        for matrix_type in count fpkm tpm; do
            for strain in "${STRAINS[@]}"; do
                local input_file="$cluster_dir/$matrix_type/${strain}_${matrix_type}_matrix.tsv"
                local output_file="$genefamily_dir/$matrix_type/${strain}_${GENEFAMILY_PREFIX}_${matrix_type}.tsv"
                
                if [[ -f "$input_file" ]]; then
                    (head -1 "$input_file"; grep "^${GENEFAMILY_PREFIX}" "$input_file") > "$output_file"
                    local gene_count=$(tail -n +2 "$output_file" | wc -l)
                    log_info "  Extracted $gene_count ${GENEFAMILY_PREFIX} genes from $strain $matrix_type matrix"
                fi
            done
        done
        
        # Extract from pan matrices
        for matrix_type in count fpkm tpm; do
            local input_file="$cluster_dir/pan_${matrix_type}_matrix.tsv"
            local output_file="$genefamily_dir/pan_${GENEFAMILY_PREFIX}_${matrix_type}.tsv"
            
            if [[ -f "$input_file" ]]; then
                (head -1 "$input_file"; grep "^${GENEFAMILY_PREFIX}" "$input_file") > "$output_file"
                local gene_count=$(tail -n +2 "$output_file" | wc -l)
                log_info "  Extracted $gene_count ${GENEFAMILY_PREFIX} genes from pan_${matrix_type}_matrix"
            fi
        done
        
        log_info "✓ Gene family extraction completed: $genefamily_dir"
    fi
    
    return 0
}

# SV pre-analysis: Extract Control group average expression
generate_sv_pre_matrices() {
    log_step "Generating SV_pre analysis matrices"
    
    # Check if gene_cluster is enabled
    if [[ "$GENE_CLUSTER_ENABLED" != "true" ]] && [[ "$GENE_CLUSTER_ENABLED" != "TRUE" ]]; then
        log_info "Gene clustering not enabled, skipping SV_pre"
        return 0
    fi
    
    # Check if sv_pre is enabled
    if [[ "$SV_PRE_ENABLED" != "true" ]] && [[ "$SV_PRE_ENABLED" != "TRUE" ]]; then
        log_info "SV_pre function not enabled, skipping"
        return 0
    fi
    
    # Check sample grouping file
    if [[ -z "$SAMPLE_GROUP_FILE" ]] || [[ ! -f "$SAMPLE_GROUP_FILE" ]]; then
        log_warn "Sample grouping file does not exist, cannot identify Control group, skipping SV_pre"
        return 0
    fi
    
    # Check gene_cluster output directory
    local cluster_dir="$OUTPUT_DIR/$GENE_CLUSTER_OUTPUT_DIR"
    if [[ ! -d "$cluster_dir/count" ]]; then
        log_warn "Gene clustering results do not exist, skipping SV_pre"
        return 0
    fi
    
    log_info "Sample grouping file: $SAMPLE_GROUP_FILE"
    
    # Create sv_pre output directory
    local sv_pre_dir="$OUTPUT_DIR/sv_pre"
    mkdir -p "$sv_pre_dir"/{count,fpkm,tpm}
    
    log_info "SV_pre output directory: $sv_pre_dir"
    
    # Check Python and calculation script
    if ! command -v python3 &> /dev/null; then
        log_error "python3 not found, cannot execute SV_pre analysis"
        return 1
    fi
    
    local script_dir=$(cd "$(dirname "$0")" && pwd)
    local calc_script="$script_dir/calculate_control_average.py"
    
    if [[ ! -f "$calc_script" ]]; then
        log_error "Control average calculation script does not exist: $calc_script"
        return 1
    fi
    
    # Determine source directory: use genefamily directory if exists, otherwise use full matrix
    local source_dir="$cluster_dir"
    if [[ -n "$GENEFAMILY_PREFIX" ]] && [[ -d "$cluster_dir/${GENEFAMILY_PREFIX}_genefamily" ]]; then
        source_dir="$cluster_dir/${GENEFAMILY_PREFIX}_genefamily"
        log_info "Using gene family matrices from: $source_dir"
    else
        log_info "Using full renamed matrices from: $source_dir"
    fi
    
    # Process each strain
    local processed_count=0
    
    for strain in "${STRAINS[@]}"; do
        log_info "Processing strain: $strain"
        
        # Extract Control sample names for this strain from grouping file
        local control_samples=()
        local line_num=0
        
        while IFS=$'\t' read -r strain_col sample_col group_col || [[ -n "$strain_col" ]]; do
            line_num=$((line_num + 1))
            
            # Skip header
            if [[ $line_num -eq 1 ]]; then
                if [[ "$strain_col" == "Strain" ]] || [[ "$strain_col" =~ ^[[:space:]]*Strain ]]; then
                    continue
                fi
            fi
            
            # Skip empty lines
            [[ -z "$strain_col" ]] && continue
            [[ -z "$sample_col" ]] && continue
            
            # Remove whitespace
            strain_col=$(echo "$strain_col" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
            sample_col=$(echo "$sample_col" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
            group_col=$(echo "$group_col" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
            
            # Only extract Control group samples for this strain
            if [[ "$strain_col" == "$strain" ]] && [[ "$group_col" == "Control" ]]; then
                control_samples+=("$sample_col")
            fi
        done < "$SAMPLE_GROUP_FILE"
        
        if [[ ${#control_samples[@]} -eq 0 ]]; then
            log_warn "Strain $strain has no Control group samples, skipping"
            continue
        fi
        
        log_info "  Found ${#control_samples[@]} Control samples"
        
        # Create temporary sample list file
        local control_list_file="$sv_pre_dir/.${strain}_control_samples.txt"
        printf "%s\n" "${control_samples[@]}" > "$control_list_file"
        
        # Process three matrix types separately
        for mtype in count fpkm tpm; do
            # Try to find matrix file with different naming patterns
            local input_matrix=""
            if [[ -n "$GENEFAMILY_PREFIX" ]] && [[ "$source_dir" == *"${GENEFAMILY_PREFIX}_genefamily"* ]]; then
                # Gene family matrix naming: Blh-1_AtTPS_count.tsv
                input_matrix="$source_dir/${mtype}/${strain}_${GENEFAMILY_PREFIX}_${mtype}.tsv"
            else
                # Full matrix naming: Blh-1_count_matrix.tsv
                input_matrix="$source_dir/${mtype}/${strain}_${mtype}_matrix.tsv"
            fi
            
            local output_file="$sv_pre_dir/${mtype}/${strain}_control_${mtype}.tsv"
            
            if [[ ! -f "$input_matrix" ]]; then
                log_warn "  Matrix file does not exist: $input_matrix"
                continue
            fi
            
            # Call Python script to calculate average
            log_info "  Calculating ${mtype} average..."
            python3 "$calc_script" \
                -i "$input_matrix" \
                -o "$output_file" \
                -c "$control_list_file" \
                -s "$strain" 2>&1 | while IFS= read -r line; do
                    log_info "    $line"
                done
        done
        
        # Clean up temporary file
        rm -f "$control_list_file"
        
        processed_count=$((processed_count + 1))
    done
    
    log_info "✓ Completed SV_pre processing for $processed_count strains"
    
    # Integrate into pan-genome matrices
    local merge_script="$script_dir/merge_renamed_matrices.py"
    
    if [[ ! -f "$merge_script" ]]; then
        log_warn "Matrix integration script does not exist: $merge_script"
        log_warn "Skipping SV_pre matrix integration"
        return 0
    fi
    
    log_info "Integrating SV_pre matrices..."
    
    # Integrate count matrix
    if [[ $(ls "$sv_pre_dir/count"/*.tsv 2>/dev/null | wc -l) -gt 0 ]]; then
        log_info "Integrating control count matrix..."
        python3 "$merge_script" \
            "$sv_pre_dir/count" \
            "$sv_pre_dir/pan_control_count.tsv" \
            "count" 2>&1 | while IFS= read -r line; do
                log_info "  $line"
            done
    fi
    
    # Integrate FPKM matrix
    if [[ $(ls "$sv_pre_dir/fpkm"/*.tsv 2>/dev/null | wc -l) -gt 0 ]]; then
        log_info "Integrating control FPKM matrix..."
        python3 "$merge_script" \
            "$sv_pre_dir/fpkm" \
            "$sv_pre_dir/pan_control_fpkm.tsv" \
            "fpkm" 2>&1 | while IFS= read -r line; do
                log_info "  $line"
            done
    fi
    
    # Integrate TPM matrix
    if [[ $(ls "$sv_pre_dir/tpm"/*.tsv 2>/dev/null | wc -l) -gt 0 ]]; then
        log_info "Integrating control TPM matrix..."
        python3 "$merge_script" \
            "$sv_pre_dir/tpm" \
            "$sv_pre_dir/pan_control_tpm.tsv" \
            "tpm" 2>&1 | while IFS= read -r line; do
                log_info "  $line"
            done
    fi
    
    log_info "✓ SV_pre analysis completed"
    
    return 0
}

# Generate strain report
generate_strain_report() {
    local strain="$1"
    local strain_dir="$2"
    local report_file="$strain_dir/strain_analysis_report.txt"
    local matrix_file="$strain_dir/${strain}_count_matrix.tsv"
    
    log_info "Generating analysis report for strain $strain"
    
    cat > "$report_file" << EOF
=======================================
    Strain $strain Transcriptome Analysis Report
=======================================

Analysis Information:
- Program Version: $VERSION
- Strain Name: $strain
- Data Type: $DATA_TYPE
- Analysis Time: $(date)

Sample List:
EOF
    
    # Add sample information
    local samples="${STRAIN_SAMPLES[$strain]}"
    read -ra sample_array <<< "$samples"
    for sample in "${sample_array[@]}"; do
        sample=$(echo "$sample" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        echo "  - $sample" >> "$report_file"
    done
    
    cat >> "$report_file" << EOF

Result Files:
- Soft-linked data: $strain_dir/linked_data/
- Sample results: $strain_dir/[SampleName]/
- featureCounts summary: $strain_dir/featureCounts/
EOF

    if [[ -f "$matrix_file" ]]; then
        local gene_count=$(tail -n +2 "$matrix_file" | wc -l)
        local sample_count=$(head -1 "$matrix_file" | awk '{print NF-1}')
        cat >> "$report_file" << EOF
- Expression matrix: $matrix_file
  · Gene count: $gene_count
  · Sample count: $sample_count
EOF
    fi
    
    # Add FPKM/TPM information
    if [[ -f "$strain_dir/${strain}_fpkm_matrix.tsv" ]]; then
        cat >> "$report_file" << EOF
- FPKM matrix: $strain_dir/${strain}_fpkm_matrix.tsv
EOF
    fi
    
    if [[ -f "$strain_dir/${strain}_tpm_matrix.tsv" ]]; then
        cat >> "$report_file" << EOF
- TPM matrix: $strain_dir/${strain}_tpm_matrix.tsv
EOF
    fi
    
    # Add DEG analysis results
    if [[ -d "$strain_dir/DEG_results" ]]; then
        cat >> "$report_file" << EOF

Differential Expression Analysis Results:
- Result directory: $strain_dir/DEG_results
EOF
        
        if [[ -f "$strain_dir/DEG_results/DEG_up.csv" ]] && [[ -f "$strain_dir/DEG_results/DEG_down.csv" ]]; then
            local up_count=$(tail -n +2 "$strain_dir/DEG_results/DEG_up.csv" | wc -l)
            local down_count=$(tail -n +2 "$strain_dir/DEG_results/DEG_down.csv" | wc -l)
            cat >> "$report_file" << EOF
  · Up-regulated genes: $up_count
  · Down-regulated genes: $down_count
EOF
        fi
    fi
    
    # Add gene family extraction results
    if [[ -n "$GENEFAMILY_PREFIX" ]] && [[ -d "$strain_dir/${GENEFAMILY_PREFIX}_matrix" ]]; then
        cat >> "$report_file" << EOF

Gene Family Extraction Results ($GENEFAMILY_PREFIX):
- Result directory: $strain_dir/${GENEFAMILY_PREFIX}_matrix
EOF
        
        if [[ -f "$strain_dir/${GENEFAMILY_PREFIX}_matrix/${GENEFAMILY_PREFIX}_count_matrix.tsv" ]]; then
            local gf_gene_count=$(tail -n +2 "$strain_dir/${GENEFAMILY_PREFIX}_matrix/${GENEFAMILY_PREFIX}_count_matrix.tsv" | wc -l)
            cat >> "$report_file" << EOF
  · Extracted gene count: $gf_gene_count
EOF
        fi
    fi

    cat >> "$report_file" << EOF

=======================================
Strain Analysis Completed!
=======================================
EOF
}

# Batch process all strains
batch_process_strains() {
    log_step "Starting batch processing of ${#STRAINS[@]} strains"
    
    local success_count=0
    local failed_count=0
    local failed_strains=()
    
    local total=${#STRAINS[@]}
    local current=0
    
    for strain in "${STRAINS[@]}"; do
        current=$((current + 1))
        
        echo "" >&2
        echo "========================================================================" >&2
        
        if process_strain "$strain" "$current" "$total"; then
            success_count=$((success_count + 1))
        else
            failed_count=$((failed_count + 1))
            failed_strains+=("$strain")
        fi
        
        echo "========================================================================" >&2
        echo "" >&2
    done
    
    # Output summary information
    log_step "Batch processing completed"
    log_info "Total strains: $total"
    log_info "Successful: $success_count"
    log_info "Failed: $failed_count"
    
    if [[ $failed_count -gt 0 ]]; then
        log_warn "Failed strains:"
        for strain in "${failed_strains[@]}"; do
            log_warn "  - $strain"
        done
    fi
    
    # Return number of successfully processed strains
    echo "$success_count"
}

# Generate pan-transcriptome summary report
generate_pan_summary() {
    log_step "Generating pan-transcriptome summary report"
    
    local summary_dir="$OUTPUT_DIR/pan_summary"
    mkdir -p "$summary_dir"
    
    local summary_file="$summary_dir/all_strains_summary.txt"
    local comparison_file="$summary_dir/strain_comparison.txt"
    
    # Generate overall summary
    cat > "$summary_file" << EOF
=======================================
    Pan-transcriptome Analysis Summary Report
=======================================

Analysis Information:
- Program Version: $VERSION
- Data Type: $DATA_TYPE
- Analysis Time: $(date)
- Total Strains: ${#STRAINS[@]}

Configuration Parameters:
- Transcriptome Directory: $TRANSCRIPTOME_DIR
- Genome Directory: $GENOME_DIR
- Annotation Directory: $ANNOTATION_DIR
- Output Directory: $OUTPUT_DIR
- Threads: $THREADS

Strain List:
EOF
    
    for strain in "${STRAINS[@]}"; do
        echo "  - $strain" >> "$summary_file"
    done
    
    cat >> "$summary_file" << EOF

Strain Result Directories:
EOF
    
    for strain in "${STRAINS[@]}"; do
        echo "  - $OUTPUT_DIR/$strain/" >> "$summary_file"
    done
    
    # Generate strain comparison table
    cat > "$comparison_file" << EOF
=======================================
    Inter-strain Statistical Comparison
=======================================

Strain Name | Sample Count | Gene Count | DEG Up | DEG Down | Gene Family
EOF
    
    echo "---------------------------------------------------------------------" >> "$comparison_file"
    
    for strain in "${STRAINS[@]}"; do
        local matrix_file="$OUTPUT_DIR/$strain/${strain}_count_matrix.tsv"
        
        local sample_count=0
        local gene_count=0
        local deg_up="-"
        local deg_down="-"
        local genefamily="-"
        
        if [[ -f "$matrix_file" ]]; then
            sample_count=$(head -1 "$matrix_file" | awk '{print NF-1}')
            gene_count=$(tail -n +2 "$matrix_file" | wc -l)
        fi
        
        # Check DEG results
        if [[ -f "$OUTPUT_DIR/$strain/DEG_results/DEG_up.csv" ]]; then
            deg_up=$(tail -n +2 "$OUTPUT_DIR/$strain/DEG_results/DEG_up.csv" | wc -l)
        fi
        if [[ -f "$OUTPUT_DIR/$strain/DEG_results/DEG_down.csv" ]]; then
            deg_down=$(tail -n +2 "$OUTPUT_DIR/$strain/DEG_results/DEG_down.csv" | wc -l)
        fi
        
        # Check gene family extraction results
        if [[ -n "$GENEFAMILY_PREFIX" ]] && [[ -f "$OUTPUT_DIR/$strain/${GENEFAMILY_PREFIX}_matrix/${GENEFAMILY_PREFIX}_count_matrix.tsv" ]]; then
            genefamily=$(tail -n +2 "$OUTPUT_DIR/$strain/${GENEFAMILY_PREFIX}_matrix/${GENEFAMILY_PREFIX}_count_matrix.tsv" | wc -l)
        fi
        
        printf "%-15s | %-6s | %-8s | %-7s | %-7s | %-8s\n" \
            "$strain" "$sample_count" "$gene_count" \
            "$deg_up" "$deg_down" "$genefamily" >> "$comparison_file"
    done
    
    cat >> "$comparison_file" << EOF

Notes:
- If strain analysis failed, corresponding statistics may show 0 or -
- "-" indicates the function is not enabled or results were not generated

=======================================
EOF
    
    # Add DEG summary information
    if [[ -d "$OUTPUT_DIR/DEG_results" ]]; then
        cat >> "$summary_file" << EOF

Differentially Expressed Genes (DEG) Summary:
-- DEG Summary Directory: $OUTPUT_DIR/DEG_results/
EOF
        
        local deg_strain_count=$(ls -d "$OUTPUT_DIR/DEG_results"/*/ 2>/dev/null | wc -l)
        echo "-- Strains with DEG: $deg_strain_count" >> "$summary_file"
        
        for strain in "${STRAINS[@]}"; do
            if [[ -d "$OUTPUT_DIR/DEG_results/$strain" ]]; then
                local strain_deg_files=$(ls "$OUTPUT_DIR/DEG_results/$strain"/DEG_*.csv 2>/dev/null | wc -l)
                if [[ $strain_deg_files -gt 0 ]]; then
                    echo "  · $strain: $strain_deg_files DEG files" >> "$summary_file"
                fi
            fi
        done
    fi
    
    # Add gene clustering and matrix integration information
    if [[ -d "$OUTPUT_DIR/$GENE_CLUSTER_OUTPUT_DIR" ]]; then
        cat >> "$summary_file" << EOF

Gene Clustering and Matrix Integration:
-- Output Directory: $OUTPUT_DIR/$GENE_CLUSTER_OUTPUT_DIR/
EOF
        
        if [[ -f "$OUTPUT_DIR/$GENE_CLUSTER_OUTPUT_DIR/pan_count_matrix.tsv" ]]; then
            local pan_gene_count=$(tail -n +2 "$OUTPUT_DIR/$GENE_CLUSTER_OUTPUT_DIR/pan_count_matrix.tsv" | wc -l)
            local pan_sample_count=$(head -1 "$OUTPUT_DIR/$GENE_CLUSTER_OUTPUT_DIR/pan_count_matrix.tsv" | awk '{print NF-1}')
            echo "-- Integrated Count Matrix: $pan_gene_count genes × $pan_sample_count samples" >> "$summary_file"
        fi
        
        if [[ -f "$OUTPUT_DIR/$GENE_CLUSTER_OUTPUT_DIR/pan_fpkm_matrix.tsv" ]]; then
            echo "-- Integrated FPKM Matrix: $OUTPUT_DIR/$GENE_CLUSTER_OUTPUT_DIR/pan_fpkm_matrix.tsv" >> "$summary_file"
        fi
        
        if [[ -f "$OUTPUT_DIR/$GENE_CLUSTER_OUTPUT_DIR/pan_tpm_matrix.tsv" ]]; then
            echo "-- Integrated TPM Matrix: $OUTPUT_DIR/$GENE_CLUSTER_OUTPUT_DIR/pan_tpm_matrix.tsv" >> "$summary_file"
        fi
        
        if [[ -d "$OUTPUT_DIR/$GENE_CLUSTER_OUTPUT_DIR/deg" ]]; then
            local deg_strains=$(ls -d "$OUTPUT_DIR/$GENE_CLUSTER_OUTPUT_DIR/deg"/*/ 2>/dev/null | wc -l)
            if [[ $deg_strains -gt 0 ]]; then
                echo "-- Clustered DEG Results: $deg_strains strains" >> "$summary_file"
            fi
        fi
    fi
    
    # Add SV pre-analysis information
    if [[ -d "$OUTPUT_DIR/sv_pre" ]]; then
        cat >> "$summary_file" << EOF

SV Pre-analysis (Control Group Average Expression):
-- Output Directory: $OUTPUT_DIR/sv_pre/
EOF
        
        if [[ -f "$OUTPUT_DIR/sv_pre/pan_control_count.tsv" ]]; then
            local sv_gene_count=$(tail -n +2 "$OUTPUT_DIR/sv_pre/pan_control_count.tsv" | wc -l)
            local sv_strain_count=$(head -1 "$OUTPUT_DIR/sv_pre/pan_control_count.tsv" | awk '{print NF-1}')
            echo "-- Control Average Count Matrix: $sv_gene_count genes × $sv_strain_count strains" >> "$summary_file"
        fi
        
        if [[ -f "$OUTPUT_DIR/sv_pre/pan_control_fpkm.tsv" ]]; then
            echo "-- Control Average FPKM Matrix: $OUTPUT_DIR/sv_pre/pan_control_fpkm.tsv" >> "$summary_file"
        fi
        
        if [[ -f "$OUTPUT_DIR/sv_pre/pan_control_tpm.tsv" ]]; then
            echo "-- Control Average TPM Matrix: $OUTPUT_DIR/sv_pre/pan_control_tpm.tsv" >> "$summary_file"
        fi
    fi
    
    log_info "Pan-transcriptome summary report generated"
    log_info "  - Overall summary: $summary_file"
    log_info "  - Strain comparison: $comparison_file"
}


# Main program
main() {
    log_info "Starting pan-transcriptome analysis - $PROGRAM_NAME v$VERSION"
    
    # Create output directory
    mkdir -p "$OUTPUT_DIR"
    
    # Set log file
    LOG_FILE="$OUTPUT_DIR/pan_analysis.log"
    exec &> >(tee -a "$LOG_FILE")
    
    # Auto-detect script directory
    local script_dir=$(cd "$(dirname "$0")" && pwd)
    log_info "Script directory: $script_dir"
    
    # Auto-detect RNAseq.sh
    RNASEQ_SCRIPT="$script_dir/RNAseq.sh"
    if [[ ! -f "$RNASEQ_SCRIPT" ]]; then
        log_error "RNAseq.sh script does not exist: $RNASEQ_SCRIPT"
        log_error "Please ensure RNAseq.sh and PanRNAseq.sh are in the same directory"
        exit 1
    fi
    log_info "Using RNAseq.sh script: $RNASEQ_SCRIPT"
    
    
    # Validate configuration
    validate_config
    
    # Parse strain list
    parse_strain_list
    
    # Batch process all strains
    local success_count=$(batch_process_strains)
    
    if [[ $success_count -gt 0 ]]; then
        # Collect DEG results
        collect_deg_results
        
        # Gene clustering and matrix integration
        cluster_and_integrate_matrices
        
        # SV_pre analysis (depends on gene_cluster, independent of DEG)
        generate_sv_pre_matrices
        
        # Generate summary report
        generate_pan_summary
    else
        log_error "No strains processed successfully, skipping summary steps"
    fi
    
    echo "" >&2
    log_info "Pan-transcriptome analysis completed! Results saved in: $OUTPUT_DIR"
    log_info "View summary report: $OUTPUT_DIR/pan_summary/all_strains_summary.txt"
    log_info "View strain comparison: $OUTPUT_DIR/pan_summary/strain_comparison.txt"
    
    for strain in "${STRAINS[@]}"; do
        if [[ -f "$OUTPUT_DIR/$strain/gene_count_matrix.tsv" ]]; then
            log_info "Strain $strain expression matrix: $OUTPUT_DIR/$strain/gene_count_matrix.tsv"
        fi
    done
    echo "" >&2
}

# ============================================================================
# Parameter parsing
# ============================================================================

# Default parameters
CONFIG_FILE=""
STRAIN_LIST=""
TRANSCRIPTOME_DIR=""
GENOME_DIR=""
ANNOTATION_DIR=""
OUTPUT_DIR=""
DATA_TYPE="PE"
THREADS=8
KEEP_INTERMEDIATE="false"
QC_BEFORE="true"
QC_AFTER="true"
TRIMMING="true"
TRIM_QUALITY="25"
TRIM_LENGTH="35"
FASTQC_OPTIONS=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -c)
            CONFIG_FILE="$2"
            shift 2
            ;;
        -h|--help)
            print_usage
            exit 0
            ;;
        -v|--version)
            echo "$PROGRAM_NAME version $VERSION"
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

if [[ -z "$CONFIG_FILE" ]]; then
    log_error "Missing required parameter: -c <configuration file>"
    echo ""
    print_usage
    exit 1
fi

if [[ ! -f "$CONFIG_FILE" ]]; then
    log_error "Configuration file does not exist: $CONFIG_FILE"
    exit 1
fi

parse_yaml "$CONFIG_FILE"

main

