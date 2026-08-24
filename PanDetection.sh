#!/bin/bash

set -euo pipefail

print_usage() {
    cat << EOF
PanDection - Pan-genome BITACORA Orchestrator

USAGE:
    $0 -c pan_genome_config.yaml [--dry-run]

OPTIONS:
    -c <FILE>        YAML config file
    --dry-run        Print resolved commands without executing
    -h, --help       Show help

CONFIG KEYS (pan_genome_config.yaml):
    name                Gene family name (maps to -N for all runs)
    genome_dir          Directory containing all genome files
    gff_dir             Directory containing all GFF files  
    protein_dir         Directory containing all protein files
    querydir_bitacora   FPDB directory for BITACORA (required if bitacora=T)
    evalue              E-value threshold (maps to -e)
    threads             Number of threads
    id_file             File containing strain IDs (one per line)
    output_dir          Output directory for all results
    clean               Optional, export CLEAN (T/F)
    addfilter           Optional, export ADDFILTER (T/F)
    filterlength        Optional, export FILTERLENGTH (int)
    
    # Method selection (at least one must be true)
    bitacora            Optional, run BITACORA analysis (T/F, default: T)
    hmmer               Optional, run HMMER analysis (T/F, default: F)
    blastp              Optional, run BLASTP analysis (T/F, default: F)
    
    # Additional files for HMMER/BLASTP methods
    seed_file           Reference gene family protein file for BLASTP (required if blastp=T)
    hmm_files           HMM files for HMMER, comma-separated (required if hmmer=T)
    evalue_hmm          E-value threshold for HMMER (default: 1e-5)

FILE NAMING CONVENTION:
    Genome files:   {id}.genome.fasta
    GFF files:      {id}.gff3
    Protein files:  {id}.pep.fasta
EOF
}

CONFIG=""
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        -c)
            CONFIG="$2"; shift 2 ;;
        --dry-run)
            DRY_RUN=true; shift ;;
        -h|--help)
            print_usage; exit 0 ;;
        *)
            echo "Unknown argument: $1" >&2
            print_usage; exit 1 ;;
    esac
done

if [[ -z "${CONFIG}" ]]; then
    echo "ERROR: Missing -c <config>" >&2
    print_usage; exit 1
fi

if [[ ! -f "${CONFIG}" ]]; then
    echo "ERROR: Config not found: ${CONFIG}" >&2
    exit 1
fi

if ! command -v yq &>/dev/null; then
    echo "ERROR: yq not found. Please install yq (https://mikefarah.gitbook.io/yq)." >&2
    exit 1
fi

echo "[PanDection] Starting Pan-genome BITACORA analysis..."
echo "[PanDection] Config file: ${CONFIG}"

NAME=$(yq -r '.name // ""' "${CONFIG}")
GENOME_DIR=$(yq -r '.genome_dir // ""' "${CONFIG}")
GFF_DIR=$(yq -r '.gff_dir // ""' "${CONFIG}")
PROTEIN_DIR=$(yq -r '.protein_dir // ""' "${CONFIG}")
QUERYDIR_BITACORA=$(yq -r '.querydir_bitacora // ""' "${CONFIG}")
SCRIPT_DIR=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
EVALUE=$(yq -r '.evalue // "1e-5"' "${CONFIG}")
THREADS=$(yq -r '.threads // "12"' "${CONFIG}")
ID_FILE=$(yq -r '.id_file // ""' "${CONFIG}")
OUTPUT_DIR=$(yq -r '.output_dir // ""' "${CONFIG}")

CLEAN_VAL=$(yq -r '.clean // ""' "${CONFIG}")
ADDFILTER_VAL=$(yq -r '.addfilter // ""' "${CONFIG}")
FILTERLENGTH_VAL=$(yq -r '.filterlength // ""' "${CONFIG}")

RUN_BITACORA=$(yq -r '.bitacora // "T"' "${CONFIG}")
RUN_HMMER=$(yq -r '.hmmer // "F"' "${CONFIG}")
RUN_BLASTP=$(yq -r '.blastp // "F"' "${CONFIG}")

SEED_FILE=$(yq -r '.seed_file // ""' "${CONFIG}")
HMM_FILES_STR=$(yq -r '.hmm_files // ""' "${CONFIG}")
EVALUE_HMM=$(yq -r '.evalue_hmm // "1e-5"' "${CONFIG}")


missing=()
[[ -z "${NAME}" || "${NAME}" == "null" ]] && missing+=("name")
[[ -z "${ID_FILE}" || "${ID_FILE}" == "null" ]] && missing+=("id_file")
[[ -z "${OUTPUT_DIR}" || "${OUTPUT_DIR}" == "null" ]] && missing+=("output_dir")


if [[ "${RUN_BITACORA}" == "T" ]]; then
    [[ -z "${GENOME_DIR}" || "${GENOME_DIR}" == "null" ]] && missing+=("genome_dir")
    [[ -z "${GFF_DIR}" || "${GFF_DIR}" == "null" ]] && missing+=("gff_dir")
    [[ -z "${QUERYDIR_BITACORA}" || "${QUERYDIR_BITACORA}" == "null" ]] && missing+=("querydir_bitacora")
fi

if [[ "${RUN_HMMER}" == "T" || "${RUN_BLASTP}" == "T" ]]; then
    [[ -z "${PROTEIN_DIR}" || "${PROTEIN_DIR}" == "null" ]] && missing+=("protein_dir")
fi

if [[ "${RUN_BLASTP}" == "T" ]]; then
    [[ -z "${SEED_FILE}" || "${SEED_FILE}" == "null" ]] && missing+=("seed_file")
fi

if [[ "${RUN_HMMER}" == "T" ]]; then
    [[ -z "${HMM_FILES_STR}" || "${HMM_FILES_STR}" == "null" ]] && missing+=("hmm_files")
fi

if [[ ${#missing[@]} -gt 0 ]]; then
    echo "ERROR: Missing required config keys: ${missing[*]}" >&2
    exit 1
fi

if [[ "${RUN_BITACORA}" != "T" && "${RUN_HMMER}" != "T" && "${RUN_BLASTP}" != "T" ]]; then
    echo "ERROR: At least one method (bitacora, hmmer, or blastp) must be enabled (set to T)" >&2
    exit 1
fi


if [[ "${RUN_BITACORA}" == "T" ]]; then
    [[ ! -d "${GENOME_DIR}" ]] && { echo "ERROR: genome_dir not found: ${GENOME_DIR}" >&2; exit 1; }
    [[ ! -d "${GFF_DIR}" ]] && { echo "ERROR: gff_dir not found: ${GFF_DIR}" >&2; exit 1; }
    [[ ! -d "${QUERYDIR_BITACORA}" ]] && { echo "ERROR: querydir_bitacora not found: ${QUERYDIR_BITACORA}" >&2; exit 1; }
fi

if [[ "${RUN_HMMER}" == "T" || "${RUN_BLASTP}" == "T" ]]; then
    [[ ! -d "${PROTEIN_DIR}" ]] && { echo "ERROR: protein_dir not found: ${PROTEIN_DIR}" >&2; exit 1; }
fi

if [[ "${RUN_BLASTP}" == "T" ]]; then
    [[ ! -f "${SEED_FILE}" ]] && { echo "ERROR: seed_file not found: ${SEED_FILE}" >&2; exit 1; }
fi

if [[ "${RUN_HMMER}" == "T" ]]; then
    IFS=',' read -r -a HMM_FILES <<< "${HMM_FILES_STR}"
    for hmm_file in "${HMM_FILES[@]}"; do
        [[ ! -f "${hmm_file}" ]] && { echo "ERROR: hmm_file not found: ${hmm_file}" >&2; exit 1; }
    done
fi

[[ ! -f "${ID_FILE}" ]] && { echo "ERROR: id_file not found: ${ID_FILE}" >&2; exit 1; }

BITACORA_SCRIPT="${SCRIPT_DIR}/changeBitacora.sh"
GENE_DETECTION_SCRIPT="${SCRIPT_DIR}/gene_detection.sh"

if [[ "${RUN_BITACORA}" == "T" ]]; then
    [[ ! -x "${BITACORA_SCRIPT}" ]] && { echo "ERROR: BITACORA script not found or not executable: ${BITACORA_SCRIPT}" >&2; exit 1; }
fi
if [[ "${RUN_HMMER}" == "T" || "${RUN_BLASTP}" == "T" ]]; then
    [[ ! -x "${GENE_DETECTION_SCRIPT}" ]] && { echo "ERROR: Gene detection script not found or not executable: ${GENE_DETECTION_SCRIPT}" >&2; exit 1; }
fi


echo "[PanDection] Configuration loaded:"
echo "[PanDection] NAME=${NAME}"
echo "[PanDection] Methods: BITACORA=${RUN_BITACORA}, HMMER=${RUN_HMMER}, BLASTP=${RUN_BLASTP}"
if [[ "${RUN_BITACORA}" == "T" ]]; then
    echo "[PanDection] GENOME_DIR=${GENOME_DIR}"
    echo "[PanDection] GFF_DIR=${GFF_DIR}"
    echo "[PanDection] QUERYDIR_BITACORA=${QUERYDIR_BITACORA}"
fi
if [[ "${RUN_HMMER}" == "T" || "${RUN_BLASTP}" == "T" ]]; then
    echo "[PanDection] PROTEIN_DIR=${PROTEIN_DIR}"
fi
echo "[PanDection] SCRIPT_DIR=${SCRIPT_DIR}"
echo "[PanDection] EVALUE=${EVALUE}"
echo "[PanDection] THREADS=${THREADS}"
echo "[PanDection] ID_FILE=${ID_FILE}"
echo "[PanDection] OUTPUT_DIR=${OUTPUT_DIR}"

# Set optional environment variables for changeBitacora.sh
if [[ -n "${CLEAN_VAL}" && "${CLEAN_VAL}" != "null" ]]; then
    export CLEAN="${CLEAN_VAL}"
    echo "[PanDection] CLEAN=${CLEAN}"
fi
if [[ -n "${ADDFILTER_VAL}" && "${ADDFILTER_VAL}" != "null" ]]; then
    export ADDFILTER="${ADDFILTER_VAL}"
    echo "[PanDection] ADDFILTER=${ADDFILTER}"
fi
if [[ -n "${FILTERLENGTH_VAL}" && "${FILTERLENGTH_VAL}" != "null" ]]; then
    export FILTERLENGTH="${FILTERLENGTH_VAL}"
    echo "[PanDection] FILTERLENGTH=${FILTERLENGTH}"
fi
if [[ -n "${THREADS}" && "${THREADS}" != "null" ]]; then
    export THREADS="${THREADS}"
fi

mkdir -p "${OUTPUT_DIR}"

echo "[PanDection] Processing strain IDs and setting up directories..."

COMMANDS=()
STRAIN_COUNT=0

while IFS= read -r strain_id; do
    [[ -z "${strain_id}" || "${strain_id}" =~ ^[[:space:]]*# ]] && continue
    
    STRAIN_COUNT=$((STRAIN_COUNT + 1))
    echo "[PanDection] Processing strain: ${strain_id}"
    
    if [[ -n "${GENOME_DIR:-}" && "${GENOME_DIR}" != "null" ]]; then
        GENOME_FILE="${GENOME_DIR}/${strain_id}.genome.fasta"
    fi
    if [[ -n "${GFF_DIR:-}" && "${GFF_DIR}" != "null" ]]; then
        GFF_FILE="${GFF_DIR}/${strain_id}.gff3"
    fi
    PROTEIN_FILE="${PROTEIN_DIR}/${strain_id}.pep.fasta"
    
    missing_files=()
    
    if [[ "${RUN_BITACORA}" == "T" ]]; then
        [[ ! -f "${GENOME_FILE}" ]] && missing_files+=("${GENOME_FILE}")
        [[ ! -f "${GFF_FILE}" ]] && missing_files+=("${GFF_FILE}")
    fi
    
    if [[ "${RUN_BITACORA}" == "T" || "${RUN_HMMER}" == "T" || "${RUN_BLASTP}" == "T" ]]; then
        [[ ! -f "${PROTEIN_FILE}" ]] && missing_files+=("${PROTEIN_FILE}")
    fi
    
    if [[ ${#missing_files[@]} -gt 0 ]]; then
        echo "ERROR: Missing files for strain ${strain_id}: ${missing_files[*]}" >&2
        exit 1
    fi
    
    STRAIN_OUTPUT_DIR="${OUTPUT_DIR}/${strain_id}"
    mkdir -p "${STRAIN_OUTPUT_DIR}"
    
    echo "[PanDection] Creating symbolic links for ${strain_id} (skip existing)..."
    
    if [[ -n "${GENOME_FILE:-}" ]]; then
        [[ -e "${STRAIN_OUTPUT_DIR}/${strain_id}.genome.fasta" ]] || ln -s "$(readlink -f "${GENOME_FILE}")" "${STRAIN_OUTPUT_DIR}/${strain_id}.genome.fasta"
    fi
    if [[ -n "${GFF_FILE:-}" ]]; then
        [[ -e "${STRAIN_OUTPUT_DIR}/${strain_id}.gff3" ]] || ln -s "$(readlink -f "${GFF_FILE}")" "${STRAIN_OUTPUT_DIR}/${strain_id}.gff3"
    fi
    [[ -e "${STRAIN_OUTPUT_DIR}/${strain_id}.pep.fasta" ]] || ln -s "$(readlink -f "${PROTEIN_FILE}")" "${STRAIN_OUTPUT_DIR}/${strain_id}.pep.fasta"
    
    if [[ "${RUN_BITACORA}" == "T" ]]; then
        CMD=("${BITACORA_SCRIPT}" 
             -N "${NAME}" 
             -g "${STRAIN_OUTPUT_DIR}/${strain_id}.genome.fasta" 
             -a "${STRAIN_OUTPUT_DIR}/${strain_id}.gff3" 
             -p "${STRAIN_OUTPUT_DIR}/${strain_id}.pep.fasta" 
             -q "${QUERYDIR_BITACORA}" 
             -e "${EVALUE}")
        COMMANDS+=("BITACORA:${CMD[*]}")
    fi
    
done < "${ID_FILE}"

echo "[PanDection] Total strains to process: ${STRAIN_COUNT}"

if [[ ${STRAIN_COUNT} -eq 0 ]]; then
    echo "ERROR: No valid strain IDs found in ${ID_FILE}" >&2
    exit 1
fi

echo "[PanDection] File preparation completed. Ready to run analysis..."

if [[ "${DRY_RUN}" == true ]]; then
    echo "[PanDection] Dry-run mode - Commands that would be executed:"
    for i in "${!COMMANDS[@]}"; do
        echo "[$((i+1))/${STRAIN_COUNT}] ${COMMANDS[$i]}"
    done
    echo "[PanDection] Dry-run completed."
    exit 0
fi

echo "[PanDection] Starting analysis for all strains..."
CURRENT=0

while IFS= read -r strain_id; do
    [[ -z "${strain_id}" || "${strain_id}" =~ ^[[:space:]]*# ]] && continue
    
    CURRENT=$((CURRENT + 1))
    STRAIN_OUTPUT_DIR="${OUTPUT_DIR}/${strain_id}"
    
    echo ""
    echo "=========================================="
    echo "[PanDection] Processing strain ${CURRENT}/${STRAIN_COUNT}: ${strain_id}"
    echo "=========================================="
    
    cd "${STRAIN_OUTPUT_DIR}"
    
    # Run BITACORA if enabled
    if [[ "${RUN_BITACORA}" == "T" ]]; then
        echo "[PanDection] Running BITACORA analysis for ${strain_id}..."
        "${BITACORA_SCRIPT}" \
            -N "${NAME}" \
            -g "${strain_id}.genome.fasta" \
            -a "${strain_id}.gff3" \
            -p "${strain_id}.pep.fasta" \
            -q "${QUERYDIR_BITACORA}" \
            -e "${EVALUE}"
        
        if [[ $? -eq 0 ]]; then
            echo "[PanDection] ✓ BITACORA analysis completed for ${strain_id}"
        else
            echo "[PanDection] ✗ BITACORA analysis failed for ${strain_id}" >&2
            exit 1
        fi
    fi
    
    # Run HMMER/BLASTP if enabled
    if [[ "${RUN_HMMER}" == "T" || "${RUN_BLASTP}" == "T" ]]; then
        echo "[PanDection] Running HMMER/BLASTP analysis for ${strain_id}..."
        
        if [[ "${RUN_BLASTP}" == "T" && "${RUN_HMMER}" == "T" ]]; then
            METHOD_PARAM="both"
        elif [[ "${RUN_BLASTP}" == "T" ]]; then
            METHOD_PARAM="blastp"
        else
            METHOD_PARAM="hmmer"
        fi
        
        GENE_DETECTION_ARGS=(
            --method "${METHOD_PARAM}"
            -i "${strain_id}.pep.fasta"
            -o "${NAME}_gene_detection"
            --threads "${THREADS}"
            --evalue-blast "${EVALUE}"
            --evalue-hmm "${EVALUE_HMM}"
        )
        
        if [[ "${RUN_BLASTP}" == "T" ]]; then
            GENE_DETECTION_ARGS+=(--seed "${SEED_FILE}")
        fi
        if [[ "${RUN_HMMER}" == "T" ]]; then
            GENE_DETECTION_ARGS+=(--hmm "${HMM_FILES_STR}")
        fi
        
        "${GENE_DETECTION_SCRIPT}" "${GENE_DETECTION_ARGS[@]}"
        
        if [[ $? -eq 0 ]]; then
            echo "[PanDection] ✓ HMMER/BLASTP analysis completed for ${strain_id}"
        else
            echo "[PanDection] ✗ HMMER/BLASTP analysis failed for ${strain_id}" >&2
            exit 1
        fi
    fi
    
    cd - > /dev/null
    
done < "${ID_FILE}"

echo ""
echo "=========================================="
echo "[PanDection] Processing and collecting results..."
echo "=========================================="

GENE_FAMILY_DIR="${OUTPUT_DIR}/gene_family"
mkdir -p "${GENE_FAMILY_DIR}"

RESULT_COUNT=0
while IFS= read -r strain_id; do
    [[ -z "${strain_id}" || "${strain_id}" =~ ^[[:space:]]*# ]] && continue
    
    RESULT_COUNT=$((RESULT_COUNT + 1))
    echo "[PanDection] Processing results for strain ${strain_id}..."
    
    if [[ "${RUN_BITACORA}" == "T" ]]; then
        STRAIN_RESULT_DIR="${OUTPUT_DIR}/${strain_id}/$(echo ${NAME} | tr '[:upper:]' '[:lower:]')"
        GFF_RESULT="${STRAIN_RESULT_DIR}/$(echo ${NAME} | tr '[:upper:]' '[:lower:]')_genomic_and_annotated_proteins_trimmed_idseqsclustered.gff3"
        FASTA_RESULT="${STRAIN_RESULT_DIR}/$(echo ${NAME} | tr '[:upper:]' '[:lower:]')_genomic_and_annotated_proteins_trimmed_idseqsclustered.fasta"
        
        if [[ -f "${GFF_RESULT}" && -f "${FASTA_RESULT}" ]]; then
            GENE_LIST_FILE="${GENE_FAMILY_DIR}/${strain_id}_bitacora_genefamily.list"
            echo "[PanDection] Extracting BITACORA gene family members for ${strain_id}..."
            grep "^>" "${FASTA_RESULT}" | sed 's/^>//' | cut -d' ' -f1 > "${GENE_LIST_FILE}"
            
            RENAMED_GFF="${GENE_FAMILY_DIR}/${strain_id}_bitacora_genomic_and_annotated_proteins_trimmed_idseqsclustered.gff3"
            RENAMED_FASTA="${GENE_FAMILY_DIR}/${strain_id}_bitacora_genomic_and_annotated_proteins_trimmed_idseqsclustered.fasta"
            
            echo "[PanDection] Creating symbolic links for BITACORA result files for ${strain_id} (skip existing)..."
            [[ -e "${RENAMED_GFF}" ]] || ln -s "$(readlink -f "${GFF_RESULT}")" "${RENAMED_GFF}"
            [[ -e "${RENAMED_FASTA}" ]] || ln -s "$(readlink -f "${FASTA_RESULT}")" "${RENAMED_FASTA}"
            
            echo "[PanDection] ✓ Processed BITACORA results for ${strain_id}: $(wc -l < "${GENE_LIST_FILE}") genes"
        else
            echo "[PanDection] Warning: BITACORA result files not found for ${strain_id}"
        fi
    fi
    
    if [[ "${RUN_HMMER}" == "T" || "${RUN_BLASTP}" == "T" ]]; then
        GENE_DETECTION_DIR="${OUTPUT_DIR}/${strain_id}/${NAME}_gene_detection"
        GENE_DETECTION_LIST="${GENE_DETECTION_DIR}/family_genes.txt"
        GENE_DETECTION_FASTA="${GENE_DETECTION_DIR}/family_protein.fa"
        
        if [[ "${RUN_HMMER}" == "T" && "${RUN_BLASTP}" == "T" ]]; then
            METHOD_SUFFIX="hmmer_blastp"
        elif [[ "${RUN_HMMER}" == "T" ]]; then
            METHOD_SUFFIX="hmmer"
        elif [[ "${RUN_BLASTP}" == "T" ]]; then
            METHOD_SUFFIX="blastp"
        fi
        
        if [[ -f "${GENE_DETECTION_LIST}" && -f "${GENE_DETECTION_FASTA}" ]]; then
            RENAMED_LIST="${GENE_FAMILY_DIR}/${strain_id}_${METHOD_SUFFIX}_genefamily.list"
            RENAMED_DETECTION_FASTA="${GENE_FAMILY_DIR}/${strain_id}_${METHOD_SUFFIX}_family_protein.fasta"
            
            echo "[PanDection] Creating symbolic links for ${METHOD_SUFFIX^^} result files for ${strain_id} (skip existing)..."
            [[ -e "${RENAMED_LIST}" ]] || ln -s "$(readlink -f "${GENE_DETECTION_LIST}")" "${RENAMED_LIST}"
            [[ -e "${RENAMED_DETECTION_FASTA}" ]] || ln -s "$(readlink -f "${GENE_DETECTION_FASTA}")" "${RENAMED_DETECTION_FASTA}"
            
            echo "[PanDection] ✓ Processed ${METHOD_SUFFIX^^} results for ${strain_id}: $(wc -l < "${RENAMED_LIST}") genes"
        else
            echo "[PanDection] Warning: ${METHOD_SUFFIX^^} result files not found for ${strain_id}"
        fi
    fi
    if [[ "${RUN_BITACORA}" == "T" && "${RUN_HMMER}" == "T" && "${RUN_BLASTP}" == "T" ]]; then
        echo "[PanDection] Integrating BITACORA + HMMER + BLASTP results for ${strain_id}..."
        
        BITACORA_LIST="${GENE_FAMILY_DIR}/${strain_id}_bitacora_genefamily.list"
        BITACORA_FASTA="${GENE_FAMILY_DIR}/${strain_id}_bitacora_genomic_and_annotated_proteins_trimmed_idseqsclustered.fasta"
        HMMER_BLASTP_LIST="${GENE_FAMILY_DIR}/${strain_id}_hmmer_blastp_genefamily.list"
        HMMER_BLASTP_FASTA="${GENE_FAMILY_DIR}/${strain_id}_hmmer_blastp_family_protein.fasta"
        
        INTEGRATED_LIST="${GENE_FAMILY_DIR}/${strain_id}_integrated_genefamily.list"
        INTEGRATED_FASTA="${GENE_FAMILY_DIR}/${strain_id}_integrated_family_protein.fasta"
        
        if [[ -f "${BITACORA_LIST}" && -f "${HMMER_BLASTP_LIST}" && -f "${BITACORA_FASTA}" && -f "${HMMER_BLASTP_FASTA}" ]]; then
            echo "[PanDection] Computing (HMMER ∩ BLASTP) ∪ BITACORA union..."
            cat "${HMMER_BLASTP_LIST}" "${BITACORA_LIST}" | sort -u > "${INTEGRATED_LIST}"
            
            if [[ -s "${INTEGRATED_LIST}" ]]; then
                echo "[PanDection] Extracting protein sequences for integrated results..."
                {
                    seqkit grep -f "${INTEGRATED_LIST}" "${HMMER_BLASTP_FASTA}" 2>/dev/null || true
                    seqkit grep -f "${INTEGRATED_LIST}" "${BITACORA_FASTA}" 2>/dev/null || true
                } | seqkit rmdup -n > "${INTEGRATED_FASTA}"
                
                INTEGRATED_COUNT=$(wc -l < "${INTEGRATED_LIST}")
                echo "[PanDection] ✓ Integrated results for ${strain_id}: ${INTEGRATED_COUNT} genes (HMMER ∩ BLASTP ∪ BITACORA)"
            else
                echo "[PanDection] Warning: No genes found in integrated results for ${strain_id}"
                touch "${INTEGRATED_FASTA}"
            fi
        else
            echo "[PanDection] Warning: Cannot integrate results for ${strain_id} - missing input files"
            echo "  Required: ${BITACORA_LIST}, ${HMMER_BLASTP_LIST}, ${BITACORA_FASTA}, ${HMMER_BLASTP_FASTA}"
        fi
    fi
    
done < "${ID_FILE}"

echo ""
echo "=========================================="
echo "[PanDection] Pan-genome analysis completed successfully!"
echo "[PanDection] Results saved in: ${OUTPUT_DIR}"
echo "[PanDection] Gene family results collected in: ${GENE_FAMILY_DIR}"
echo "[PanDection] Total strains processed: ${STRAIN_COUNT}"
echo "[PanDection] Total result sets collected: ${RESULT_COUNT}"
echo "=========================================="

