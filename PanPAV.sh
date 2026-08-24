#!/bin/bash

set -euo pipefail

print_usage() {
    cat << EOF
PanPAV - A simple PAV analysis pipeline controller

USAGE:
    $0 -c config.yaml

OPTIONS:
    -c <FILE>        YAML config file
    -h, --help       Show this help message

CONFIG KEYS (in config.yaml):
    name:                Gene family name/prefix (e.g., "TaTPS")
    representatives:     Path to representative gene sequences FASTA file
    protein_dir:         Directory containing individual strain protein files
    id_file:             File containing strain IDs, one per line
    output_dir:          Directory for all results and intermediate files
    identity_threshold:  Identity threshold for presence (optional, default: 90)
    core_threshold:      Threshold for classifying core genes (optional, default: 0.95)
    blast_threads:       Number of threads for BLAST alignment (optional, default: 8)

NOTE ON PROTEIN FILES:
    Protein files in the 'protein_dir' must follow the naming convention:
    {strain_id}_genefamily_pep.fa, where {strain_id} is from your id_file.
EOF
}


check_dependencies() {
    echo "[INFO] Checking for required tools..."
    local required_tools=("yq" "makeblastdb" "blastp" "Rscript")
    local missing_tools=()

    for tool in "${required_tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            missing_tools+=("$tool")
        fi
    done

    if [ ${#missing_tools[@]} -ne 0 ]; then
        echo "[ERROR] Missing required tools: ${missing_tools[*]}" >&2
        echo "[ERROR] Please install them and ensure they are in your system's PATH." >&2
        exit 1
    fi
    echo "[INFO] All required tools are available."
}

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

CONFIG=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -c)
            CONFIG="$2"; shift 2 ;;
        -h|--help)
            print_usage; exit 0 ;;
        *)
            echo "[ERROR] Unknown argument: $1" >&2
            print_usage; exit 1 ;;
    esac
done

if [[ -z "${CONFIG}" ]]; then
    echo "[ERROR] Missing required argument: -c <config.yaml>" >&2
    print_usage; exit 1
fi

if [[ ! -f "${CONFIG}" ]]; then
    echo "[ERROR] Config file not found: ${CONFIG}" >&2
    exit 1
fi

check_dependencies

echo "[INFO] Starting PanPAV analysis pipeline..."
echo "[INFO] Using config file: ${CONFIG}"


NAME=$(yq -r '.name // ""' "${CONFIG}")
REPRESENTATIVES=$(yq -r '.representatives // ""' "${CONFIG}")
PROTEIN_DIR=$(yq -r '.protein_dir // ""' "${CONFIG}")
ID_FILE=$(yq -r '.id_file // ""' "${CONFIG}")
OUTPUT_DIR=$(yq -r '.output_dir // ""' "${CONFIG}")
IDENTITY_THRESHOLD=$(yq -r '.identity_threshold // 90' "${CONFIG}")
CORE_THRESHOLD=$(yq -r '.core_threshold // 0.95' "${CONFIG}")

error_found=0
[[ -z "${NAME}" || "${NAME}" == "null" ]] && { echo "[ERROR] 'name' is missing in config"; error_found=1; }
[[ -z "${REPRESENTATIVES}" || "${REPRESENTATIVES}" == "null" ]] && { echo "[ERROR] 'representatives' is missing in config"; error_found=1; }
[[ -z "${PROTEIN_DIR}" || "${PROTEIN_DIR}" == "null" ]] && { echo "[ERROR] 'protein_dir' is missing in config"; error_found=1; }
[[ -z "${ID_FILE}" || "${ID_FILE}" == "null" ]] && { echo "[ERROR] 'id_file' is missing in config"; error_found=1; }
[[ -z "${OUTPUT_DIR}" || "${OUTPUT_DIR}" == "null" ]] && { echo "[ERROR] 'output_dir' is missing in config"; error_found=1; }
[[ $error_found -eq 1 ]] && exit 1

[[ ! -f "${REPRESENTATIVES}" ]] && { echo "[ERROR] Representatives file not found: ${REPRESENTATIVES}"; exit 1; }
[[ ! -d "${PROTEIN_DIR}" ]] && { echo "[ERROR] Protein directory not found: ${PROTEIN_DIR}"; exit 1; }
[[ ! -f "${ID_FILE}" ]] && { echo "[ERROR] ID file not found: ${ID_FILE}"; exit 1; }

echo "[INFO] Configuration loaded successfully."

DB_DIR="${OUTPUT_DIR}/blast_db"
BLAST_RESULT_DIR="${OUTPUT_DIR}/blast_results"
PAV_MATRIX_DIR="${OUTPUT_DIR}/pav_matrices"
FINAL_RESULTS_DIR="${OUTPUT_DIR}/final_results"

mkdir -p "${DB_DIR}" "${BLAST_RESULT_DIR}" "${PAV_MATRIX_DIR}" "${FINAL_RESULTS_DIR}"

echo "Step 1: Creating BLAST databases"
echo ""
echo "[STEP 1/5] Creating BLAST databases..."
while IFS= read -r ind; do
    [[ -z "${ind}" || "${ind}" =~ ^[[:space:]]*# ]] && continue
    protein_file="${PROTEIN_DIR}/${ind}_genefamily_pep.fa"
    if [[ ! -f "$protein_file" ]]; then
        echo "[ERROR] Protein file not found for strain '${ind}': ${protein_file}" >&2
        exit 1
    fi
    echo "  -> Processing: ${ind}"
    makeblastdb -in "${protein_file}" -dbtype prot -out "${DB_DIR}/${ind}_${NAME}" &> "${DB_DIR}/${ind}.log"
done < "${ID_FILE}"

echo ""
echo "[STEP 2/5] Running BLASTP alignments..."
while IFS= read -r ind; do
    [[ -z "${ind}" || "${ind}" =~ ^[[:space:]]*# ]] && continue
    echo "  -> Aligning against: ${ind}"
    blastp -query "${REPRESENTATIVES}" \
           -db "${DB_DIR}/${ind}_${NAME}" \
           -outfmt 6 \
           -out "${BLAST_RESULT_DIR}/${ind}.blastp.csv"
done < "${ID_FILE}"

echo ""
echo "[STEP 3/5] Generating individual PAV matrices..."
while IFS= read -r ind; do
    [[ -z "${ind}" || "${ind}" =~ ^[[:space:]]*# ]] && continue
    echo "  -> Generating matrix for: ${ind}"
    Rscript "${SCRIPT_DIR}/pav_analysis.R" \
        -i "${BLAST_RESULT_DIR}/${ind}.blastp.csv" \
        -o "${PAV_MATRIX_DIR}/${ind}_pav.csv" \
        -s "${ind}" \
        -q "${NAME}" \
        -t "${IDENTITY_THRESHOLD}"
done < "${ID_FILE}"

echo ""
echo "[STEP 4/5] Merging all PAV matrices..."
Rscript "${SCRIPT_DIR}/merge_pav.R" \
    -dir "${PAV_MATRIX_DIR}" \
    -o "${FINAL_RESULTS_DIR}/final_merged_pav_matrix.tsv"

echo ""
echo "[STEP 5/5] Classifying genes (core, dispensable, private)..."
Rscript "${SCRIPT_DIR}/pavgene_classification.R" \
    -i "${FINAL_RESULTS_DIR}/final_merged_pav_matrix.tsv" \
    --core_threshold "${CORE_THRESHOLD}" \
    -o "${FINAL_RESULTS_DIR}/gene_categories.tsv"

echo ""
echo "================================================="
echo "[SUCCESS] PanPAV analysis completed."
echo "Results are located in: ${FINAL_RESULTS_DIR}"
echo "Key output files:"
echo "  - final_merged_pav_matrix.tsv (The final combined PAV matrix)"
echo "  - gene_categories.tsv (Core, dispensable, and private gene classification)"
echo "================================================="
