#!/bin/bash

# ============================================================================
# Visualize.sh - Auto-matching Visualization Dispatcher
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
VIS_DIR="${SCRIPT_DIR}/../Visualization"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
BOLD='\033[1m'
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

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

print_banner() {
    echo -e "${CYAN}"
    echo "========================================================================"
    echo "  __     ___                 _ _         "
    echo "  \ \   / (_)___ _   _  __ _| (_)_______ "
    echo "   \ \ / /| / __| | | |/ _\` | | |_  / _ \\"
    echo "    \ V / | \__ \ |_| | (_| | | |/ /  __/"
    echo "     \_/  |_|___/\__,_|\__,_|_|_/___\___|"
    echo "                                          "
    echo "        Auto-matching Visualization Module"
    echo "========================================================================"
    echo -e "${NC}"
}

print_usage() {
    cat << EOF
${BOLD}USAGE:${NC}
    bash Visualize.sh -i <input_dir> -o <output_dir> [options]
    Or call via PanFamily.sh:
    bash PanFamily.sh -V -i <input_dir> -o <output_dir> [options]

${BOLD}DESCRIPTION:${NC}
    Auto-visualization module - Scans R scripts, finds matching files and executes

${BOLD}SUPPORTED VISUALIZATIONS:${NC}
    cis_data           Cis-regulatory element heatmap (cis_data.csv)
    expression         Expression heatmap (expression.csv)
    tree               Phylogenetic tree (tree.nwk)
    treeclass          Annotated phylogenetic tree (treeclass_tree.nwk + treeclass_class.csv)
    pav                PAV heatmap (pav.csv)
    kaks               Ka/Ks analysis plot (kaks.csv)
    enrichment         Enrichment analysis (enrichment.csv)
    variation_strain   Variation-strain statistics (variation_strain.csv)
    variation_genefamily Variation-gene family (variation_genefamily.csv)
    gene_structure     Gene chromosome distribution (gene_structure.csv + chromosome_info.csv)
    collinearity       Collinearity circos plot (collinearity.csv)

${BOLD}OPTIONS:${NC}
    -i, --input <dir>       Input directory containing data files (required)
    -o, --output <dir>      Output directory (required)
    -w, --width <num>       Image width, default 10 inches
    -H, --height <num>      Image height, default 8 inches
    -d, --dpi <num>         Resolution, default 300
    -f, --format <fmt>      Output format: pdf,png,svg,tiff (default: pdf,png)
    --list                  List all available visualization scripts
    -h, --help              Show this help message

${BOLD}EXAMPLES:${NC}
    # Visualize all matching files in directory
    bash Visualize.sh -i ./results/ -o ./plots/

    # Specify image parameters
    bash Visualize.sh -i ./results/ -o ./plots/ -w 12 -H 10 -d 600 -f pdf,png,svg

    # Call via PanFamily.sh
    bash PanFamily.sh -V -i ./results/ -o ./plots/

EOF
}

list_available_scripts() {
    echo -e "\n${BOLD}Available Visualization Scripts:${NC}"
    echo "========================================"
    
    if [[ ! -d "$VIS_DIR" ]]; then
        log_warn "Visualization directory not found: $VIS_DIR"
        return
    fi
    
    local count=0
    for script in "$VIS_DIR"/V_*.R; do
        if [[ -f "$script" && "$(basename "$script")" != "V_template.R" ]]; then
            local script_name=$(basename "$script" .R)
            local file_mapping=$(get_file_mapping "$script_name")
            echo -e "  ${GREEN}●${NC} ${script_name}.R  =>  ${file_mapping}"
            ((count++))
        fi
    done
    
    if [[ $count -eq 0 ]]; then
        echo -e "  ${YELLOW}(No visualization scripts available)${NC}"
    else
        echo ""
        echo "Total $count visualization scripts available"
    fi
    echo "========================================"
}

get_file_mapping() {
    local script_name="$1"
    
    case "$script_name" in
        "V_cis_data")
            echo "cis_data.csv/tsv/xlsx"
            ;;
        "V_expression")
            echo "expression.csv/tsv/xlsx"
            ;;
        "V_tree")
            echo "tree.nwk"
            ;;
        "V_treeclass")
            echo "treeclass_tree.nwk + treeclass_class.csv"
            ;;
        "V_pav")
            echo "pav.csv/tsv/xlsx"
            ;;
        "V_kaks")
            echo "kaks.csv/tsv/xlsx"
            ;;
        "V_enrichment")
            echo "enrichment.csv/tsv/xlsx"
            ;;
        "V_variation_strain")
            echo "variation_strain.csv and/or variation_matrix.csv"
            ;;
        "V_variation_genefamily")
            echo "variation_genefamily.csv/tsv/xlsx"
            ;;
        "V_collinearity")
            echo "collinearity.csv/tsv/xlsx"
            ;;
        "V_gene_structure")
            echo "gene_structure.csv + chromosome_info.csv"
            ;;
        *)
            echo "Unknown mapping"
            ;;
    esac
}

check_input_file_exists() {
    local input_dir="$1"
    local script_name="$2"
    
    case "$script_name" in
        "V_cis_data")
            find_file "$input_dir" "cis_data"
            ;;
        "V_expression")
            find_file "$input_dir" "expression"
            ;;
        "V_tree")
            find_tree_file "$input_dir" "tree"
            ;;
        "V_treeclass")
            find_tree_file "$input_dir" "treeclass_tree"
            ;;
        "V_pav")
            find_file "$input_dir" "pav"
            ;;
        "V_kaks")
            find_file "$input_dir" "kaks"
            ;;
        "V_enrichment")
            find_file "$input_dir" "enrichment"
            ;;
        "V_variation_strain")
            local f1=$(find_file "$input_dir" "variation_strain")
            local f2=$(find_file "$input_dir" "variation_matrix")
            if [[ -n "$f1" || -n "$f2" ]]; then
                echo "$input_dir"
            else
                echo ""
            fi
            ;;
        "V_variation_genefamily")
            find_file "$input_dir" "variation_genefamily"
            ;;
        "V_collinearity")
            find_file "$input_dir" "collinearity"
            ;;
        "V_gene_structure")
            find_file "$input_dir" "gene_structure"
            ;;
        *)
            echo ""
            ;;
    esac
}

find_file() {
    local dir="$1"
    local base_name="$2"
    local extensions=("csv" "tsv" "xlsx" "xls" "txt")
    
    for ext in "${extensions[@]}"; do
        if [[ -f "${dir}/${base_name}.${ext}" ]]; then
            echo "${dir}/${base_name}.${ext}"
            return 0
        fi
    done
    
    echo ""
    return 0
}

find_tree_file() {
    local dir="$1"
    local base_name="$2"
    local extensions=("nwk" "newick" "tree" "treefile" "txt")
    
    for ext in "${extensions[@]}"; do
        if [[ -f "${dir}/${base_name}.${ext}" ]]; then
            echo "${dir}/${base_name}.${ext}"
            return 0
        fi
    done
    
    echo ""
    return 0
}

run_visualization() {
    local input_dir="$1"
    local script_path="$2"
    local output_dir="$3"
    local width="$4"
    local height="$5"
    local dpi="$6"
    local format="$7"
    
    local script_name=$(basename "$script_path" .R)
    
    log_info "Processing: $script_name"
    log_info "  Script: $script_path"
    log_info "  Input:  $input_dir"
    log_info "  Output: $output_dir"
    
    if [[ ! -f "$script_path" ]]; then
        log_error "Script not found: $script_path"
        return 1
    fi
    
    local vis_output="${output_dir}/${script_name}"
    mkdir -p "$vis_output"
    
    if Rscript "$script_path" \
        --input "$input_dir" \
        --output "$vis_output" \
        --width "$width" \
        --height "$height" \
        --dpi "$dpi" \
        --format "$format" 2>&1 | tee "${vis_output}/vis_log.txt"; then
        log_success "Visualization completed: $script_name"
        return 0
    else
        log_error "Visualization failed: $script_name"
        return 1
    fi
}

generate_report() {
    local output_dir="$1"
    local success_count="$2"
    local fail_count="$3"
    local skip_count="$4"
    shift 4
    local results=("$@")
    
    local report_file="${output_dir}/visualization_report.txt"
    
    cat > "$report_file" << EOF
================================================================================
                      PanFamily Visualization Report
================================================================================
Generated: $(date '+%Y-%m-%d %H:%M:%S')
Output Directory: $output_dir

SUMMARY
-------
  ✓ Success: $success_count
  ✗ Failed:  $fail_count
  - Skipped: $skip_count
  Total:     $((success_count + fail_count + skip_count))

DETAILS
-------
EOF
    
    for result in "${results[@]}"; do
        echo "$result" >> "$report_file"
    done
    
    cat >> "$report_file" << EOF

================================================================================
                             End of Report
================================================================================
EOF
    
    log_info "Report saved: $report_file"
}

main() {
    print_banner
    
    local INPUT_DIR=""
    local OUTPUT_DIR=""
    local WIDTH=10
    local HEIGHT=8
    local DPI=300
    local FORMAT="pdf,png"
    
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -i|--input)
                INPUT_DIR="$2"
                shift 2
                ;;
            -o|--output)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            -w|--width)
                WIDTH="$2"
                shift 2
                ;;
            -H|--height)
                HEIGHT="$2"
                shift 2
                ;;
            -d|--dpi)
                DPI="$2"
                shift 2
                ;;
            -f|--format)
                FORMAT="$2"
                shift 2
                ;;
            --list)
                list_available_scripts
                exit 0
                ;;
            -h|--help)
                print_usage
                exit 0
                ;;
            *)
                log_error "Unknown option: $1"
                print_usage
                exit 1
                ;;
        esac
    done
    
    if [[ -z "$INPUT_DIR" ]]; then
        log_error "Input directory is required (-i)"
        print_usage
        exit 1
    fi
    
    if [[ -z "$OUTPUT_DIR" ]]; then
        log_error "Output directory is required (-o)"
        print_usage
        exit 1
    fi
    
    if [[ ! -d "$INPUT_DIR" ]]; then
        log_error "Input directory not found: $INPUT_DIR"
        exit 1
    fi
    
    mkdir -p "$OUTPUT_DIR"
    
    if [[ ! -d "$VIS_DIR" ]]; then
        log_warn "Visualization directory not found, creating: $VIS_DIR"
        mkdir -p "$VIS_DIR"
    fi
    
    log_info "========================================"
    log_info "Starting Auto-Visualization"
    log_info "========================================"
    log_info "Input Directory:  $INPUT_DIR"
    log_info "Output Directory: $OUTPUT_DIR"
    log_info "Script Directory: $VIS_DIR"
    log_info "Image Size:       ${WIDTH}x${HEIGHT} inch @ ${DPI} dpi"
    log_info "Output Format:    $FORMAT"
    log_info "========================================"
    
    local success_count=0
    local fail_count=0
    local skip_count=0
    local results=()
    local processed_count=0
    
    for script in "$VIS_DIR"/V_*.R; do
        if [[ ! -f "$script" ]]; then
            continue
        fi
        
        local script_name=$(basename "$script" .R)
        
        if [[ "$script_name" == "V_template" ]]; then
            continue
        fi
        
        processed_count=$((processed_count + 1))
        
        local input_file=$(check_input_file_exists "$INPUT_DIR" "$script_name")
        
        if [[ -z "$input_file" ]]; then
            local file_mapping=$(get_file_mapping "$script_name")
            log_warn "No input file for: $script_name (expected: $file_mapping)"
            results+=("  - SKIPPED: $script_name (no input file)")
            skip_count=$((skip_count + 1))
            continue
        fi
        
        log_info "Found input for: $script_name"
        
        if run_visualization "$INPUT_DIR" "$script" "$OUTPUT_DIR" \
            "$WIDTH" "$HEIGHT" "$DPI" "$FORMAT"; then
            results+=("  ✓ SUCCESS: $script_name")
            success_count=$((success_count + 1))
        else
            results+=("  ✗ FAILED:  $script_name")
            fail_count=$((fail_count + 1))
        fi
        
        echo ""
    done
    
    if [[ $processed_count -eq 0 ]]; then
        log_warn "No visualization scripts found in $VIS_DIR"
        exit 0
    fi
    
    generate_report "$OUTPUT_DIR" "$success_count" "$fail_count" "$skip_count" "${results[@]}"
    
    echo ""
    log_info "========================================"
    log_info "Visualization Summary"
    log_info "========================================"
    echo -e "  ${GREEN}✓ Success:${NC} $success_count"
    echo -e "  ${RED}✗ Failed:${NC}  $fail_count"
    echo -e "  ${YELLOW}- Skipped:${NC} $skip_count"
    log_info "========================================"
    
    if [[ $fail_count -gt 0 ]]; then
        exit 1
    fi
}

main "$@"
