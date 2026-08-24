#!/bin/bash

# ============================================================================
# ProteinProperty.sh - 单基因组蛋白理化性质分析模块
# 计算蛋白质的分子量、等电点、亲水性、脂肪族指数、不稳定指数等
# ============================================================================

set -euo pipefail

# 日志函数
log_info() {
    echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_error() {
    echo "[ERROR] $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2
}

log_warn() {
    echo "[WARN] $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

# 获取脚本所在目录
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
# Python 分析脚本路径 (与本脚本同目录)
PYTHON_SCRIPT="${SCRIPT_DIR}/protein_analyzer.py"

# 显示帮助信息
show_help() {
    cat << EOF
ProteinProperty.sh - 单基因组蛋白理化性质分析模块

用法:
    bash ProteinProperty.sh -c <config.yaml>
    bash ProteinProperty.sh -i <fasta> -o <output_dir> [options]

配置文件参数 (YAML):
    input_fasta      - 输入的蛋白质 FASTA 文件 (必需)
    output_dir       - 输出目录 (必需)
    output_prefix    - 输出文件前缀 (可选，默认: protein_properties)

命令行参数:
    -c, --config     配置文件路径
    -i, --input      输入蛋白质 FASTA 文件
    -o, --output     输出目录
    -p, --prefix     输出文件前缀
    -h, --help       显示帮助信息

输出文件:
    <prefix>.csv     - 蛋白理化性质表格
    
输出列说明:
    ID                    - 蛋白ID
    Protein_Length        - 蛋白长度 (aa)
    Molecular_Weight_kDa  - 分子量 (kDa)
    Isoelectric_Point     - 等电点 (pI)
    Hydrophilicity        - 亲水性 (GRAVY)
    Aliphatic_Index       - 脂肪族指数
    Instability_Index     - 不稳定指数

示例:
    # 使用配置文件
    bash ProteinProperty.sh -c protein_analysis.yaml

    # 使用命令行参数
    bash ProteinProperty.sh -i proteins.fa -o ./results -p my_proteins

EOF
    exit 0
}

# 解析 YAML 配置文件
parse_yaml() {
    local yaml_file="$1"
    
    if [[ ! -f "$yaml_file" ]]; then
        log_error "配置文件不存在: $yaml_file"
        exit 1
    fi
    
    # 解析 YAML 参数
    INPUT_FASTA=$(grep -E "^input_fasta:" "$yaml_file" | sed 's/input_fasta:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "")
    OUTPUT_DIR=$(grep -E "^output_dir:" "$yaml_file" | sed 's/output_dir:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "")
    OUTPUT_PREFIX=$(grep -E "^output_prefix:" "$yaml_file" | sed 's/output_prefix:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "protein_properties")
}

# 检查依赖
check_dependencies() {
    log_info "检查依赖..."
    
    # 检查 Python
    if ! command -v python3 &> /dev/null && ! command -v python &> /dev/null; then
        log_error "未找到 Python，请先安装 Python3"
        exit 1
    fi
    
    # 确定 Python 命令
    if command -v python3 &> /dev/null; then
        PYTHON_CMD="python3"
    else
        PYTHON_CMD="python"
    fi
    
    # 检查 Biopython
    if ! $PYTHON_CMD -c "from Bio import SeqIO" 2>/dev/null; then
        log_error "未找到 Biopython，请先安装: pip install biopython"
        exit 1
    fi
    
    # 检查 Python 分析脚本
    if [[ ! -f "$PYTHON_SCRIPT" ]]; then
        log_error "未找到分析脚本: $PYTHON_SCRIPT"
        exit 1
    fi
    
    log_info "依赖检查通过"
}

# 验证输入
validate_inputs() {
    if [[ -z "$INPUT_FASTA" ]]; then
        log_error "未指定输入 FASTA 文件"
        exit 1
    fi
    
    if [[ ! -f "$INPUT_FASTA" ]]; then
        log_error "输入文件不存在: $INPUT_FASTA"
        exit 1
    fi
    
    if [[ -z "$OUTPUT_DIR" ]]; then
        log_error "未指定输出目录"
        exit 1
    fi
    
    # 检查序列数量
    local seq_count=$(grep -c "^>" "$INPUT_FASTA" || echo "0")
    if [[ "$seq_count" -eq 0 ]]; then
        log_error "输入文件中没有序列"
        exit 1
    fi
    log_info "输入序列数量: $seq_count"
}

# 运行蛋白理化性质分析
run_analysis() {
    log_info "=========================================="
    log_info "开始蛋白理化性质分析"
    log_info "=========================================="
    log_info "输入文件: $INPUT_FASTA"
    log_info "输出目录: $OUTPUT_DIR"
    log_info "输出前缀: $OUTPUT_PREFIX"
    
    # 创建输出目录
    mkdir -p "$OUTPUT_DIR"
    
    local output_csv="${OUTPUT_DIR}/${OUTPUT_PREFIX}.csv"
    
    # 运行 Python 分析脚本
    log_info "执行分析..."
    $PYTHON_CMD "$PYTHON_SCRIPT" --fasta "$INPUT_FASTA" --csv "$output_csv"
    
    local exit_code=$?
    if [[ $exit_code -ne 0 ]]; then
        log_error "分析失败，退出码: $exit_code"
        exit 1
    fi
    
    # 统计结果
    local result_count=$(($(wc -l < "$output_csv") - 1))
    log_info "分析完成，共 $result_count 条蛋白记录"
    log_info "结果文件: $output_csv"
}

# 生成分析报告
generate_report() {
    log_info "生成分析报告..."
    
    local report_file="${OUTPUT_DIR}/${OUTPUT_PREFIX}_report.txt"
    local output_csv="${OUTPUT_DIR}/${OUTPUT_PREFIX}.csv"
    
    cat > "$report_file" << EOF
================================================================================
Protein Physicochemical Properties Analysis Report
================================================================================
Generated: $(date '+%Y-%m-%d %H:%M:%S')

Input Parameters:
-----------------
Input FASTA: $INPUT_FASTA
Output Directory: $OUTPUT_DIR

Results Summary:
----------------
Total Proteins Analyzed: $(($(wc -l < "$output_csv") - 1))

Output Files:
-------------
- ${OUTPUT_PREFIX}.csv : Protein properties table

Column Descriptions:
--------------------
- ID                    : Protein identifier
- Protein_Length        : Protein length in amino acids
- Molecular_Weight_kDa  : Molecular weight in kilodaltons
- Isoelectric_Point     : Isoelectric point (pI)
- Hydrophilicity        : Grand average of hydropathicity (GRAVY)
                          Negative = hydrophilic, Positive = hydrophobic
- Aliphatic_Index       : Aliphatic index (thermostability indicator)
- Instability_Index     : Instability index
                          <40 = stable, >40 = unstable

================================================================================
EOF
    
    log_info "报告已保存至: $report_file"
}

# 主函数
main() {
    # 默认值
    INPUT_FASTA=""
    OUTPUT_DIR=""
    OUTPUT_PREFIX="protein_properties"
    CONFIG_FILE=""
    
    # 解析命令行参数
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -h|--help)
                show_help
                ;;
            -c|--config)
                CONFIG_FILE="$2"
                shift 2
                ;;
            -i|--input)
                INPUT_FASTA="$2"
                shift 2
                ;;
            -o|--output)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            -p|--prefix)
                OUTPUT_PREFIX="$2"
                shift 2
                ;;
            *)
                log_error "未知参数: $1"
                show_help
                ;;
        esac
    done
    
    # 如果提供了配置文件，解析它
    if [[ -n "$CONFIG_FILE" ]]; then
        parse_yaml "$CONFIG_FILE"
    fi
    
    # 检查依赖
    check_dependencies
    
    # 验证输入
    validate_inputs
    
    # 运行分析
    run_analysis
    
    # 生成报告
    generate_report
    
    log_info "=========================================="
    log_info "蛋白理化性质分析完成!"
    log_info "结果目录: $OUTPUT_DIR"
    log_info "=========================================="
}

# 错误处理
trap 'log_error "脚本执行失败，行号: $LINENO"' ERR

main "$@"

