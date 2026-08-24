#!/bin/bash

# ============================================================================
# PanProteinProperty.sh - 泛基因组蛋白理化性质分析模块
# 对多个品系的蛋白序列进行批量理化性质分析
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
PanProteinProperty.sh - 泛基因组蛋白理化性质分析模块

用法:
    bash PanProteinProperty.sh -c <config.yaml>

配置文件参数 (YAML):
    strain_list      - 品系名称列表文件，每行一个品系名 (必需)
    input_dir        - 输入序列文件目录 (必需)
    input_pattern    - 输入文件命名模式，使用 {strain} 作为占位符
                       例如: "{strain}.protein.fa" 或 "{strain}.pep.fa"
    output_dir       - 输出目录 (必需)
    threads          - 并行处理的品系数 (可选，默认: 4)
    merge_results    - 是否合并所有品系的结果 (可选，默认: true)

命令行参数:
    -c, --config     配置文件路径 (必需)
    -h, --help       显示帮助信息

目录结构示例:
    input_dir/
    ├── Strain1.protein.fa
    ├── Strain2.protein.fa
    ├── Strain3.protein.fa
    └── ...

    strain_list.txt:
    Strain1
    Strain2
    Strain3

输出结构:
    output_dir/
    ├── Strain1_properties.csv
    ├── Strain2_properties.csv
    ├── ...
    ├── merged_properties.csv      # 合并的结果
    └── pan_protein_report.txt     # 分析报告

示例:
    bash PanProteinProperty.sh -c panprotein.yaml

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
    
    log_info "解析配置文件: $yaml_file"
    
    # 解析 YAML 参数
    STRAIN_LIST=$(grep -E "^strain_list:" "$yaml_file" | sed 's/strain_list:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "")
    INPUT_DIR=$(grep -E "^input_dir:" "$yaml_file" | sed 's/input_dir:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "")
    INPUT_PATTERN=$(grep -E "^input_pattern:" "$yaml_file" | sed 's/input_pattern:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "{strain}.protein.fa")
    OUTPUT_DIR=$(grep -E "^output_dir:" "$yaml_file" | sed 's/output_dir:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "")
    THREADS=$(grep -E "^threads:" "$yaml_file" | sed 's/threads:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "4")
    MERGE_RESULTS=$(grep -E "^merge_results:" "$yaml_file" | sed 's/merge_results:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "true")
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
    
    # 检查 parallel
    if ! command -v parallel &> /dev/null; then
        log_warn "未找到 GNU parallel，将使用串行处理"
        USE_PARALLEL=false
    else
        USE_PARALLEL=true
    fi
    
    log_info "依赖检查通过"
}

# 验证输入
validate_inputs() {
    if [[ -z "$STRAIN_LIST" ]]; then
        log_error "未指定品系列表文件 (strain_list)"
        exit 1
    fi
    
    if [[ ! -f "$STRAIN_LIST" ]]; then
        log_error "品系列表文件不存在: $STRAIN_LIST"
        exit 1
    fi
    
    if [[ -z "$INPUT_DIR" ]]; then
        log_error "未指定输入目录 (input_dir)"
        exit 1
    fi
    
    if [[ ! -d "$INPUT_DIR" ]]; then
        log_error "输入目录不存在: $INPUT_DIR"
        exit 1
    fi
    
    if [[ -z "$OUTPUT_DIR" ]]; then
        log_error "未指定输出目录 (output_dir)"
        exit 1
    fi
    
    # 统计品系数量
    local strain_count=$(wc -l < "$STRAIN_LIST")
    log_info "品系数量: $strain_count"
}

# 单个品系的分析
run_single_analysis() {
    local strain="$1"
    local input_file="${INPUT_DIR}/${INPUT_PATTERN/\{strain\}/$strain}"
    local output_csv="${OUTPUT_DIR}/${strain}_properties.csv"
    
    # 检查输入文件
    if [[ ! -f "$input_file" ]]; then
        log_warn "[$strain] 输入文件不存在: $input_file，跳过"
        return 1
    fi
    
    # 检查序列数量
    local seq_count=$(grep -c "^>" "$input_file" || echo "0")
    if [[ "$seq_count" -eq 0 ]]; then
        log_warn "[$strain] 序列数量为 0，跳过"
        return 1
    fi
    
    log_info "[$strain] 开始分析 ($seq_count 条序列)..."
    
    # 运行分析
    if $PYTHON_CMD "$PYTHON_SCRIPT" --fasta "$input_file" --csv "$output_csv" 2>"${output_csv}.log"; then
        log_info "[$strain] 分析完成"
        return 0
    else
        log_error "[$strain] 分析失败，查看日志: ${output_csv}.log"
        return 1
    fi
}

# 导出函数供 parallel 使用
export -f run_single_analysis log_info log_warn log_error

# 批量运行分析
run_batch_analysis() {
    log_info "=========================================="
    log_info "开始泛基因组蛋白理化性质分析"
    log_info "=========================================="
    log_info "品系列表: $STRAIN_LIST"
    log_info "输入目录: $INPUT_DIR"
    log_info "输入模式: $INPUT_PATTERN"
    log_info "输出目录: $OUTPUT_DIR"
    log_info "并行数: $THREADS"
    
    # 创建输出目录
    mkdir -p "$OUTPUT_DIR"
    
    # 导出必要变量供子进程使用
    export INPUT_DIR INPUT_PATTERN OUTPUT_DIR PYTHON_CMD PYTHON_SCRIPT
    
    local success_count=0
    local fail_count=0
    local total_count=$(wc -l < "$STRAIN_LIST")
    
    if [[ "$USE_PARALLEL" == true && "$THREADS" -gt 1 ]]; then
        log_info "使用 GNU parallel 并行处理 ($THREADS 个进程)"
        
        # 使用 parallel 并行处理
        cat "$STRAIN_LIST" | parallel -j "$THREADS" --halt never \
            "run_single_analysis {}"
        
        # 统计成功/失败数
        while IFS= read -r strain; do
            local output_csv="${OUTPUT_DIR}/${strain}_properties.csv"
            if [[ -f "$output_csv" ]]; then
                ((success_count++))
            else
                ((fail_count++))
            fi
        done < "$STRAIN_LIST"
    else
        log_info "使用串行处理"
        
        while IFS= read -r strain; do
            if run_single_analysis "$strain"; then
                ((success_count++))
            else
                ((fail_count++))
            fi
        done < "$STRAIN_LIST"
    fi
    
    log_info "=========================================="
    log_info "批量分析完成"
    log_info "成功: $success_count / $total_count"
    log_info "失败: $fail_count / $total_count"
    log_info "=========================================="
}

# 合并结果
merge_results() {
    if [[ "$MERGE_RESULTS" != "true" ]]; then
        log_info "跳过结果合并"
        return
    fi
    
    log_info "合并分析结果..."
    
    local merged_file="${OUTPUT_DIR}/merged_properties.csv"
    local first_file=true
    
    # 写入表头（添加 Strain 列）
    echo "Strain,ID,Protein_Length,Molecular_Weight_kDa,Isoelectric_Point,Hydrophilicity,Aliphatic_Index,Instability_Index" > "$merged_file"
    
    # 遍历每个品系的结果
    while IFS= read -r strain; do
        local csv_file="${OUTPUT_DIR}/${strain}_properties.csv"
        
        if [[ -f "$csv_file" ]]; then
            # 跳过表头，添加品系名列
            tail -n +2 "$csv_file" | while IFS= read -r line; do
                echo "${strain},${line}" >> "$merged_file"
            done
        fi
    done < "$STRAIN_LIST"
    
    local total_records=$(($(wc -l < "$merged_file") - 1))
    log_info "合并完成，共 $total_records 条记录"
    log_info "合并文件: $merged_file"
}

# 生成统计摘要
generate_summary() {
    log_info "生成统计摘要..."
    
    local merged_file="${OUTPUT_DIR}/merged_properties.csv"
    local summary_file="${OUTPUT_DIR}/properties_summary.csv"
    
    if [[ ! -f "$merged_file" ]]; then
        log_warn "合并文件不存在，跳过统计摘要"
        return
    fi
    
    # 使用 awk 计算每个品系的统计信息
    cat > "$summary_file" << EOF
Strain,Protein_Count,Avg_Length,Avg_MW_kDa,Avg_pI,Avg_GRAVY,Avg_Aliphatic,Avg_Instability
EOF
    
    tail -n +2 "$merged_file" | awk -F',' '
    {
        strain[$1]++
        length_sum[$1] += $3
        mw_sum[$1] += $4
        pi_sum[$1] += $5
        gravy_sum[$1] += $6
        aliphatic_sum[$1] += $7
        instability_sum[$1] += $8
    }
    END {
        for (s in strain) {
            n = strain[s]
            printf "%s,%d,%.2f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
                s, n,
                length_sum[s]/n,
                mw_sum[s]/n,
                pi_sum[s]/n,
                gravy_sum[s]/n,
                aliphatic_sum[s]/n,
                instability_sum[s]/n
        }
    }
    ' | sort >> "$summary_file"
    
    log_info "统计摘要: $summary_file"
}

# 生成分析报告
generate_report() {
    log_info "生成分析报告..."
    
    local report_file="${OUTPUT_DIR}/pan_protein_report.txt"
    
    cat > "$report_file" << EOF
================================================================================
Pan-genome Protein Physicochemical Properties Analysis Report
================================================================================
Generated: $(date '+%Y-%m-%d %H:%M:%S')

Analysis Parameters:
--------------------
Strain List: $STRAIN_LIST
Input Directory: $INPUT_DIR
Input Pattern: $INPUT_PATTERN
Parallel Jobs: $THREADS

Results Summary:
----------------
EOF
    
    local success_count=0
    local fail_count=0
    
    while IFS= read -r strain; do
        local output_csv="${OUTPUT_DIR}/${strain}_properties.csv"
        if [[ -f "$output_csv" ]]; then
            local count=$(($(wc -l < "$output_csv") - 1))
            echo "  [✓] $strain - $count proteins" >> "$report_file"
            ((success_count++))
        else
            echo "  [✗] $strain - FAILED" >> "$report_file"
            ((fail_count++))
        fi
    done < "$STRAIN_LIST"
    
    cat >> "$report_file" << EOF

Statistics:
-----------
Total Strains: $((success_count + fail_count))
Successful: $success_count
Failed: $fail_count

Output Files:
-------------
- <strain>_properties.csv : Individual strain results
- merged_properties.csv   : Combined results (all strains)
- properties_summary.csv  : Statistical summary per strain

Column Descriptions:
--------------------
- Strain                : Strain name
- ID                    : Protein identifier
- Protein_Length        : Protein length in amino acids
- Molecular_Weight_kDa  : Molecular weight in kilodaltons
- Isoelectric_Point     : Isoelectric point (pI)
- Hydrophilicity        : Grand average of hydropathicity (GRAVY)
- Aliphatic_Index       : Aliphatic index
- Instability_Index     : Instability index (<40 stable, >40 unstable)

================================================================================
EOF
    
    log_info "报告已保存至: $report_file"
}

# 主函数
main() {
    # 默认值
    STRAIN_LIST=""
    INPUT_DIR=""
    INPUT_PATTERN="{strain}.protein.fa"
    OUTPUT_DIR=""
    THREADS="4"
    MERGE_RESULTS="true"
    USE_PARALLEL=false
    CONFIG_FILE=""
    PYTHON_CMD=""
    
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
            *)
                log_error "未知参数: $1"
                show_help
                ;;
        esac
    done
    
    # 必须提供配置文件
    if [[ -z "$CONFIG_FILE" ]]; then
        log_error "请提供配置文件: -c <config.yaml>"
        show_help
    fi
    
    # 解析配置文件
    parse_yaml "$CONFIG_FILE"
    
    # 检查依赖
    check_dependencies
    
    # 验证输入
    validate_inputs
    
    # 批量运行分析
    run_batch_analysis
    
    # 合并结果
    merge_results
    
    # 生成统计摘要
    generate_summary
    
    # 生成报告
    generate_report
    
    log_info "=========================================="
    log_info "泛基因组蛋白理化性质分析完成!"
    log_info "结果目录: $OUTPUT_DIR"
    log_info "=========================================="
}

# 错误处理
trap 'log_error "脚本执行失败，行号: $LINENO"' ERR

main "$@"

