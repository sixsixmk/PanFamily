#!/bin/bash

# ============================================================================
# PanMEME.sh - 泛基因组 Motif 分析模块
# 对多个品系的基因家族序列进行批量 MEME motif 分析
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

# 显示帮助信息
show_help() {
    cat << EOF
PanMEME.sh - 泛基因组 Motif 分析模块

用法:
    bash PanMEME.sh -c <config.yaml>

配置文件参数 (YAML):
    strain_list      - 品系名称列表文件，每行一个品系名 (必需)
    input_dir        - 输入序列文件目录 (必需)
    input_pattern    - 输入文件命名模式，使用 {strain} 作为占位符
                       例如: "{strain}.tps.cds.fa" 或 "{strain}_protein.fa"
    output_dir       - 输出目录 (必需)
    seq_type         - 序列类型: protein 或 dna (默认: protein)
    nmotifs          - 每个品系发现的 motif 数量 (默认: 10)
    minw             - motif 最小宽度 (默认: 6)
    maxw             - motif 最大宽度 (默认: 50)
    meme_options     - MEME 额外参数 (可选)
    threads          - 并行分析的品系数 (默认: 4)
    merge_results    - 是否合并所有品系的结果 (默认: true)

命令行参数:
    -c, --config     配置文件路径 (必需)
    -h, --help       显示帮助信息

目录结构示例:
    input_dir/
    ├── Strain1.tps.cds.fa
    ├── Strain2.tps.cds.fa
    ├── Strain3.tps.cds.fa
    └── ...

    strain_list.txt:
    Strain1
    Strain2
    Strain3

输出结构:
    output_dir/
    ├── Strain1_meme/
    │   ├── meme.html
    │   ├── meme.txt
    │   └── ...
    ├── Strain2_meme/
    │   └── ...
    ├── merged_motifs/
    │   ├── all_motifs_summary.txt
    │   └── motif_comparison.txt
    └── pan_motif_report.txt

示例:
    bash PanMEME.sh -c panmotif.yaml

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
    INPUT_PATTERN=$(grep -E "^input_pattern:" "$yaml_file" | sed 's/input_pattern:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "{strain}.fa")
    OUTPUT_DIR=$(grep -E "^output_dir:" "$yaml_file" | sed 's/output_dir:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "")
    SEQ_TYPE=$(grep -E "^seq_type:" "$yaml_file" | sed 's/seq_type:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "protein")
    NMOTIFS=$(grep -E "^nmotifs:" "$yaml_file" | sed 's/nmotifs:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "10")
    MINW=$(grep -E "^minw:" "$yaml_file" | sed 's/minw:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "6")
    MAXW=$(grep -E "^maxw:" "$yaml_file" | sed 's/maxw:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "50")
    MEME_OPTIONS=$(grep -E "^meme_options:" "$yaml_file" | sed 's/meme_options:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "")
    THREADS=$(grep -E "^threads:" "$yaml_file" | sed 's/threads:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "4")
    MERGE_RESULTS=$(grep -E "^merge_results:" "$yaml_file" | sed 's/merge_results:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo "true")
}

# 检查依赖
check_dependencies() {
    log_info "检查依赖软件..."
    
    if ! command -v meme &> /dev/null; then
        log_error "未找到 MEME 软件，请先安装 MEME Suite"
        log_error "安装方法: conda install -c bioconda meme"
        exit 1
    fi
    
    if ! command -v parallel &> /dev/null; then
        log_warn "未找到 GNU parallel，将使用串行处理"
        USE_PARALLEL=false
    else
        USE_PARALLEL=true
    fi
    
    log_info "MEME 版本: $(meme -version 2>&1 | head -1 || echo 'unknown')"
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
    
    # 检查序列类型
    if [[ "$SEQ_TYPE" != "protein" && "$SEQ_TYPE" != "dna" ]]; then
        log_warn "未知序列类型 '$SEQ_TYPE'，使用默认值 'protein'"
        SEQ_TYPE="protein"
    fi
    
    # 统计品系数量
    local strain_count=$(wc -l < "$STRAIN_LIST")
    log_info "品系数量: $strain_count"
}

# 单个品系的 MEME 分析
run_single_meme() {
    local strain="$1"
    local input_file="${INPUT_DIR}/${INPUT_PATTERN/\{strain\}/$strain}"
    local output_subdir="${OUTPUT_DIR}/${strain}_meme"
    
    # 检查输入文件
    if [[ ! -f "$input_file" ]]; then
        log_warn "[$strain] 输入文件不存在: $input_file，跳过"
        return 1
    fi
    
    # 检查序列数量
    local seq_count=$(grep -c "^>" "$input_file" || echo "0")
    if [[ "$seq_count" -lt 2 ]]; then
        log_warn "[$strain] 序列数量不足 ($seq_count)，MEME 需要至少 2 条序列，跳过"
        return 1
    fi
    
    log_info "[$strain] 开始分析 ($seq_count 条序列)..."
    
    # 构建 MEME 命令
    local meme_cmd="meme $input_file -oc $output_subdir -nmotifs $NMOTIFS -minw $MINW -maxw $MAXW"
    
    if [[ "$SEQ_TYPE" == "protein" ]]; then
        meme_cmd="$meme_cmd -protein"
    else
        meme_cmd="$meme_cmd -dna"
    fi
    
    if [[ -n "$MEME_OPTIONS" ]]; then
        meme_cmd="$meme_cmd $MEME_OPTIONS"
    fi
    
    # 运行 MEME
    if eval "$meme_cmd" > "${output_subdir}.log" 2>&1; then
        log_info "[$strain] 分析完成"
        return 0
    else
        log_error "[$strain] 分析失败，查看日志: ${output_subdir}.log"
        return 1
    fi
}

# 导出函数供 parallel 使用
export -f run_single_meme log_info log_warn log_error

# 批量运行 MEME
run_batch_meme() {
    log_info "=========================================="
    log_info "开始泛基因组 Motif 分析"
    log_info "=========================================="
    log_info "品系列表: $STRAIN_LIST"
    log_info "输入目录: $INPUT_DIR"
    log_info "输入模式: $INPUT_PATTERN"
    log_info "输出目录: $OUTPUT_DIR"
    log_info "序列类型: $SEQ_TYPE"
    log_info "Motif 数量: $NMOTIFS"
    log_info "Motif 宽度: $MINW - $MAXW"
    log_info "并行数: $THREADS"
    
    # 创建输出目录
    mkdir -p "$OUTPUT_DIR"
    
    # 导出必要变量供子进程使用
    export INPUT_DIR INPUT_PATTERN OUTPUT_DIR SEQ_TYPE NMOTIFS MINW MAXW MEME_OPTIONS
    
    local success_count=0
    local fail_count=0
    local total_count=$(wc -l < "$STRAIN_LIST")
    
    if [[ "$USE_PARALLEL" == true && "$THREADS" -gt 1 ]]; then
        log_info "使用 GNU parallel 并行处理 ($THREADS 个进程)"
        
        # 使用 parallel 并行处理
        cat "$STRAIN_LIST" | parallel -j "$THREADS" --halt never \
            "run_single_meme {}"
        
        # 统计成功/失败数
        while IFS= read -r strain; do
            local output_subdir="${OUTPUT_DIR}/${strain}_meme"
            if [[ -f "${output_subdir}/meme.html" ]]; then
                ((success_count++))
            else
                ((fail_count++))
            fi
        done < "$STRAIN_LIST"
    else
        log_info "使用串行处理"
        
        while IFS= read -r strain; do
            if run_single_meme "$strain"; then
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
    
    local merged_dir="${OUTPUT_DIR}/merged_motifs"
    mkdir -p "$merged_dir"
    
    local summary_file="${merged_dir}/all_motifs_summary.txt"
    local comparison_file="${merged_dir}/motif_comparison.csv"
    
    # 创建汇总文件头
    echo "Strain,Motif_ID,Width,Sites,E-value,Consensus" > "$comparison_file"
    
    cat > "$summary_file" << EOF
================================================================================
Pan-genome Motif Analysis Summary
================================================================================
Generated: $(date '+%Y-%m-%d %H:%M:%S')

Strains Analyzed:
EOF
    
    # 遍历每个品系的结果
    while IFS= read -r strain; do
        local meme_txt="${OUTPUT_DIR}/${strain}_meme/meme.txt"
        
        if [[ -f "$meme_txt" ]]; then
            echo "  - $strain: SUCCESS" >> "$summary_file"
            
            # 提取 motif 信息
            # 从 meme.txt 中提取 motif 信息
            grep -A5 "^MOTIF" "$meme_txt" 2>/dev/null | while read -r line; do
                if [[ "$line" =~ ^MOTIF ]]; then
                    local motif_info=$(echo "$line" | awk '{print $2}')
                    echo "${strain},${motif_info}" >> "$comparison_file"
                fi
            done || true
        else
            echo "  - $strain: FAILED" >> "$summary_file"
        fi
    done < "$STRAIN_LIST"
    
    cat >> "$summary_file" << EOF

Output Directory: $OUTPUT_DIR
Merged Results: $merged_dir

================================================================================
EOF
    
    log_info "汇总文件: $summary_file"
    log_info "比较文件: $comparison_file"
}

# 生成总报告
generate_report() {
    log_info "生成分析报告..."
    
    local report_file="${OUTPUT_DIR}/pan_motif_report.txt"
    
    cat > "$report_file" << EOF
================================================================================
Pan-genome MEME Motif Analysis Report
================================================================================
Generated: $(date '+%Y-%m-%d %H:%M:%S')

Analysis Parameters:
--------------------
Strain List: $STRAIN_LIST
Input Directory: $INPUT_DIR
Input Pattern: $INPUT_PATTERN
Sequence Type: $SEQ_TYPE
Number of Motifs per Strain: $NMOTIFS
Motif Width Range: $MINW - $MAXW
Parallel Jobs: $THREADS

Results Summary:
----------------
EOF
    
    local success_count=0
    local fail_count=0
    
    while IFS= read -r strain; do
        local output_subdir="${OUTPUT_DIR}/${strain}_meme"
        if [[ -f "${output_subdir}/meme.html" ]]; then
            echo "  [✓] $strain - ${output_subdir}/meme.html" >> "$report_file"
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

Output Structure:
-----------------
$(ls -la "$OUTPUT_DIR" 2>/dev/null | head -20 || echo "No files found")

================================================================================
EOF
    
    log_info "报告已保存至: $report_file"
}

# 主函数
main() {
    # 默认值
    STRAIN_LIST=""
    INPUT_DIR=""
    INPUT_PATTERN="{strain}.fa"
    OUTPUT_DIR=""
    SEQ_TYPE="protein"
    NMOTIFS="10"
    MINW="6"
    MAXW="50"
    MEME_OPTIONS=""
    THREADS="4"
    MERGE_RESULTS="true"
    USE_PARALLEL=false
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
    
    # 批量运行 MEME
    run_batch_meme
    
    # 合并结果
    merge_results
    
    # 生成报告
    generate_report
    
    log_info "=========================================="
    log_info "泛基因组 Motif 分析完成!"
    log_info "结果目录: $OUTPUT_DIR"
    log_info "=========================================="
}

# 错误处理
trap 'log_error "脚本执行失败，行号: $LINENO"' ERR

main "$@"

