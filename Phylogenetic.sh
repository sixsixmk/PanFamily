#!/bin/bash

while getopts "c:h" opt; do
    case $opt in
        c)
            CONFIG_FILE="$OPTARG"
            ;;
        h)
            echo "Usage: $0 -c config_file"
            echo "Options:"
            echo "  -c    Configuration file path"
            echo "  -h    Show this help message"
            exit 0
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            echo "Usage: $0 -c config_file"
            exit 1
            ;;
    esac
done

if [ -z "$CONFIG_FILE" ]; then
    echo "Error: Configuration file is required!"
    echo "Usage: $0 -c config_file"
    exit 1
fi

if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Config file '$CONFIG_FILE' not found!"
    exit 1
fi

source "$CONFIG_FILE"
if [ -z "$INPUT_FILE" ] || [ -z "$OUTPUT_PREFIX" ] || [ -z "$FASTTREE_RUN" ] || [ -z "$TRIMAL_RUN" ]; then
    echo "Error: Required parameters missing in config file!"
    echo "Required: INPUT_FILE, OUTPUT_PREFIX, FASTTREE_RUN, TRIMAL_RUN"
    exit 1
fi

echo "=== Phylogenetic Analysis Pipeline ==="
echo "Input file: $INPUT_FILE"
echo "Output prefix: $OUTPUT_PREFIX"
echo "Output directory: $OUTPUT_DIR"
echo "MAFFT method: auto (automatically select best method)"
echo "Run FastTree: $FASTTREE_RUN"
echo "Run trimAl: $TRIMAL_RUN"
echo

THREADS=${THREADS:-12}

OUTPUT_DIR="./phylogenetic"
mkdir -p "$OUTPUT_DIR"
ALIGN_FILE="${OUTPUT_DIR}/${OUTPUT_PREFIX}_align.fasta"
TRIMMED_FILE="${OUTPUT_DIR}/${OUTPUT_PREFIX}_agl_trimed.fasta"
FASTTREE_FILE="${OUTPUT_DIR}/${OUTPUT_PREFIX}_fasttree.nwk"
IQTREE_FILE="${OUTPUT_DIR}/${OUTPUT_PREFIX}_iqtree.treefile"

echo "Step 1: Multiple sequence alignment with MAFFT"
echo "Running: mafft --auto --thread $THREADS $INPUT_FILE > $ALIGN_FILE"

mafft --auto --thread $THREADS "$INPUT_FILE" > "$ALIGN_FILE"

if [ $? -ne 0 ]; then
    echo "Error: MAFFT alignment failed!"
    exit 1
fi
echo "MAFFT alignment completed: $ALIGN_FILE"
echo

if [ "$TRIMAL_RUN" = "T" ] || [ "$TRIMAL_RUN" = "true" ] || [ "$TRIMAL_RUN" = "TRUE" ]; then
    echo "Step 2: Trimming alignment with trimAl"
    echo "Running: trimal -in $ALIGN_FILE -out $TRIMMED_FILE -automated1"
    
    trimal -in "$ALIGN_FILE" -out "$TRIMMED_FILE" -automated1
    
    if [ $? -ne 0 ]; then
        echo "Error: trimAl failed!"
        exit 1
    fi
    echo "Alignment trimming completed: $TRIMMED_FILE"
    TREE_INPUT_FILE="$TRIMMED_FILE"
else
    echo "Step 2: Skipping trimAl (TRIMAL_RUN=F)"
    echo "Using MAFFT alignment directly for tree building"
    TREE_INPUT_FILE="$ALIGN_FILE"
fi
echo

if [ "$FASTTREE_RUN" = "T" ] || [ "$FASTTREE_RUN" = "true" ] || [ "$FASTTREE_RUN" = "TRUE" ]; then
    echo "Step 3a: Building phylogenetic tree with FastTree"
    echo "Running: FastTree $TRIMMED_FILE > $FASTTREE_FILE"
    
    FastTree "$TREE_INPUT_FILE" > "$FASTTREE_FILE"
    
    if [ $? -ne 0 ]; then
        echo "Error: FastTree failed!"
        exit 1
    fi
    echo "FastTree analysis completed: $FASTTREE_FILE"
    echo
fi

echo "Step 3b: Building phylogenetic tree with IQ-TREE"
echo "Running: iqtree -s $TREE_INPUT_FILE -m mfp -b 1000 -nt auto"

cd "$OUTPUT_DIR"
iqtree -s "$(basename "$TREE_INPUT_FILE")" -m mfp -b 1000 -nt auto
cd ..

if [ $? -ne 0 ]; then
    echo "Error: IQ-TREE failed!"
    exit 1
fi
echo "IQ-TREE analysis completed"
echo

echo "=== Pipeline completed successfully! ==="
echo "Output files (all in $OUTPUT_DIR/):"
echo "  - Alignment: $(basename "$ALIGN_FILE")"
if [ "$TRIMAL_RUN" = "T" ] || [ "$TRIMAL_RUN" = "true" ] || [ "$TRIMAL_RUN" = "TRUE" ]; then
    echo "  - Trimmed alignment: $(basename "$TRIMMED_FILE")"
    TREE_SUFFIX="$(basename "$TREE_INPUT_FILE" .fasta)"
else
    echo "  - Trimmed alignment: Skipped"
    TREE_SUFFIX="$(basename "$TREE_INPUT_FILE" .fasta)"
fi
if [ "$FASTTREE_RUN" = "T" ] || [ "$FASTTREE_RUN" = "true" ] || [ "$FASTTREE_RUN" = "TRUE" ]; then
    echo "  - FastTree: $(basename "$FASTTREE_FILE")"
fi
echo "  - IQ-TREE: ${TREE_SUFFIX}.treefile"
echo "  - Bootstrap support: ${TREE_SUFFIX}.ufboot"