#!/bin/bash
# Analysis: template enumeration for number of unique precursors (sampled products).
# Run from template_enumeration/; expects 00_run_template_enumeration.sh env (or set LOG_DIR, DATA_DIR, etc.).

# Config and paths
NUM_SAMPLES=10000
INPUT_FILE="$TEMPL_REL_PROCESSED_DATA_PATH/processed_test.csv"

if [ ! -f "$INPUT_FILE" ]; then
    echo "Input file not found, using raw csv file with id and rxn_smiles columns"
    INPUT_FILE="$DATA_DIR/raw_test.csv"
fi

SAMPLED_PRODUCTS_FILE="$DATA_DIR/sampled_products_${NUM_SAMPLES}.csv"
OUTPUT_FILE="$RESULTS_DIR/num_unique_precs_${NUM_SAMPLES}.csv"
LOG_SUBDIR="$LOG_DIR/uniq_precs"
mkdir -p "$LOG_SUBDIR"
TIMESTAMP=$(date +"%y%m%d-%H%Mh")
LOG_FILE="$LOG_SUBDIR/num_unique_precs_${NUM_SAMPLES}_${TIMESTAMP}.log"

# Run
echo "Running template enumeration for unique precursors ($NUM_SAMPLES samples)..."
python templ_enum_uniq_precs.py \
    --template_path "$TEMPL_PATH" \
    --filtered_template_path "$FILTERED_TEMPL_PATH" \
    --input_file "$INPUT_FILE" \
    --sampled_products_file "$SAMPLED_PRODUCTS_FILE" \
    --output_file "$OUTPUT_FILE" \
    --num_samples "$NUM_SAMPLES" \
    --num_cores "$NUM_CORES" \
    --log_file "$LOG_FILE"
echo "Template enumeration for unique precursors completed"
