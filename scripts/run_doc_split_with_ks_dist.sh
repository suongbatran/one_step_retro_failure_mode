#!/bin/bash
# Script to run adaptive split with KS-test

set -e 

# Paths
DATASET="uspto"
INPUT_FILE="data_${DATASET}/${DATASET}_preprocessed.json"

FEATURES_FILE="data_${DATASET}/complexity_features/rxn_complexity.csv"

OUTPUT_DIR="data_${DATASET}/75_5_20_same_complex_doc_split"

# Parameters
TEST_VAL_SPLIT=0.25
TEST_SPLIT=0.8
MAX_ITERATIONS=100
TARGET_PVALUE=0.05
RANDOM_SEED=42

# Log
LOG_DIR="./logs"
mkdir -p "$LOG_DIR"

# Run the command
echo "Running adaptive split with KS-test..."
echo "Dataset: $DATASET"
echo "Input file: $INPUT_FILE"
echo "Features file: $FEATURES_FILE"
echo "Output directory: $OUTPUT_DIR"
echo ""

python data_preprocessing/document_split_with_ks_dist.py \
    --input "$INPUT_FILE" \
    --features "$FEATURES_FILE" \
    --output "$OUTPUT_DIR" \
    --test_val_split "$TEST_VAL_SPLIT" \
    --test_split "$TEST_SPLIT" \
    --max_iterations "$MAX_ITERATIONS" \
    --random_seed "$RANDOM_SEED" \
    --log_dir "$LOG_DIR"

echo "Script completed successfully!"

