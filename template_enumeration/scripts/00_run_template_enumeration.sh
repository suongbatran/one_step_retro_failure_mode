#!/bin/bash
# Run template enumeration pipeline for top-inf calculation.
# MUST BE RUN FROM template_enumeration/ directory.
#
# Automatically detects mode based on whether template-relevance predictions exist:
# - Standard mode: If predictions exist at $TEMPL_REL_RESULTS_PATH/predictions.csv,
#                  merges enumeration results with model predictions
# - Standalone mode: If predictions don't exist, runs enumeration only using INPUT_CSV
#
# Environment variables:
#   DATA_NAME           - Dataset name (default: "uspto")
#   TEMPL_PATH          - Path to templates.jsonl (auto-detected if not set)
#   FILTERED_TEMPL_PATH - Path to filtered_templates.jsonl (auto-detected if not set)
#   TEMPL_REL_PATH      - Path to template_relevance directory (default: "../template_relevance")
#   NUM_CORES           - Number of CPU cores to use (default: 8)
#   RAW_CSV_PATH           - CSV file with 'id' and 'rxn_smiles' for standalone mode
#                         (default: $DATA_DIR/raw_test.csv)
#
# Usage:
#   sh scripts/00_run_template_enumeration.sh

set -e  # Exit on error

# Configuration
export DATA_NAME="${DATA_NAME:-uspto}"
export TEMPL_REL_PATH="${TEMPL_REL_PATH:-../template_relevance}"
export TEMPL_REL_PROCESSED_DATA_PATH="$TEMPL_REL_PATH/data/$DATA_NAME/processed"
export TEMPL_REL_RESULTS_PATH="$TEMPL_REL_PATH/results/$DATA_NAME"
export TEMPL_PATH="${TEMPL_PATH:-$TEMPL_REL_PROCESSED_DATA_PATH/templates.jsonl}"
export FILTERED_TEMPL_PATH="${FILTERED_TEMPL_PATH:-$TEMPL_REL_PROCESSED_DATA_PATH/filtered_templates.jsonl}"


export LOG_DIR="./logs/$DATA_NAME"
export DATA_DIR="./data/$DATA_NAME"
export RESULTS_DIR="./results/$DATA_NAME"

# csv file to default to if prediction file does not exist
# requires id and rxn_smiles columns
export RAW_CSV_PATH="$DATA_DIR/raw_test.csv"
export NUM_CORES="${NUM_CORES:-8}"

mkdir -p "$LOG_DIR" "$DATA_DIR" "$RESULTS_DIR"

# Check if template-relevance predictions exist
PREDICTIONS_FILE="$TEMPL_REL_RESULTS_PATH/predictions.csv"

if [ -f "$PREDICTIONS_FILE" ]; then
    # ==================== STANDARD MODE ====================
    echo "=============================================="
    echo "Template Enumeration - Standard Mode"
    echo "(Template-relevance predictions found)"
    echo "=============================================="
    echo "Data name: $DATA_NAME"
    echo "Predictions: $PREDICTIONS_FILE"
    echo "Template file: $TEMPL_PATH"
    echo "Filtered template file: $FILTERED_TEMPL_PATH"
    echo "Output directory: $RESULTS_DIR"
    echo "Number of cores: $NUM_CORES"
    echo "=============================================="

    sh scripts/01_run_templ_rel_scorer.sh
    sh scripts/02_run_enumerator_default.sh
    sh scripts/03_run_enumerator_filtered.sh
    sh scripts/04_merge_enumeration_results.sh

else
    # ==================== STANDALONE MODE ====================
    export STANDALONE=1
    export INPUT_CSV="${INPUT_CSV:-$RAW_CSV_PATH}"

    echo "=============================================="
    echo "Template Enumeration - Standalone Mode"
    echo "(No template-relevance predictions found)"
    echo "=============================================="
    echo "Data name: $DATA_NAME"
    echo "Input CSV: $INPUT_CSV"
    echo "Template file: $TEMPL_PATH"
    echo "Output directory: $RESULTS_DIR"
    echo "=============================================="

    sh scripts/02_run_enumerator_default.sh
    sh scripts/03_run_enumerator_filtered.sh
    sh scripts/04_merge_enumeration_results.sh
fi

echo ""
echo "=============================================="
echo "Template enumeration completed!"
echo "Final results: $RESULTS_DIR/template_enumeration_final_results.csv"
echo "=============================================="

# Analysis scripts that can be run independently
# (uncomment to run)
sh scripts/run_enumerator_unique_precs.sh

# Note: second pass enumeration requires template-relevance predictions 
# sh scripts/run_enumerator_second_pass.sh

