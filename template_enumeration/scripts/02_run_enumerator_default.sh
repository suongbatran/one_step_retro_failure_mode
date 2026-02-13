#!/bin/bash
# Step 2: Template enumeration with the default template set.
#
# Two modes of operation:
# 1. Standard mode: Extracts rows from template-relevance predictions
# 2. Standalone mode (STANDALONE=1): Uses provided input CSV directly
#
# For standalone mode, set:
#   export STANDALONE=1
#   export INPUT_CSV="path/to/input.csv"  # CSV with 'id' and 'rxn_smiles' columns

# Defaults (when not run from 00_run_template_enumeration.sh)
if [ -z "$DATA_NAME" ]; then export DATA_NAME="uspto"; fi
if [ -z "$TEMPL_REL_RESULTS_PATH" ]; then
    export TEMPL_REL_PATH="${TEMPL_REL_PATH:-../template_relevance}"
    export TEMPL_REL_RESULTS_PATH="$TEMPL_REL_PATH/results/$DATA_NAME"
fi
if [ -z "$DATA_DIR" ]; then export DATA_DIR="./data/$DATA_NAME"; fi
if [ -z "$RESULTS_DIR" ]; then export RESULTS_DIR="./results/$DATA_NAME"; fi
if [ -z "$LOG_DIR" ]; then export LOG_DIR="./logs/$DATA_NAME"; fi

# Determine input file based on mode
# Output always goes to default_set_reactions_enumerated.csv for consistency
OUTPUT_FILE="$RESULTS_DIR/default_set_reactions_enumerated.csv"

if [ "$STANDALONE" = "1" ]; then
    # Standalone mode: use provided input CSV
    if [ -z "$INPUT_CSV" ]; then
        echo "Error: STANDALONE=1 but INPUT_CSV not set"
        echo "Set INPUT_CSV to path of CSV file with 'id' and 'rxn_smiles' columns"
        exit 1
    fi
    if [ ! -f "$INPUT_CSV" ]; then
        echo "Error: Input file not found: $INPUT_CSV"
        exit 1
    fi
    INPUT_FILE="$INPUT_CSV"
    echo "Running in STANDALONE mode with input: $INPUT_FILE"
else
    # Standard mode: extract from template-relevance results
    PROCESSED_ROWS_FILE="$DATA_DIR/default_set_reactions.csv"
    INPUT_FILE="$PROCESSED_ROWS_FILE"
fi

LOG_SUBDIR="$LOG_DIR/default"
mkdir -p "$LOG_SUBDIR" "$RESULTS_DIR"
TIMESTAMP=$(date +"%y%m%d-%H%Mh")
LOG_FILE="$LOG_SUBDIR/enumerator_${DATA_NAME}_default_${TIMESTAMP}.log"

# Extract rows (if needed and not in standalone mode)
if [ "$STANDALONE" != "1" ]; then
    if [ -f "$PROCESSED_ROWS_FILE" ]; then
        echo "Processed rows file exists, skipping extraction"
    else
        echo "Processed rows file does not exist, extracting rows"
        python extract_templates_and_rows.py extract-rows-default-set \
            --templ_rel_results_path "$TEMPL_REL_RESULTS_PATH" \
            --data_name "$DATA_NAME" \
            --output_file "$PROCESSED_ROWS_FILE" || exit 1
        echo "Done extracting rows"
    fi
fi

# Run enumeration
echo "Starting template enumeration for default template set"
python templ_enumerator.py \
    --template_file "$TEMPL_PATH" \
    --input_file "$INPUT_FILE" \
    --output_file "$OUTPUT_FILE" \
    --num_cores $NUM_CORES \
    --log_file "$LOG_FILE" || exit 1
echo "Template enumeration completed for default set"
echo "Output saved to: $OUTPUT_FILE"



