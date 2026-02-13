#!/bin/bash
# Step 3: Template enumeration with the filtered template set (min_freq).
#
# This script extracts rows where the default enumeration failed (exact_match=False)
# and runs enumeration with the filtered (low-frequency) templates.
#
# Works with both standard and standalone modes - uses the output from step 2.
# In standalone mode, set DEFAULT_ENUM_FILE to override the default results path.

# Config
MIN_FREQ=5
echo "Using FILTERED templates (min frequency $MIN_FREQ)"

# Default set results file (output from step 2)
DEFAULT_SET_RESULTS_FILE="${DEFAULT_ENUM_FILE:-$RESULTS_DIR/default_set_reactions_enumerated.csv}"

# Paths
PROCESSED_ROWS_FILE="$DATA_DIR/filtered_set_reactions.csv"
INPUT_FILE="$PROCESSED_ROWS_FILE"
INPUT_BASENAME=$(basename "$INPUT_FILE" .csv)
OUTPUT_FILE="$RESULTS_DIR/${INPUT_BASENAME}_enumerated.csv"

LOG_SUBDIR="$LOG_DIR/filtered"
mkdir -p "$LOG_SUBDIR" "$RESULTS_DIR" "$DATA_DIR"
TIMESTAMP=$(date +"%y%m%d-%H%Mh")
LOG_FILE="$LOG_SUBDIR/enumerator_${DATA_NAME}_filtered_${TIMESTAMP}.log"


if [ -f "$PROCESSED_ROWS_FILE" ]; then
    echo "Processed rows file exists, skipping extraction"
else
    echo "Processed rows file does not exist, extracting rows"
    if [ ! -f "$DEFAULT_SET_RESULTS_FILE" ]; then
        echo "Error: default set results not found at $DEFAULT_SET_RESULTS_FILE."
        echo "Run step 2 (02_run_enumerator_default.sh) first."
        exit 1
    fi
    python extract_templates_and_rows.py extract-rows-filtered-set \
        --input_path "$DEFAULT_SET_RESULTS_FILE" \
        --output_path "$PROCESSED_ROWS_FILE"
    echo "Processed rows file extracted"
fi

if [ ! -f "$PROCESSED_ROWS_FILE" ]; then
    echo "Error: processed rows file not found at $PROCESSED_ROWS_FILE. Cannot run enumeration."
    exit 1
fi

# Run enumeration
echo "Starting template enumeration for filtered template set"
python templ_enumerator.py \
    --template_file "$FILTERED_TEMPL_PATH" \
    --input_file "$INPUT_FILE" \
    --output_file "$OUTPUT_FILE" \
    --num_cores $NUM_CORES \
    --log_file "$LOG_FILE" || exit 1
echo "Template enumeration completed for filtered set"
echo "Output saved to: $OUTPUT_FILE"
