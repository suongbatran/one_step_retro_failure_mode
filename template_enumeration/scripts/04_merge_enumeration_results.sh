#!/bin/bash
# Step 4: Merge template enumeration results into final summary.
#
# Two modes of operation:
# 1. Standard mode: Merges enumeration results with template-relevance predictions
# 2. Standalone mode (STANDALONE=1): Merges enumeration results only

# Paths
OUTPUT_FILE="$RESULTS_DIR/template_enumeration_final_results.csv"
LOG_SUBDIR="$LOG_DIR/merge"
mkdir -p "$LOG_SUBDIR" "$RESULTS_DIR"
TIMESTAMP=$(date +"%y%m%d-%H%Mh")
LOG_FILE="$LOG_SUBDIR/merge_${DATA_NAME}_${TIMESTAMP}.log"

if [ "$STANDALONE" = "1" ]; then
    echo "Running merge in STANDALONE mode"
    python merge_templ_enum_results.py \
        --results_dir "$RESULTS_DIR" \
        --output_file "$OUTPUT_FILE" \
        --include_filtered \
        --standalone \
        --log_file "$LOG_FILE" || exit 1
else
    # Standard mode: merge with template-relevance predictions
    python merge_templ_enum_results.py \
        --templ_rel_results_path "$TEMPL_REL_RESULTS_PATH" \
        --results_dir "$RESULTS_DIR" \
        --output_file "$OUTPUT_FILE" \
        --include_filtered \
        --log_file "$LOG_FILE" || exit 1
fi

echo "Template enumeration results merged into final summary"
echo "Output saved to: $OUTPUT_FILE"
