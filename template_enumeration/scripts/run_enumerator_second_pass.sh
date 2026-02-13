#!/bin/bash
# Analysis: template enumeration second pass for candidates 1..5.
# Run from template_enumeration/; expects 00 env (DATA_DIR, RESULTS_DIR, LOG_DIR, TEMPL_REL_RESULTS_PATH, etc.).

# Defaults when not run from 00
if [ -z "$DATA_NAME" ]; then export DATA_NAME="uspto"; fi
if [ -z "$TEMPL_REL_RESULTS_PATH" ]; then
    export TEMPL_REL_PATH="${TEMPL_REL_PATH:-../template_relevance}"
    export TEMPL_REL_RESULTS_PATH="$TEMPL_REL_PATH/results/$DATA_NAME"
fi
if [ -z "$DATA_DIR" ]; then export DATA_DIR="./data/$DATA_NAME"; fi
if [ -z "$RESULTS_DIR" ]; then export RESULTS_DIR="./results/$DATA_NAME"; fi
if [ -z "$LOG_DIR" ]; then export LOG_DIR="./logs/$DATA_NAME"; fi

mkdir -p "$DATA_DIR/second_run" "$RESULTS_DIR/second_run"

LOG_SUBDIR="$LOG_DIR/second_pass"
mkdir -p "$LOG_SUBDIR"


for CAND_NUM in 1 2 3 4 5; do

    INPUT_FILE="$DATA_DIR/second_run/raw_cand_${CAND_NUM}_test.csv"
    INPUT_BASENAME=$(basename "$INPUT_FILE" .csv)
    OUTPUT_FILE="$RESULTS_DIR/second_run/${INPUT_BASENAME}_enumerated.csv"
    TIMESTAMP=$(date +"%y%m%d-%H%Mh")
    LOG_FILE="$LOG_SUBDIR/enumerator_${DATA_NAME}_cand${CAND_NUM}_${TIMESTAMP}.log"

    if [ ! -f "$INPUT_FILE" ]; then
        echo "Second-pass data not found, running get_2nd_pass_data.py for candidate $CAND_NUM..."
        python get_2nd_pass_data.py \
            --templ_rel_results_path "$TEMPL_REL_RESULTS_PATH" \
            --data_dir "$DATA_DIR" \
            --candidate_num "$CAND_NUM" || exit 1
        echo "Done generating second-pass data"
    fi

    echo "Running template enumeration second pass, candidate $CAND_NUM..."
    python templ_enumerator.py \
        --template_file "$TEMPL_PATH" \
        --input_file "$INPUT_FILE" \
        --output_file "$OUTPUT_FILE" \
        --num_cores $NUM_CORES \
        --log_file "$LOG_FILE"
    echo "Template enumeration completed for candidate $CAND_NUM"
done
