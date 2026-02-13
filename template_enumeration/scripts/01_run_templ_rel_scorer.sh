#!/bin/bash
# Step 1: Run template-relevance scorer to add accuracy to predictions.
# Reads: $TEMPL_REL_RESULTS_PATH/predictions.csv, ground truth from processed data.
# Writes: $TEMPL_REL_RESULTS_PATH/predictions_with_accuracy.csv

# Config
TOPK=50

# Run scorer (if predictions_with_accuracy.csv does not exist)
if [ -f "$TEMPL_REL_RESULTS_PATH/predictions_with_accuracy.csv" ]; then
    echo "Predictions with accuracy file exists, skipping template-relevance scorer"
else
    echo "Predictions with accuracy file does not exist, running template-relevance scorer"
    python templ_rel_scorer.py \
        --processed_data_path "$TEMPL_REL_PROCESSED_DATA_PATH" \
        --test_output_path "$TEMPL_REL_RESULTS_PATH" \
        --topk "$TOPK" \
        --num_cores "$NUM_CORES" || exit 1
    echo "Template-relevance scorer completed"
fi



