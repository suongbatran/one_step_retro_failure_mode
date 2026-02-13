#!/bin/bash

set -e

DATASET="uspto"

# 1. Assign predictions to test reactions

for model in TR AT G2S; do
    python accuracy_metrics/assign_prediction_to_test.py \
        --test_file data_${DATASET}/75_5_20_same_complex_doc_split/raw_test.csv \
        --prediction_file data_${DATASET}/75_5_20_same_complex_doc_split_results/${model}/predictions.csv
done

# 2. Calculate exact match, superset, stereochemistry-agnostic accuracy

for model in TR AT G2S; do
    python accuracy_metrics/cal_exact_superset_stereo.py \
        --predictions data_${DATASET}/75_5_20_same_complex_doc_split_results/${model}/assigned_predictions.csv \
        --accuracy_type exact_match superset stereochemistry_agnostic \
        --num_processes 32 \
        --n_best 50
done

# 3. Calculate synthon accuracy

## Generate mapped reactions from predicted precursors

for model in TR AT G2S; do
    python accuracy_metrics/map_prediction.py \
        --test_set_w_rxnmapper data_${DATASET}/75_5_20_same_complex_doc_split/raw_test.csv \
        --predictions data_${DATASET}/75_5_20_same_complex_doc_split_results/${model}/assigned_predictions.csv \
        --n_best 50 \
        --num_processes 32
done

## then calculate synthon accuracy
for model in TR AT G2S; do
    python accuracy_metrics/cal_synthon.py \
        --superset_accuracy_file data_${DATASET}/75_5_20_same_complex_doc_split_results/${model}/accuracy_superset.csv \
        --predictions_with_mapping data_${DATASET}/75_5_20_same_complex_doc_split_results/${model}/predictions_with_mapping.csv \
        --num_processes 32 \
        --n_best 50
done

# 4. Calculate two-step superset accuracy

# assigning 2nd run prediction only need to do for AT, G2S, because TR is evaluated at top-inf
for model in AT G2S; do
    python accuracy_metrics/assign_prediction_to_test_2nd_run.py \
        --first_run_assigned_predictions_file data_${DATASET}/75_5_20_same_complex_doc_split_results/${model}/assigned_predictions.csv \
        --second_run_predictions_path data_${DATASET}/75_5_20_same_complex_doc_split_results/${model}/second_run \
        --n_best_second_run 5
    done

# then calculate two-step superset accuracy

python accuracy_metrics/cal_two_step_superset.py \
    --superset_accuracy_file data_${DATASET}/75_5_20_same_complex_doc_split_results/AT/accuracy_superset.csv \
    --second_run_assigned_predictions_path data_${DATASET}/75_5_20_same_complex_doc_split_results/AT/second_run \
    --n_best_two_step 5
    
python accuracy_metrics/cal_two_step_superset_TR.py \
    --superset_accuracy_file data_${DATASET}/75_5_20_same_complex_doc_split_results/TR/accuracy_superset.csv \
    --second_run_assigned_predictions_path data_${DATASET}/75_5_20_same_complex_doc_split_results/TR/second_run \
    --n_best_two_step 5 \
    --num_processes 32