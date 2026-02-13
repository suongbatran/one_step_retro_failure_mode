set -e

## Assign predictions to test reactions
python accuracy_metrics/assign_prediction_to_test.py \
    --test_file examples/test.csv \
    --prediction_file examples/results/predictions.csv

## Calculate exact, superset, and stereochemistry-agnostic accuracy
python accuracy_metrics/cal_exact_superset_stereo.py \
    --predictions examples/results/assigned_predictions.csv \
    --accuracy_type exact_match superset stereochemistry_agnostic \
    --num_processes 8 \
    --n_best 50


## Calculate synthon accuracy

#generate mapped reaction of predictions
python accuracy_metrics/map_prediction.py \
    --test_set_w_rxnmapper examples/test.csv \
    --predictions examples/results/assigned_predictions.csv \
    --n_best 50 \
    --num_processes 8

#calculate
python accuracy_metrics/cal_synthon.py \
    --superset_accuracy_file examples/results/accuracy_superset.csv \
    --predictions_with_mapping examples/results/predictions_with_mapping.csv \
    --num_processes 8 \
    --n_best 50

## Calculate two-step superset accuracy

#assign 2nd run prediction
python accuracy_metrics/assign_prediction_to_test_2nd_run.py \
    --first_run_assigned_predictions_file examples/results/assigned_predictions.csv \
    --second_run_predictions_path examples/results/second_run \
    --n_best_second_run 5

#calculate
python accuracy_metrics/cal_two_step_superset.py \
    --superset_accuracy_file examples/results/accuracy_superset.csv \
    --second_run_assigned_predictions_path examples/results/second_run \
    --n_best_two_step 5


### Complexity calculation
python complexity_features/cal_product_complexity.py \
    --input examples/test.csv \
    --output examples/complexity_features/product_complexity.csv

python complexity_features/cal_rxn_complexity.py \
    --input examples/test.csv \
    --output examples/complexity_features/rxn_complexity.csv

