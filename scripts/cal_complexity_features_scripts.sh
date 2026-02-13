#!/bin/bash
set -e

DATASET="uspto"
python complexity_features/cal_product_complexity.py \
    --input data_${DATASET}/${DATASET}_preprocessed.json \
    --output data_${DATASET}/complexity_features/product_complexity.csv

python complexity_features/cal_rxn_complexity.py \
    --input data_${DATASET}/${DATASET}_preprocessed.json \
    --output data_${DATASET}/complexity_features/rxn_complexity.csv

python complexity_features/cal_stereochem.py \
    --input data_${DATASET}/${DATASET}_preprocessed.json \
    --output data_${DATASET}/complexity_features/stereochem.csv
