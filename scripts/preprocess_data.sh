#!/bin/bash
set -e

DATASET="uspto"

python data_preprocessing/preprocess_one_step.py \
    --input data_${DATASET}/${DATASET}.reactions.json \
    --output data_${DATASET}/${DATASET}_preprocessed.json