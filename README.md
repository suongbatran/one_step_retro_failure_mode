# Quantifying the Failure Modes of Current One-step Retrosynthesis Models

This repository contains code and scripts for analyzing and quantifying failure modes of one-step retrosynthesis models, as described in our paper (link TBD). Please refer to the paper for detailed methodology and results.

---

## 1. Environment setup

Create the conda environment for the analysis:

```bash
conda env create -f environment.yml
conda activate one_step_retro_failure_mode
```

---

## 2. Data

The USPTO-Full dataset used in this project is available [here](https://zenodo.org/records/18623634). We cannot publicly release CAS or Pistachio datasets.

**Split USPTO-Full data:**

```bash
wget https://zenodo.org/records/18623634/files/data_uspto.zip?download=1 -O data_uspto.zip
unzip data_uspto.zip
```

The `data_uspto/` directory contains:
- Clean & processed USPTO data (with original and RXNMapper-remapped reaction SMILES)
- `complexity_features/`: Computed complexity features
- `75_5_20_same_complex_doc_split/`: Document-based split data was used to train and evaluate model
- `75_5_20_random_split/`: Random split data was used to train and evaluate model
- `75_5_20_same_complex_doc_split_results/`:
    - Model evaluation and complementary metrics for further inspection
    - Top-k accuracy summaries

The `models_uspto/` directory contains:
- Model checkpoints and vocabulary files for the Augmented Transformer (AT) and Graph2SMILES (G2S) models
- Model checkpoints and template files for the Template-relevance (TR) model

---

## 3. Retrosynthesis models (ASKCOS)

- One-step models evaluated in this study are from ASKCOS: for installation, training, or testing, see the [ASKCOS repository](https://gitlab.com/mlpds_mit/askcosv2/retro).
- Navigate `template_enumeration/` folder for tasks related to template enumeration.

---

## 4. Usage

### Data preprocessing 

To create your own split of the dataset, follow these steps:
1. **Download the USPTO dataset**

  Download the USPTO dataset `uspto.reactions.json.tgz` file from [this link](https://chemrxiv.org/engage/chemrxiv/article-details/60c741240f50dbe270395a6d).

  ```bash
  mkdir -p data_uspto
  wget https://chemrxiv.org/doi/suppl/10.26434/chemrxiv.7949024.v1/suppl_file/uspto.reactions.json.tgz -O data_uspto/uspto.reactions.json.tgz
  tar -xzvf data_uspto/uspto.reactions.json.tgz -C data_uspto
  ```

2. **Clean and preprocess reaction data**

   Run the provided preprocessing script to clean the raw reaction data

   ```bash
   sh scripts/preprocess_data.sh
   ```

   This will generate a cleaned reaction file:
   - Output: `data_uspto/uspto_preprocessed.json`

Note: Any code related to NextMove's NameRXN tools has been removed, so your processed files may differ if you do not have access to this software. The details of our preprocessing workflow are provided in the SI link TBD.

3. **Dataset splitting (random or document-based)**
To create splits of 75% train, 5% validation, and 20% test:
  - for document-based split, run:

   ```bash
   sh scripts/run_doc_split_with_ks_dist.sh
   ```

  - for random split, run:

   ```bash
   python data_preprocessing/random_split.py --input data_uspto/uspto_preprocessed.json --output data/75_5_20_random_split --test_val_split 0.25 --test_split 0.8
   ```

   - Outputs:
     - Document-based split: `data_uspto/75_5_20_same_complex_doc_split`
     - Random split: `data_uspto/75_5_20_random_split`
### Calculate complexity metrics

To compute producte and reaction complexity features:

```bash
sh complexity_features/cal_complexity_features_scripts.sh
```
### Calculate accuracy metrics

To compute all accuracies (exact match, stereochemistry-ignore, superset, two-step superset) of the model, first run:

```bash
sh scripts/accuracy_cal_script.sh
```

and summarize the complementary accuracy metrics for each top-k value, run:

```bash
python accuracy_metrics/cal_top_k_accuracy.py --data uspto
```

## 5. Example

If you have a test set and predictions, you can compute all complementary metrics.
We provide an example with a random set of 15 USPTO reactions, and results.
- `examples/test.csv`
- `examples/results/predictions.csv`
- `examples/results/second_run/cand_k_predictions.csv` for k from 1 to 5

Then run:
```bash
sh examples/example_analysis.sh
```

Use `examples/example_view_prediction.ipynb` to view per-entry predictions/accuracy.

---

## Citation

@article{xxx,
  title={Quantifying the Failure Modes of Current One-step Retrosynthesis Models},
  author={...},
  journal={...},
  year={...},
  note={link TBD}
}

