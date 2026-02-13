# Template Enumeration

This module computes **top-∞ accuracy** metrics for template-based retrosynthesis models by exhaustively applying all templates to products and checking whether the ground truth precursors can be recovered. This analysis helps quantify upper-bound model performance and identify failure modes related to template coverage.

**Note:** This module requires a different conda environment than the main repository. Create the environment using:

```bash
conda env create -f environment.yml
conda activate one_step_templ
```

---

## Overview

Template enumeration evaluates three match types:
- **Exact match**: Predicted precursors exactly match the ground truth
- **Exact set match**: Predicted precursor set matches ground truth
- **Superset match**: Ground truth is a subset of predicted precursors

The pipeline processes reactions that the template-relevance model failed to predict correctly at top-k, then applies templates exhaustively to determine if the correct precursors are theoretically recoverable.

---

## Directory Structure

```
template_enumeration/
├── scripts/                          # Shell scripts for running the pipeline
│   ├── 00_run_template_enumeration.sh   # Main entry point
│   ├── 01_run_templ_rel_scorer.sh       # Step 1: Score predictions
│   ├── 02_run_enumerator_default.sh     # Step 2: Default template set
│   ├── 03_run_enumerator_filtered.sh    # Step 3: Filtered template set
│   ├── 04_merge_enumeration_results.sh  # Step 4: Merge results
│   ├── run_enumerator_second_pass.sh    # Second-pass enumeration
│   └── run_enumerator_unique_precs.sh   # Count unique precursors
├── utils/                            # Utility modules
│   ├── chem_utils.py                    # SMILES canonicalization
│   ├── file_utils.py                    # File I/O utilities
│   └── templ_enum_utils.py              # Template application utilities
├── templ_enumerator.py               # Main template enumeration class
├── templ_rel_scorer.py               # Template-relevance prediction scorer
├── templ_enum_uniq_precs.py          # Unique precursor counter
├── extract_templates_and_rows.py     # Extract templates and rows for enumeration
├── get_2nd_pass_data.py              # Generate second-pass input data
└── merge_templ_enum_results.py       # Merge enumeration results
```

---

## Usage

### Prerequisites

**Template files:** Running the template relevance preprocessing will automatically generate the `templates.jsonl` file. Alternatively, the processed template files (`templates.jsonl` and `filtered_templates.jsonl`) can be downloaded from [this link](https://doi.org/10.5281/zenodo.18623634). 

**Standard mode:**
1. Template-relevance model predictions must exist (from the `template_relevance` module)
2. Processed template files (default and filtered) must be available

**Standalone mode:**
1. A CSV file with `id` and `rxn_smiles` columns
2. Processed template files (default and optionally filtered) 

### Running the Full Pipeline

Navigate to the `template_enumeration/` directory and run:

```bash
sh scripts/00_run_template_enumeration.sh
```

The script **automatically detects** which mode to use:
- **Standard mode**: If template-relevance predictions exist at `$TEMPL_REL_RESULTS_PATH/predictions.csv`, merges enumeration results with model predictions
- **Standalone mode**: If no predictions exist, runs enumeration only using `INPUT_CSV` (defaults to `$DATA_DIR/raw_test.csv`)

**Environment variables:**

| Variable | Default | Description |
|----------|---------|-------------|
| `DATA_NAME` | `uspto` | Dataset name |
| `TEMPL_PATH` | Auto-detected | Path to templates.jsonl |
| `FILTERED_TEMPL_PATH` | Auto-detected | Path to filtered_templates.jsonl |
| `TEMPL_REL_PATH` | `../template_relevance` | Path to template_relevance directory |
| `NUM_CORES` | `8` | Number of CPU cores |
| `INPUT_CSV` | `$DATA_DIR/raw_test.csv` | Input CSV for standalone mode |

**Examples:**

```bash
# Standard mode (if predictions exist) or standalone mode (if not)
sh scripts/00_run_template_enumeration.sh

# Standalone mode with custom input
export INPUT_CSV="data/my_reactions.csv"
sh scripts/00_run_template_enumeration.sh

# Custom dataset name
export DATA_NAME="my_dataset"
sh scripts/00_run_template_enumeration.sh
```

**Input CSV format (standalone mode):**

```csv
id,rxn_smiles
1,CC(=O)O.CCO>>CCOC(C)=O
2,c1ccccc1Br.CC>>c1ccccc1C
```

**Pipeline steps:**

Standard mode:
1. **Score predictions**: Adds accuracy metrics to template-relevance predictions
2. **Default enumeration**: Applies templates with min_freq ≥ 5
3. **Filtered enumeration**: Applies remaining low-frequency templates
4. **Merge results**: Combines with model predictions into final summary

Standalone mode:
1. **Default enumeration**: Applies all default templates to input reactions
2. **Filtered enumeration** (if available): Applies low-frequency templates to failures
3. **Merge results**: Combines enumeration results into final summary

### Running Individual Steps

Each step can be run independently after configuring environment variables:

```bash
export DATA_NAME="uspto"
export TEMPL_REL_PATH="../template_relevance"
export TEMPL_REL_PROCESSED_DATA_PATH="$TEMPL_REL_PATH/data/$DATA_NAME/processed"
export TEMPL_REL_RESULTS_PATH="$TEMPL_REL_PATH/results/$DATA_NAME"
export TEMPL_PATH="$TEMPL_REL_PROCESSED_DATA_PATH/templates.jsonl"
export FILTERED_TEMPL_PATH="$TEMPL_REL_PROCESSED_DATA_PATH/filtered_templates.jsonl"
export DATA_DIR="./data/$DATA_NAME"
export RESULTS_DIR="./results/$DATA_NAME"
export LOG_DIR="./logs/$DATA_NAME"
export NUM_CORES=8

# For standalone mode (no predictions), also set:
# export STANDALONE=1
# export INPUT_CSV="data/my_reactions.csv"

sh scripts/02_run_enumerator_default.sh
sh scripts/03_run_enumerator_filtered.sh
sh scripts/04_merge_enumeration_results.sh
```

---

## Core Scripts

### `templ_enumerator.py`

Main template enumeration engine. Applies all templates to products and checks for matches.

```bash
python templ_enumerator.py \
    --template_file <path_to_templates.jsonl> \
    --input_file <path_to_input.csv> \
    --output_file <path_to_output.csv> \
    --num_cores 8 \
    --log_file <path_to_log>
```

**Input CSV columns**: `id`, `product`, `reactants` (or `rxn_smiles`)

**Output columns**: Original columns + `exact_match`, `exact_set_match`, `superset_match`, `exact_template`, etc.

### `templ_rel_scorer.py`

Scores template-relevance model predictions against ground truth.

```bash
python templ_rel_scorer.py \
    --processed_data_path <path_to_processed_data> \
    --test_output_path <path_to_predictions> \
    --topk 50 \
    --num_cores 8
```

### `extract_templates_and_rows.py`

Extracts rows and templates for different enumeration stages:

```bash
# Extract rows where template-relevance model fails (for default-set enumeration)
python extract_templates_and_rows.py extract-rows-default-set \
    --templ_rel_results_path <path> \
    --data_name uspto

# Extract rows where default enumeration fails (for filtered-set enumeration)
python extract_templates_and_rows.py extract-rows-filtered-set \
    --input_path <default_results.csv> \
    --output_path <output.csv>

# Extract low-frequency templates
python extract_templates_and_rows.py extract-filtered-templates \
    --processed_data_path <path> \
    --data_name uspto \
    --min_freq 5
```

### `merge_templ_enum_results.py`

Merges enumeration results back into template-relevance predictions, or runs in standalone mode.

**Standard mode** (merge with template-relevance predictions):
```bash
python merge_templ_enum_results.py \
    --templ_rel_results_path <path_to_predictions> \
    --results_dir <path_to_enumeration_results> \
    --output_file <final_summary.csv> \
    --include_filtered
```

**Standalone mode** (enumeration results only):
```bash
python merge_templ_enum_results.py \
    --default_enum_file <path_to_default_enumerated.csv> \
    --output_file <final_summary.csv> \
    --include_filtered \
    --standalone
```

### `templ_enum_uniq_precs.py`

Counts unique precursors generated by template application (for analyzing template coverage).

```bash
python templ_enum_uniq_precs.py \
    --template_path <templates.jsonl> \
    --filtered_template_path <filtered_templates.jsonl> \
    --input_file <input.csv> \
    --sampled_products_file <sampled.csv> \
    --output_file <output.csv> \
    --num_samples 1000
```

---

## Output

The final output (`template_enumeration_final_results.csv`) contains:

| Column | Description |
|--------|-------------|
| `id` | Reaction identifier |
| `rxn_smiles` | Reaction SMILES |
| `prod_smi` | Product SMILES |
| `precursor` | Ground truth precursor |
| `exact_match` | True if exact match found |
| `exact_set_match` | True if exact set match found |
| `superset_match` | True if superset match found |


## Notes

- Template enumeration is computationally expensive. The pipeline parallelizes across rows (first pass) and across templates (second pass for timed-out rows).
- Rows that time out in the first pass (600s limit per row) are retried with per-template parallelization.
- The filtered template set contains templates with frequency < 5 that were excluded from the default template-relevance model.
