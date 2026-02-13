#!/usr/bin/env python3
"""
Script for template and row extraction operations.

1. Extract rows for default set enumeration - extracts rows where template-relevance model fails
2. Extract rows for filtered set enumeration - extracts rows where default template enumeration fails
3. Extract filtered templates - extracts templates filtered out due to min_freq threshold
"""

import argparse
import csv
import json
import logging
import os
import pandas as pd
from collections import defaultdict
from typing import Dict, List, Tuple

# Functions for extracting rows for default set enumeration

def extract_rows_for_default_set_enumeration(args):
    """
    Extracts rows where the template-relevance model fails to predict the correct top-N precursors.
    Only these rows are used for template enumeration for computational efficiency.
    The results are merged back into the original predictions file for top-inf accuracy (min freq = 5).

    Example:
    python extract_templates_and_rows.py extract-rows-default-set uspto 500

    Output:
    data/uspto/rows_for_default_set_enum.csv
    """
    # input 
    data_name = args.data_name
    
    
    input_predictions_path = f"{args.templ_rel_results_path}/predictions.csv"
    input_prediction_accuracy_path = f"{args.templ_rel_results_path}/predictions_with_accuracy.csv"
    output_path = f"data/{args.data_name}"

    predictions = pd.read_csv(input_prediction_accuracy_path)
    prediction_cols = [
        col for col in predictions.columns.tolist() 
        if col.startswith('cand_precursor_')
    ]

    max_col = sorted(prediction_cols, key=lambda x: int(x.split('_')[-1]))[-1]

    print("Done reading predictions files, start filtering with num rows:", len(predictions))

    # Filter out rows where the maximum candidate precursor (cand_precursor_i for max i) is 0
    predictions = predictions[
        (predictions[max_col].eq(0)) & 
        (predictions['precursor'].notna()) &
        (predictions['prod_smi'].notna())
    ]

    # Filtered data ready for output
    print(f"Filtered data has {len(predictions)} rows")

    predictions = predictions.drop(columns=prediction_cols)

    # Use provided output file or default
    if args.output_file:
        output_file = args.output_file
    else:
        output_file = f"{output_path}/rows_for_default_set_enum.csv"

    # Saving
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    predictions.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")


# Functions for extracting rows for filtered set enumeration

def extract_rows_for_filtered_set_enumeration(args):
    """
    Extracts rows where template enumeration (with the default templates) fails to recover recorded precursors.
    Only these rows are used for template enumeration with the filtered templates for computational efficiency.
    The results are merged back into the original predictions file for top-inf accuracy (all).
    """
    summary_df = pd.read_csv(args.input_path)
    print("Done reading predictions files, start filtering with num rows:", len(summary_df))

    summary_df = summary_df[summary_df["exact_match"] == False]

    # Filtered data ready for output
    print(f"Filtered data has {len(summary_df)} rows")

    cols_to_keep = ["id", "rxn_smiles", "product", "reactants"]
    summary_df = summary_df[cols_to_keep]

    # Saving
    os.makedirs(os.path.dirname(args.output_path), exist_ok=True)
    summary_df.to_csv(args.output_path, index=False)
    print(f"Results saved to {args.output_path}")


# Functions for extracting filtered templates

def count_templates_from_csv(reaction_file: str) -> Tuple[Dict[str, int], Dict[str, Dict[str, bool]]]:
    """
    Count the frequency of each template from the CSV file and collect metadata.
    
    Args:
        reaction_file: Path to reactions_with_templates.csv file
        
    Returns:
        Tuple of (template_counts, template_metadata)
    """
    templates = defaultdict(int)
    template_metadata = defaultdict(lambda: {"intra_only": False, "dimer_only": False})
    
    with open(reaction_file, "r") as f:
        csv_reader = csv.DictReader(f)
        
        for row in csv_reader:
            # Check if template extraction was successful
            if not row.get("success", "").lower() == "true":
                continue
                
            canon_reaction_smarts = row.get("canon_reaction_smarts", "")
            if not canon_reaction_smarts:
                continue
                
            templates[canon_reaction_smarts] += 1
            
            # Update metadata - if any reaction has these flags as True, keep them True
            # CSV stores these as "True"/"False" strings, so we need to parse them
            if row.get("intra_only", "").lower() == "true":
                template_metadata[canon_reaction_smarts]["intra_only"] = True
            if row.get("dimer_only", "").lower() == "true":
                template_metadata[canon_reaction_smarts]["dimer_only"] = True
    
    return dict(templates), dict(template_metadata)


def get_filtered_templates(templates: Dict[str, int], min_freq: int) -> Dict[str, int]:
    """
    Get templates that have frequency < min_freq (filtered out templates).
    
    Args:
        templates: Dictionary mapping template SMILES to frequency count
        min_freq: Minimum frequency threshold
        
    Returns:
        Dictionary of filtered out templates
    """
    filtered_templates = {
        template: count 
        for template, count in templates.items() 
        if count < min_freq
    }
    return filtered_templates


def save_filtered_templates(
    filtered_templates: Dict[str, int], 
    template_metadata: Dict[str, Dict[str, bool]],
    output_file: str, 
    data_name: str,
    min_freq: int
) -> None:
    """
    Save filtered templates to a JSONL file with metadata.
    
    Args:
        filtered_templates: Dictionary of filtered templates
        template_metadata: Dictionary of template metadata
        output_file: Output file path
        data_name: Name of the dataset
        min_freq: Minimum frequency threshold that was used
    """
    id_prefix = data_name.replace(" ", "")
    
    with open(output_file, "w") as f:
        for i, (canon_reaction_smarts, count) in enumerate(filtered_templates.items()):
            # Get the actual metadata for this template
            metadata = template_metadata.get(canon_reaction_smarts, {})
            
            template_obj = {
                "index": i,
                "reaction_smarts": canon_reaction_smarts,
                "count": count,
                "necessary_reagent": "",
                "intra_only": metadata.get("intra_only", False),
                "dimer_only": metadata.get("dimer_only", False),
                "template_set": data_name,
                "references": [],
                "attributes": {
                    "ring_delta": 1.0,
                    "chiral_delta": 0
                },
                "_id": f"{id_prefix}_filtered_{i}",
                "filtered": True,
                "original_min_freq": min_freq
            }
            f.write(f"{json.dumps(template_obj)}\n")


def extract_filtered_templates(args):
    """
    Extract templates that were filtered out due to min_freq threshold.
    """
    log_file = f"logs/extract_filtered_templates_{args.data_name}.log"
    
    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    
    # Construct file paths
    reaction_file = os.path.join(
        args.processed_data_path,
        f"reactions_with_templates.{args.data_name}.csv"
    )
    
    if args.output_file is None:
        args.output_file = os.path.join(
            args.processed_data_path,
            "filtered_templates.jsonl"
        )
    
    # Check if input file exists
    if not os.path.exists(reaction_file):
        logging.error(f"Reaction file not found: {reaction_file}")
        return
    
    logging.info(f"Reading reactions from: {reaction_file}")
    logging.info(f"Using min_freq threshold: {args.min_freq}")
    
    # Count all templates and collect metadata
    templates, template_metadata = count_templates_from_csv(reaction_file)
    logging.info(f"Total unique templates found: {len(templates)}")
    
    # Get filtered templates
    filtered_templates = get_filtered_templates(templates, args.min_freq)
    logging.info(f"Filtered templates (freq < {args.min_freq}): {len(filtered_templates)}")
    
    # Print some statistics
    if filtered_templates:
        min_count = min(filtered_templates.values())
        max_count = max(filtered_templates.values())
        avg_count = sum(filtered_templates.values()) / len(filtered_templates)
        logging.info(f"Filtered template frequency range: {min_count} - {max_count}")
        logging.info(f"Average frequency of filtered templates: {avg_count:.2f}")
        
        # Count templates with special flags
        intra_count = sum(1 for template in filtered_templates.keys() 
                         if template_metadata.get(template, {}).get("intra_only", False))
        dimer_count = sum(1 for template in filtered_templates.keys() 
                         if template_metadata.get(template, {}).get("dimer_only", False))
        logging.info(f"Filtered templates with intra_only=True: {intra_count}")
        logging.info(f"Filtered templates with dimer_only=True: {dimer_count}")
    
    # Save filtered templates
    save_filtered_templates(filtered_templates, template_metadata, args.output_file, args.data_name, args.min_freq)
    logging.info(f"Saved filtered templates to: {args.output_file}")



def main():
    parser = argparse.ArgumentParser(
        description="Unified script for template and row extraction operations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", required=True, help="Command to run")

    # --- extract-rows-default-set ---
    p_default = subparsers.add_parser(
        "extract-rows-default-set",
        help="Extract rows where template-relevance model fails (for default-set enumeration)",
    )
    p_default.add_argument(
        "--templ_rel_results_path",
        type=str,
        required=True,
        help="Path to template-relevance results directory (contains predictions_with_accuracy.csv)",
    )
    p_default.add_argument(
        "--data_name",
        type=str,
        required=True,
        help="Dataset name (e.g. uspto), used for default output path",
    )
    p_default.add_argument(
        "--output_file",
        type=str,
        default=None,
        help="Output CSV path (default: data/{data_name}/rows_for_default_set_enum.csv)",
    )
    p_default.set_defaults(func=extract_rows_for_default_set_enumeration)

    # --- extract-rows-filtered-set ---
    p_filtered = subparsers.add_parser(
        "extract-rows-filtered-set",
        help="Extract rows where default template enumeration failed (for filtered-set enumeration)",
    )
    p_filtered.add_argument("--input_path", type=str, required=True, help="Input summary CSV path")
    p_filtered.add_argument("--output_path", type=str, required=True, help="Output CSV path")
    p_filtered.set_defaults(func=extract_rows_for_filtered_set_enumeration)

    # --- extract-filtered-templates ---
    p_templates = subparsers.add_parser(
        "extract-filtered-templates",
        help="Extract templates filtered out by min_freq threshold",
    )
    p_templates.add_argument(
        "--processed_data_path",
        type=str,
        required=True,
        help="Path to processed data directory",
    )
    p_templates.add_argument("--data_name", type=str, required=True, help="Dataset name")
    p_templates.add_argument(
        "--min_freq",
        type=int,
        default=5,
        help="Minimum frequency threshold used for filtering (default: 5)",
    )
    p_templates.add_argument(
        "--output_file",
        type=str,
        default=None,
        help="Output JSONL path (default: filtered_templates.jsonl in processed_data_path)",
    )
    p_templates.set_defaults(func=extract_filtered_templates)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()

