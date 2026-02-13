import os
import csv
import json
import logging
import argparse
import multiprocessing
import numpy as np
from rdkit import Chem
from tqdm import tqdm


global G_predictions

"""
Adapted from https://gitlab.com/mlpds_mit/askcosv2/retro/template_relevance/-/blob/main/templ_rel_scorer.py
to save the accuracy of template-relevance predictions to a CSV file. This file can replace the original templ_rel_scorer.py file
so that the accuracy can be saved when computing the initial top-k accuracies without needing to re-run accuracy calculation.
"""

def canonicalize_smiles(smiles: str, remove_atom_number: bool = True):
    """Adapted from Molecular Transformer"""
    smiles = "".join(smiles.split())
    cano_smiles = ""

    mol = Chem.MolFromSmiles(smiles)

    if mol is not None:
        if remove_atom_number:
            [a.ClearProp('molAtomMapNumber') for a in mol.GetAtoms()]

        cano_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        # Sometimes stereochem takes another canonicalization... (just in case)
        mol = Chem.MolFromSmiles(cano_smiles)
        if mol is not None:
            cano_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)

    return cano_smiles


def csv2kv(_args):
    prediction_row, n_best = _args
    k = canonicalize_smiles(prediction_row["prod_smi"])
    v = []

    for i in range(n_best):
        try:
            prediction = prediction_row[f"cand_precursor_{i + 1}"]
        except KeyError:
            break

        if not prediction or prediction == "9999":          # padding
            break

        prediction = canonicalize_smiles(prediction)
        v.append(prediction)

    return k, v


def match_results(_args):
    global G_predictions
    test_line, n_best = _args
    predictions = G_predictions

    accuracy = np.zeros(n_best, dtype=np.float32)
    test_data = json.loads(test_line.strip())
    rxn_smiles = test_data["rxn_smiles"]
    gt, reagent, prod = rxn_smiles.strip().split(">")
    k = canonicalize_smiles(prod)

    if k not in predictions:
        logging.info(f"Product {prod} not found in predictions (after canonicalization), skipping")
        return accuracy, test_data

    gt = canonicalize_smiles(gt)
    for j, prediction in enumerate(predictions[k]):
        if prediction == gt:
            accuracy[j:] = 1.0
            break
    
    if "prod_smi" not in test_data:
        test_data["prod_smi"] = k
    if "precursor" not in test_data:
        test_data["precursor"] = gt

    return accuracy, test_data


def score_main(args):
    """
        Adapted to save the accuracy of template-relevance predictions to a CSV file.
    """
    global G_predictions
    n_best = args.topk

    # Load predictions and transform into a huge table {cano_prod: [cano_cand, ...]}
    args.prediction_file = os.path.join(args.test_output_path, "predictions.csv")
    logging.info(f"Loading predictions from {args.prediction_file}")
    predictions = {}
    p = multiprocessing.Pool(args.num_cores)

    with open(args.prediction_file, "r") as prediction_csv:
        prediction_reader = csv.DictReader(prediction_csv)
        for result in tqdm(p.imap(csv2kv,
                                  ((prediction_row, n_best) for prediction_row in prediction_reader)),
                                  desc="Processing predictions"):
            k, v = result
            predictions[k] = v

    G_predictions = predictions

    p.close()
    p.join()
    p = multiprocessing.Pool(args.num_cores)        # re-initialize to see the global variable

    # Results matching
    test_rxns_with_template_file = os.path.join(
        args.processed_data_path, "test_rxns_with_template.jsonl"
    )
    logging.info(f"Matching against ground truth from {test_rxns_with_template_file}")
    with open(test_rxns_with_template_file, "r") as f:
        lines = f.readlines()

    # Collect all results first
    logging.info("Processing predictions...")
    all_results = list(tqdm(p.imap(match_results,
                    ((line, n_best) for line in lines)),
                    total=len(lines)))
    
    if not all_results:
        logging.error("No results were generated. Check if predictions match the test data.")
        return

    # Write results to CSV
    output_file = os.path.join(args.test_output_path, "predictions_with_accuracy.csv")
    fieldnames = ["id", "rxn_smiles", "prod_smi", "precursor"] + [f"cand_precursor_{i+1}" for i in range(n_best)]
    
    with open(output_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        
        for accuracy, test_data in all_results:
            row = {
                "id": test_data.get("id", ""),
                "rxn_smiles": test_data.get("rxn_smiles", ""),
                "prod_smi": test_data.get("prod_smi", ""),
                "precursor": test_data.get("precursor", "")
            }
            
            # Add candidate precursors with their accuracies
            for i in range(n_best):
                row[f"cand_precursor_{i+1}"] = accuracy[i]
            
            writer.writerow(row)

    p.close()
    p.join()

    # Log statistics
    accuracies = np.stack([acc for acc, _ in all_results])
    mean_accuracies = np.mean(accuracies, axis=0)
    for n in range(n_best):
        logging.info(f"Top {n+1} accuracy: {mean_accuracies[n]}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Score predictions with template-relevance model.")
    parser.add_argument("--processed_data_path", type=str, required=True, help="Path to processed data.")
    parser.add_argument("--test_output_path", type=str, required=True, help="Path to test output.")
    parser.add_argument("--topk", type=int, default=50, help="Number of top predictions to consider.")
    parser.add_argument("--num_cores", type=int, default=8, help="Number of cores to use for multiprocessing.")
    args = parser.parse_args()
    score_main(args)