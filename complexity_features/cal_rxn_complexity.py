# Turn off RDKit logging to ignore the canonicalization errors (https://github.com/rdkit/rdkit/issues/2683)
import logging
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

import sys
import os
sys.path.append("../")
import argparse
from rxn_complexity import *
from tqdm import tqdm
import csv
from multiprocessing import Pool
import json
import pandas as pd 


def process_row(row):
    try:
        row['num_react_atom_forward'] = count_num_react_atom(row['rxn_smiles'])
        row['num_react_atom_retro'] = count_num_react_atom_retro(row['rxn_smiles'])
        row['num_change_ring'] = count_num_change_ring(row['rxn_smiles'])

        row['num_react_atom_forward_rxnmapper'] = count_num_react_atom(row['rxn_smiles_rxnmapper'])
        row['num_react_atom_retro_rxnmapper'] = count_num_react_atom_retro(row['rxn_smiles_rxnmapper'])
        row['num_change_ring_rxnmapper'] = count_num_change_ring(row['rxn_smiles_rxnmapper'])
        return row

    except Exception:
        print(f"Error processing row {row['id']}: {row['rxn_smiles']}")
        return None

def parse_arguments():
    parser = argparse.ArgumentParser(description='Extract complexity features from reaction SMILES')
    parser.add_argument('--input', type=str, required=True,
                       help='Input JSON or CSV file containing reaction data')
    parser.add_argument('--output', type=str, required=True,
                       help='Output CSV file to save complexity features')
    parser.add_argument('--processes', type=int, default=32,
                       help='Number of processes for multiprocessing (default: 32)')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_arguments()
    
    input_file = args.input
    output_file = args.output
    num_processes = args.processes
    
    # Detect input file type based on extension
    if input_file.endswith(".json"):
        with open(input_file, 'r') as f:
            total_rows = sum(1 for _ in f)

        with open(input_file, "r") as f:
            train_data = [json.loads(line) for line in f]
        train_data = pd.DataFrame(train_data)
        train_data = train_data[["id", "rxn_smiles", "rxn_smiles_rxnmapper"]]
        rows = train_data.to_dict(orient="records")
    elif input_file.endswith(".csv"):
        # For CSV, let pandas read the file directly
        train_data = pd.read_csv(input_file)
        # Compute total rows (excluding header)
        total_rows = len(train_data)
        # Try to select the columns if available, else use all
        wanted_columns = [col for col in ["id", "rxn_smiles", "rxn_smiles_rxnmapper"] if col in train_data.columns]
        if wanted_columns:
            train_data = train_data[wanted_columns]
        rows = train_data.to_dict(orient="records")
    else:
        raise ValueError("Unsupported input file type: {}".format(input_file))

    print("Start getting features!")
    results = []
    with Pool(processes=num_processes) as pool:
        for updated_row in tqdm(pool.imap(process_row, rows), total=len(rows), desc="Processing"):
            if updated_row is not None:
                results.append(updated_row)

    if results:
        with open(output_file, 'w', newline='') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)

    print(f"Results are saved at {output_file}!")
    print(f"Processed {len(results)} out of {len(rows)} rows successfully.")

