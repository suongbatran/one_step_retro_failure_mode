# Turn off RDKit logging to ignore the canonicalization errors (https://github.com/rdkit/rdkit/issues/2683)
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import argparse
from complexity_features.stereochem import *
from tqdm import tqdm
import csv
from multiprocessing import Pool
import json
import pandas as pd 


def process_row(row):

    row['has_p_atom_stereo'], row['has_p_bond_stereo'], row['atom_changes'], row['bond_changes'] = check_stereochemistry_changes(row['rxn_smiles'])
    
    return row

def parse_arguments():
    parser = argparse.ArgumentParser(description='Get stereochemistry feactures from reaction SMILES')
    parser.add_argument('--input', type=str, required=True,
                       help='Input JSON or CSV file containing reaction SMILES')
    parser.add_argument('--output', type=str, required=True,
                       help='Output CSV file to save product complexity features')
    parser.add_argument('--processes', type=int, default=32)
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
        train_data = train_data[["id", "rxn_smiles"]]
        rows = train_data.to_dict(orient="records")
    elif input_file.endswith(".csv"):
        # For CSV, let pandas read the file directly
        train_data = pd.read_csv(input_file)
        # Compute total rows (excluding header)
        total_rows = len(train_data)
        # Try to select the columns if available, else use all
        wanted_columns = [col for col in ["id", "rxn_smiles"] if col in train_data.columns]
        if wanted_columns:
            train_data = train_data[wanted_columns]
        rows = train_data.to_dict(orient="records")
    else:
        raise ValueError("Unsupported input file type: {}".format(input_file))

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    results = []
    first_row_written = False
    with open(output_file, 'w', newline='') as outfile:
        writer = None
        with Pool(processes=num_processes) as pool:
            for updated_row in tqdm(pool.imap(process_row, rows), total=len(rows), desc="Processing"):
                if updated_row is not None:
                    if not first_row_written:
                        writer = csv.DictWriter(outfile, fieldnames=updated_row.keys())
                        writer.writeheader()
                        first_row_written = True
                    writer.writerow(updated_row)
                    results.append(updated_row)

    print(f"Results are saved at {output_file}!")

