import argparse
import pandas as pd
import numpy as np
import sys
import os
from tqdm import tqdm
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from accuracy_metrics.accuracy_utils import exact_match, superset, stereochemistry_agnostic
from multiprocessing import Pool

def cal_accuracy_row(args):
    """Calculate accuracy for a single row."""
    row_id, rxn_smiles, cand_list, n_best, accuracy_type = args

    accuracy_func_map = {
        "exact_match": exact_match,
        "superset": superset,
        "stereochemistry_agnostic": stereochemistry_agnostic
    }
    result = accuracy_func_map[accuracy_type](rxn_smiles, cand_list, n_best=n_best)
    return [row_id, rxn_smiles] + result

def main():
    parser = argparse.ArgumentParser(description="Calculate various accuracy metrics for predictions.")
    parser.add_argument("--predictions",type=str,required=True, help="Path to CSV file with predictions.")
    parser.add_argument("--accuracy_type", nargs="+", required=True, choices=["exact_match", "superset", "stereochemistry_agnostic"], help="Which accuracy metric to calculate.")
    parser.add_argument("--num_processes",type=int,default=32, help="Number of processes to use for multiprocessing.")
    parser.add_argument("--n_best",type=int,default=50, help="Number of top predictions to consider.")
    args = parser.parse_args()

    first_pass_pred = pd.read_csv(args.predictions)
    n_best = args.n_best
    output_path = os.path.dirname(args.predictions)

    cols = ["id", "rxn_smiles"] + [f"cand_precursor_{i+1}" for i in range(n_best)]

    num_processes = args.num_processes


    for accuracy_type in args.accuracy_type:
        print(f"Calculating {accuracy_type.replace('_', ' ')} accuracy...")
        out_file = os.path.join(output_path, f"accuracy_{accuracy_type}.csv")
        with open(out_file, "w") as f:
            header = ",".join(cols)
            f.write(header + "\n")

            inputs = [
                (row["id"], row["rxn_smiles"], [row[f"cand_precursor_{i+1}"] for i in range(n_best)], n_best, accuracy_type)
                for _, row in first_pass_pred[cols].iterrows()
            ]

            with Pool(num_processes) as pool:
                for res in tqdm(pool.imap(cal_accuracy_row, inputs), total=len(inputs)):
                    out_line = ",".join(map(str, res))
                    f.write(out_line + "\n")
            print(f"Saved to {out_file}")
    


if __name__ == "__main__":
    main()
