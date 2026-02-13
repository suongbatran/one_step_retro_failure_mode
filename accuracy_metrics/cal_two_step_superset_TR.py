'''
Calculate two-step superset accuracy for TR model.
For second run of template-based models, we enumerated all filtered templates.
So just only look at where there is a superset match or not.
'''

import argparse
import pandas as pd
import numpy as np
import sys
import os
from tqdm import tqdm
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def two_step_superset_TR_row(superset_accuracy_list, second_pass_superset_list, n_best_two_step=5):
    accuracy = np.zeros(n_best_two_step, dtype=float)
    for k in range(n_best_two_step):
        if superset_accuracy_list[k] == 1:
            accuracy[k:] = 1.0
            return accuracy.tolist()
        else:
            if second_pass_superset_list[k] == True:
                accuracy[k:] = 1.0
                break
    return accuracy.tolist()

def process_row(args):
    """Process a single row to calculate two-step superset accuracy."""
    row, n_best_two_step = args
    superset_accuracy_list = [row[f"cand_precursor_{i+1}"] for i in range(n_best_two_step)]
    second_pass_superset_list = [row.get(f"cand_{i+1}_superset_match", False) for i in range(n_best_two_step)]
    
    acc = two_step_superset_TR_row(superset_accuracy_list, second_pass_superset_list, n_best_two_step)
    row_data = [row['id'], row['rxn_smiles']] + list(acc)
    row_data_str = [str(x) for x in row_data]
    return ",".join(row_data_str) + "\n"

def main():
    parser = argparse.ArgumentParser(description="Calculate various accuracy metrics for predictions.")
    parser.add_argument("--superset_accuracy_file", type=str, required=True, help="Path to CSV file with calculated superset accuracy.")
    parser.add_argument("--second_run_assigned_predictions_path",type=str,required=True, help="Path to CSV file with second run assigned predictions.")
    parser.add_argument("--n_best_two_step",type=int,default=5, help="Number of top predictions to consider.")
    parser.add_argument("--num_processes",type=int,default=32, help="Number of parallel jobs. Default: 32.")
    args = parser.parse_args()

    superset_accuracy = pd.read_csv(args.superset_accuracy_file)    
    output_path = os.path.dirname(args.second_run_assigned_predictions_path)
    n_best_two_step = args.n_best_two_step
    num_processes = args.num_processes

    # Merge second-pass superset match columns
    for i in range(1, n_best_two_step + 1):
        second_pass_superset_match_i = pd.read_csv(os.path.join(args.second_run_assigned_predictions_path, f"cand_{i}_assigned_predictions.csv"))
        second_pass_superset_match_i = second_pass_superset_match_i[["id", "superset_match"]]
        second_pass_superset_match_i = second_pass_superset_match_i.rename(columns={"superset_match": f"cand_{i}_superset_match"})
        superset_accuracy = superset_accuracy.merge(second_pass_superset_match_i, on="id", how="left")
        
    # Select relevant columns
    two_step_superset_accuracy = superset_accuracy[["id", "rxn_smiles"] + [f"cand_precursor_{i+1}" for i in range(n_best_two_step)] + [f"cand_{i+1}_superset_match" for i in range(n_best_two_step)]].copy()

    acc_cols = [f"cand_precursor_{i+1}" for i in range(n_best_two_step)]
    output_path = os.path.join(output_path, "accuracy_two_step_superset.csv")

    header = ['id', 'rxn_smiles'] + acc_cols
    with open(output_path, 'w') as f_out:
        f_out.write(",".join(header) + "\n")

    # Convert DataFrame to list of dictionaries for multiprocessing
    rows_list = two_step_superset_accuracy.to_dict('records')
    
    # Create tuples of (row, n_best_two_step) for multiprocessing
    rows_with_params = [(row, n_best_two_step) for row in rows_list]
    
    # Process rows in parallel
    with Pool(num_processes) as pool:
        results = list(tqdm(
            pool.imap(process_row, rows_with_params),
            total=len(rows_with_params),
            desc="Calculating and writing two-step superset accuracy"
        ))
    
    # Write all results to file
    with open(output_path, 'a') as f_out:
        f_out.writelines(results)
    
    print(f"Saved to {output_path}")

if __name__ == "__main__":
    main()
