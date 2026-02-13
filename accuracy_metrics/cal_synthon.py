import argparse
import pandas as pd
import numpy as np
import sys
import os
from tqdm import tqdm
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from accuracy_metrics.accuracy_utils import synthon
from multiprocessing import Pool


def cal_synthon_accuracy_row(args):
    """Calculate accuracy for a single row."""
    row_id, rxn_smiles, rxn_smiles_rxnmapper, superset_accuracy_list, predicted_mapped_rxn_list, n_best = args

    result = synthon(rxn_smiles_rxnmapper, superset_accuracy_list, predicted_mapped_rxn_list, n_best)
    return [row_id, rxn_smiles, rxn_smiles_rxnmapper] + result


def main():
    parser = argparse.ArgumentParser(description="Calculate various accuracy metrics for predictions.")
    parser.add_argument("--superset_accuracy_file", type=str, required=True, help="Path to CSV file with calculated superset accuracy.")
    parser.add_argument("--predictions_with_mapping",type=str,required=True, help="Path to CSV file with predictions with mapping.")
    parser.add_argument("--num_processes",type=int,default=32, help="Number of processes to use for multiprocessing.")
    parser.add_argument("--n_best",type=int,default=50, help="Number of top predictions to consider.")
    args = parser.parse_args()

    print("CAUTION: Synthon accuracy need atom-mapping of predicted reactions!")

    predictions_with_mapping = pd.read_csv(args.predictions_with_mapping)
    superset_accuracy = pd.read_csv(args.superset_accuracy_file)    
    output_path = os.path.dirname(args.predictions_with_mapping)

    num_processes = args.num_processes

    # Join superset_accuracy columns as a list for each row
    n_best = args.n_best
    superset_acc_cols = [f"cand_precursor_{i+1}" for i in range(n_best)]
    superset_accuracy["superset_accuracy"] = superset_accuracy[superset_acc_cols].values.tolist()
    superset_accuracy = superset_accuracy[["id", "rxn_smiles", "superset_accuracy"]]
    predictions_with_mapping = predictions_with_mapping.merge(superset_accuracy, on=["id", "rxn_smiles"], how="left")


    print(f"Calculating synthon accuracy...")
    out_file = os.path.join(output_path, f"accuracy_synthon.csv")
    with open(out_file, "w") as f:
        headers = ["id", "rxn_smiles", "rxn_smiles_rxnmapper"] + [f"cand_precursor_{i+1}" for i in range(n_best)]
        header = ",".join(headers)
        f.write(header + "\n")

        cols = ["id", "rxn_smiles", "rxn_smiles_rxnmapper", "superset_accuracy"] + [f"mapped_predicted_rxn_{i+1}" for i in range(n_best)]
        inputs = [
            (row["id"], row["rxn_smiles"], row["rxn_smiles_rxnmapper"], row["superset_accuracy"], [row[f"mapped_predicted_rxn_{i+1}"] for i in range(n_best)], n_best)
            for _, row in predictions_with_mapping[cols].iterrows()
        ]

        with Pool(num_processes) as pool:
            for res in tqdm(pool.imap(cal_synthon_accuracy_row, inputs), total=len(inputs)):
                out_line = ",".join(map(str, res))
                f.write(out_line + "\n")
        print(f"Saved to {out_file}")



if __name__ == "__main__":
    main()
