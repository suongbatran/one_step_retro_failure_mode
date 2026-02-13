import pandas as pd
import sys
import os
import argparse
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from accuracy_metrics.accuracy_utils import canonicalize_smiles, sort_fragments_by_heavy_atom

def split_prediction_to_fragments(cand_precursor_smiles):
    try:
        cand_precursor_smiles = canonicalize_smiles(cand_precursor_smiles)
        if cand_precursor_smiles.count(".") > 0:
            cand_precursor_smiles_fragements = cand_precursor_smiles.split(".")
            cand_precursor_smiles_fragements = sort_fragments_by_heavy_atom(cand_precursor_smiles_fragements)
            largest_fragment = cand_precursor_smiles_fragements[0]
            smaller_fragments = ".".join(cand_precursor_smiles_fragements[1:])
            return largest_fragment, smaller_fragments
        else:
            return cand_precursor_smiles, ""
    except:
        return "", ""

def main():
    parser = argparse.ArgumentParser(description="Assign predictions to test reactions by product SMILES.")
    parser.add_argument("--first_run_assigned_predictions_file", type=str, required=True, help="Path to the first run assigned predictions CSV file")
    parser.add_argument("--second_run_predictions_path", type=str, required=True, help="Path to the second run predictions CSV file")
    parser.add_argument("--n_best_second_run", type=int, required=True, help="Number of best second run predictions to consider")
    args = parser.parse_args()

    first_run_assigned_predictions = pd.read_csv(args.first_run_assigned_predictions_file)

    for i in range(1, args.n_best_second_run + 1):
        # split cand_precursor_{i} into largest and smaller fragments and assign as new columns.
        first_run_assigned_predictions[[f"cand_{i}_largest_frag", f"cand_{i}_smaller_frags"]] = first_run_assigned_predictions.apply(
            lambda row: pd.Series(split_prediction_to_fragments(row[f"cand_precursor_{i}"])), axis=1
        )

        second_run_prediction = pd.read_csv(os.path.join(args.second_run_predictions_path, f"cand_{i}_predictions.csv"))
        second_run_prediction = second_run_prediction.drop_duplicates(subset=["prod_smi"])

        second_run_assign_prediction = first_run_assigned_predictions[["id", "rxn_smiles", f"cand_{i}_largest_frag", f"cand_{i}_smaller_frags"]].merge(
            second_run_prediction, left_on=f"cand_{i}_largest_frag", right_on="prod_smi", how="left"
        )
        second_run_assign_prediction.to_csv(os.path.join(args.second_run_predictions_path, f"cand_{i}_assigned_predictions.csv"), index=False)
        print(f"Assigned predictions for cand_{i} to test reactions")

if __name__ == "__main__":
    main()