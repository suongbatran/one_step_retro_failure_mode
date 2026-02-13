import pandas as pd
import sys
import os
import argparse
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from accuracy_metrics.accuracy_utils import canonicalize_smiles

def main():
    parser = argparse.ArgumentParser(description="Assign predictions to test reactions by product SMILES.")
    parser.add_argument("--test_file", type=str, required=True, help="Path to the test CSV file")
    parser.add_argument("--prediction_file", type=str, required=True, help="Path to the predictions CSV file")
    args = parser.parse_args()

    test = pd.read_csv(args.test_file)
    test["prod_smi"] = test["rxn_smiles"].apply(lambda x: canonicalize_smiles(x.split(">")[2]))

    prediction = pd.read_csv(args.prediction_file)
    prediction_unique = prediction.drop_duplicates(subset=["prod_smi"])
    assign_prediction = test[["id", "rxn_smiles", "prod_smi"]].merge(prediction_unique, on="prod_smi", how="left")
    output_path = os.path.dirname(args.prediction_file)
    assign_prediction.to_csv(os.path.join(output_path, "assigned_predictions.csv"), index=False)

if __name__ == "__main__":
    main()