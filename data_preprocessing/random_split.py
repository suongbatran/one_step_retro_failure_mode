import pandas as pd
from sklearn.model_selection import train_test_split
import json
import os
import argparse


def random_split(input_file, output_dir, test_val_split, test_split):
    
    if input_file.endswith('.csv'):
        df = pd.read_csv(input_file)
    elif input_file.endswith('.json'):
        with open(input_file, "r") as f:
            data = [json.loads(line) for line in f]
            df = pd.DataFrame(data)
    else:
        raise ValueError(f"Unsupported input file type: {input_file}")

    # Train vs (Val + Test)
    df_train, df_test_val = train_test_split(
        df, test_size=test_val_split, random_state=42, shuffle=True
    )

    # Val vs Test (from the test_val set)
    df_val, df_test = train_test_split(
        df_test_val, test_size=test_split, random_state=42, shuffle=True
    )

    print(f"Split sizes - Train: {len(df_train)}, Val: {len(df_val)}, Test: {len(df_test)}")

    df_train = df_train[["id", "rxn_smiles"]]
    df_val = df_val[["id", "rxn_smiles"]]
    df_test = df_test[["id", "rxn_smiles"]]

    os.makedirs(output_dir, exist_ok=True)
    df_train.to_csv(os.path.join(output_dir, "raw_train.csv"), index=False)
    df_val.to_csv(os.path.join(output_dir, "raw_val.csv"), index=False)
    df_test.to_csv(os.path.join(output_dir, "raw_test.csv"), index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform random reaction-based splitting.")
    parser.add_argument('--input', type=str, required=True, help='Path to input JSON or CSVfile with reactions')
    parser.add_argument('--output', type=str, required=True, help='Output directory to save split CSV files')
    parser.add_argument('--test_val_split', type=float, required=True, help='Fraction of reactions for validation + test')
    parser.add_argument('--test_split', type=float, required=True, help='Fraction of validation + test set for test set')
    
    args = parser.parse_args()
    
    random_split(args.input, args.output, args.test_val_split, args.test_split)

# 75% train, 5% val, 20% test
# then use --test_val_split 0.25 --test_split 0.8
