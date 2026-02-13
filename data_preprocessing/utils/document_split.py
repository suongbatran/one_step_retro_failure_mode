#########
# Document-based splitting for reaction data
##########
import os
import json
import logging
import argparse
logging.getLogger('rdkit').setLevel(logging.ERROR)

import pandas as pd
from sklearn.model_selection import train_test_split

def document_split(input_file, output_dir, test_val_split, test_split, random_state=42, return_dataframes=False):
    """
    Perform random document-based splitting by splitting canIDs, not individual reactions.
    
    Args:
        input_file: Path to input JSON file with reactions
        output_dir: Directory to save split CSV files
        test_val_split: Fraction of documents for validation + test
        test_split: Fraction of validation + test set for test set
        random_state: Random seed for reproducibility
    """
    
    with open(input_file, "r") as f:
        train_data = [json.loads(line) for line in f]
        train_data = pd.DataFrame(train_data)

    # Get unique document IDs
    doc_ids = train_data['canID'].unique()

    # First split: Train vs (Val+Test)
    train_docs, valtest_docs = train_test_split(
        doc_ids, 
        test_size=test_val_split, 
        random_state=random_state
    )

    # Second split: Validation vs Test
    val_docs, test_docs = train_test_split(
        valtest_docs, 
        test_size=test_split, 
        random_state=random_state
    )

    print(f"Split sizes (docs) - Train: {len(train_docs)}, Val: {len(val_docs)}, Test: {len(test_docs)}")

    # Create final dataframes (all reactions belonging to selected documents)
    df_train = train_data[train_data["canID"].isin(train_docs)][["id", "rxn_smiles"]]
    df_val = train_data[train_data["canID"].isin(val_docs)][["id", "rxn_smiles"]]
    df_test = train_data[train_data["canID"].isin(test_docs)][["id", "rxn_smiles"]]

    print(f"Split sizes (reactions) - Train: {len(df_train)}, Val: {len(df_val)}, Test: {len(df_test)}")

    # Save to CSV
    os.makedirs(output_dir, exist_ok=True)

    df_train.to_csv(os.path.join(output_dir, "raw_train.csv"), index=False)
    df_val.to_csv(os.path.join(output_dir, "raw_val.csv"), index=False)
    df_test.to_csv(os.path.join(output_dir, "raw_test.csv"), index=False)

    if return_dataframes:
        return df_train, df_val, df_test
    else:
        return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform document-based splitting.")
    parser.add_argument('--input', type=str, required=True, help='Path to input JSON file with reactions')
    parser.add_argument('--output', type=str, required=True, help='Output directory to save split CSV files')
    parser.add_argument('--test_val_split', type=float, required=True, help='Fraction of reactions for validation + test')
    parser.add_argument('--test_split', type=float, required=True, help='Fraction of validation + test set for test set')
    
    args = parser.parse_args()
    
    document_split(args.input, args.output, args.test_val_split, args.test_split)