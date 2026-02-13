import pandas as pd
import argparse
from tqdm import tqdm
import os 
from rdkit import Chem

def sort_fragments_by_heavy_atom(fragments):
    # Sort fragments by heavy atom count
    fragments = sorted(fragments, key=lambda x: len(Chem.MolFromSmiles(x).GetAtoms()), reverse=True)
    return fragments

def main(args):
    templ_rel_results_path = args.templ_rel_results_path
    data_dir = args.data_dir
    candidate_num = args.candidate_num

    input_prediction_path = os.path.join(templ_rel_results_path, "predictions_with_accuracy.csv")
    input_prediction_mol_path = os.path.join(templ_rel_results_path, "predictions.csv")
    output_dir = os.path.join(data_dir, "second_run")

    if not os.path.exists(input_prediction_path):
        raise FileNotFoundError(f"Predictions with accuracy not found: {input_prediction_path}")
    if not os.path.exists(input_prediction_mol_path):
        raise FileNotFoundError(f"Predictions not found: {input_prediction_mol_path}")

    prediction = pd.read_csv(input_prediction_path, low_memory=False)
    prediction_mol = pd.read_csv(input_prediction_mol_path, low_memory=False)
    print("Done reading prediction files, start filtering with num rows:", len(prediction))

    # Filter
    filter_cols = [f"cand_precursor_{i}" for i in range(1, candidate_num + 1)]
    # Filter out rows where all candidate precursors are 0 
    prediction = prediction[
        (prediction[filter_cols].eq(0).all(axis=1)) & 
        (prediction['precursor'].notna()) &
        (prediction['prod_smi'].notna())
    ]

    # New data
    new_data = {
        "id": [],
        "rxn_smiles": [],
        "candidate_small_frag": []
    }

    
    for _, row in tqdm(prediction.iterrows(), total=len(prediction), desc="Processing"):
        id = row["id"]
        prod = row["prod_smi"]
        precursor = row["precursor"]
        row_mol = prediction_mol[prediction_mol["prod_smi"] == prod]
        candidate = row_mol[f"cand_precursor_{candidate_num}"].iloc[0]
        if candidate == 9999 or candidate == "9999":
            continue
        
        small_frag =""

        try:
            
            if candidate.count(".") > 0:
            # Split the candidate into fragments
                fragments = candidate.split(".")
                # Sort fragments by length 
                fragments = sort_fragments_by_heavy_atom(fragments)
                candidate = fragments[0]  # biggest part becomes the main candidate
                # Add smaller parts to the original precursor
                small_frag += ".".join(fragments[1:])

            new_rxn_smiles = precursor + ">>" + candidate
            
        except Exception as e:
            # Error generated when cannot find mol in prediction file (very small number of this case)
            print(f"Error generating for id: {id}, precursor: {precursor}, candidate: {candidate}, error: {e}")
        
        new_data["id"].append(id)
        new_data["rxn_smiles"].append(new_rxn_smiles)
        new_data["candidate_small_frag"].append(small_frag)

    new_df = pd.DataFrame(new_data)

    # Saving
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"raw_cand_{candidate_num}_test.csv")
    new_df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate second-pass input from template-relevance predictions.",
    )
    parser.add_argument(
        "--templ_rel_results_path",
        type=str,
        required=True,
        help="Path to template-relevance results dir (predictions_with_accuracy.csv, predictions.csv)",
    )
    parser.add_argument(
        "--data_dir",
        type=str,
        required=True,
        help="Path to data dir; output written to <data_dir>/second_run/raw_cand_<N>_test.csv",
    )
    parser.add_argument(
        "--candidate_num",
        type=int,
        required=True,
        help="Candidate precursor number (1–5)",
    )
    args = parser.parse_args()
    main(args)