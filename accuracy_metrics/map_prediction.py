import argparse
import pandas as pd
import numpy as np
import math
import sys
import os
from tqdm import tqdm
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from accuracy_metrics.accuracy_utils import *
from accuracy_metrics.rxnmapper.rxnmapper import BatchedMapper
from multiprocessing import Pool


def map_predicted_rxn(rxn_smiles_rxnmapper, predicted_precursor_list, n_best, device=None):
    prod = canonicalize_smiles(rxn_smiles_rxnmapper.split(">")[2])
    rxn_mapper = BatchedMapper(batch_size=16, device=device)
    to_map_rxn = []
    for k in range(n_best):
        pred = predicted_precursor_list[k]
        pred = canonicalize_smiles(pred)
        if pred is None:
            continue
        to_map_rxn.append(pred + ">>" + prod)
    mapped_prediction_rxn = list(rxn_mapper.map_reactions(to_map_rxn))
    return mapped_prediction_rxn

def map_predicted_rxn_row(args):
    row_id, rxn_smiles, rxn_smiles_rxnmapper, predicted_precursor_list, n_best, device = args

    result = map_predicted_rxn(rxn_smiles_rxnmapper, predicted_precursor_list, n_best, device)
    return [row_id, rxn_smiles, rxn_smiles_rxnmapper] + result

def main():
    parser = argparse.ArgumentParser(description="Map predicted reactions.")
    parser.add_argument("--test_set_w_rxnmapper",type=str,required=True, help="Path to CSV file with test set including id, rxn_smiles and rxn_smiles_rxnmapper column.")
    parser.add_argument("--predictions",type=str,required=True, help="Path to CSV file with predictions.")
    parser.add_argument("--n_best",type=int,default=50, help="Number of top predictions to consider.")
    parser.add_argument("--num_processes",type=int,default=32, help="Number of processes to use for multiprocessing.")
    parser.add_argument("--chunk_size",type=int,default=1000, help="Chunk size for processing.")
    parser.add_argument("--device",type=str,default=None, help="Device to use for rxnmapper (e.g., 'cuda', 'cuda:0', 'cuda:1', 'cpu'). If not specified, auto-detects.")
    args = parser.parse_args()
    predictions = pd.read_csv(args.predictions)
    test_set_w_rxnmapper = pd.read_csv(args.test_set_w_rxnmapper)
    n_best = args.n_best
    num_processes = args.num_processes
    output_path = os.path.dirname(args.predictions)
    out_file = os.path.join(output_path, f"predictions_with_mapping.csv")
    
    predictions = predictions.merge(test_set_w_rxnmapper, on=["id", "rxn_smiles"], how="left")

    # To avoid explosion of GPU memory, process a chunk of 1000 rows at a time, then clear the memory manually.
    cols = ["id", "rxn_smiles", "rxn_smiles_rxnmapper"] + [f"cand_precursor_{i+1}" for i in range(n_best)]
    chunk_size = args.chunk_size
    total_rows = len(predictions)
    num_chunks = math.ceil(total_rows / chunk_size)
    
    def clear_gpu_memory():
        try:
            import torch
            torch.cuda.empty_cache()
        except ImportError:
            pass
        try:
            import gc
            gc.collect()
        except Exception:
            pass

    with open(out_file, "w") as f:
        header = ",".join(["id", "rxn_smiles", "rxn_smiles_rxnmapper"] + [f"mapped_predicted_rxn_{i+1}" for i in range(n_best)])
        f.write(header + "\n")
        for chunk_idx in range(num_chunks):
            chunk_start = chunk_idx * chunk_size
            chunk_end = min((chunk_idx + 1) * chunk_size, total_rows)
            input_chunk = [
                (row["id"], row["rxn_smiles"], row["rxn_smiles_rxnmapper"], [row[f"cand_precursor_{i+1}"] for i in range(n_best)], n_best, args.device)
                for _, row in predictions[cols].iloc[chunk_start:chunk_end].iterrows()
            ]
            with Pool(num_processes) as pool:
                for res in tqdm(pool.imap(map_predicted_rxn_row, input_chunk), total=len(input_chunk), desc=f"Chunk {chunk_idx+1}/{num_chunks}"):
                    out_line = ",".join(map(str, res))
                    f.write(out_line + "\n")
            clear_gpu_memory()
        print(f"Saved to {out_file}")

if __name__ == "__main__":
    main()
