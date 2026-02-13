import sys, os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

import json
from utils.clean_rxnsmi_utils import *
from utils.smiles_utils import remove_mapping_rxnsmiles
from tqdm import tqdm
import argparse
import pandas as pd
from multiprocessing import Pool, cpu_count

def process_entry(entry):

    rxn_smiles = clean_rxn_smiles(entry['rxn_smiles'])
    if rxn_smiles is None:
        return None, 'step clean reaction SMILES'

    rxn_smiles_no_map = remove_mapping_rxnsmiles(rxn_smiles)
    return {
        'id': entry['id'],
        'canID': entry['canID'],
        'rxn_smiles': rxn_smiles,
        'rxn_smiles_no_map': rxn_smiles_no_map,
    }, None

def process_parallel(route_data, output_file, n_jobs=None):
    if n_jobs is None:
        n_jobs = cpu_count()

    print(f"Running in parallel using {n_jobs} processes...")
    with Pool(n_jobs) as pool:
        results = list(tqdm(pool.imap(process_entry, route_data), total=len(route_data)))

    seen_rxns = set()
    cleaned_data = []

    for result in results:
        if result is None:
            continue
        rxn_no_map = result['rxn_smiles_no_map']
        if rxn_no_map in seen_rxns:
            continue
        seen_rxns.add(rxn_no_map)
        cleaned_data.append({
            'id': result['id'],
            'canID': result['canID'],
            'rxn_smiles': result['rxn_smiles'],
        })

    with open(output_file, "w") as f:
        for item in cleaned_data:
            json.dump(item, f)
            f.write('\n')
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Preprocess one-step reactions.")
    parser.add_argument('--input', type=str, required=True, help='Path to input json file')
    parser.add_argument('--output', type=str, required=True, help='Path to output preprocessed json file')
    args = parser.parse_args()

    route_file = args.input
    output_one_step_file = args.output

    with open(route_file, "r") as f:
        data_dedup = [json.loads(line) for line in f]

    df = pd.DataFrame(data_dedup)
    print(f"Number of reaction in input file: {df.shape[0]}")

    # Remove duplicate reactions
    df = df.drop_duplicates(subset=['rxn_smiles'])
    # remove rxn not having 2 '>'
    df = df[df['rxn_smiles'].str.count('>') == 2]

    print(f"Number of reaction after remove duplication: {df.shape[0]}")
    data_dedup = df.to_dict(orient='records')

    print(f"Start preprocessing...")

    process_parallel(data_dedup, output_one_step_file, n_jobs=32)