import argparse
import pandas as pd
from tqdm import tqdm
import os
import sys

def get_topk_acc(data, topk):
    col = f'cand_precursor_{topk}'
    if col not in data.columns:
        return None
    acc = (data[col] == 1).mean() * 100
    return acc

def load_and_eval(path, topk):
    if not os.path.exists(path):
        return None
    try:
        df = pd.read_csv(path)
        return get_topk_acc(df, topk)
    except Exception as e:
        print(f"Error loading or processing {path}: {e}", file=sys.stderr)
        return None

def main():
    parser = argparse.ArgumentParser(description="Compute top-k accuracy summary.")
    parser.add_argument('--data', type=str, required=True, help="cas, pistachio, uspto")
    args = parser.parse_args()

    base_data_dir = f"data_{args.data}"

    models = ["AT", "G2S", "TR"]
    topk_list = [int(i) for i in range(1, 51)]
    accuracy_columns = ["exact_match", "superset", "stereochemistry_agnostic", "synthon", "two_step_superset"]
    results = []
    for model in models:
        for acc_type in accuracy_columns:
            for k in tqdm(topk_list, desc=f"  {acc_type}", leave=False):
                path = os.path.join(base_data_dir, f"75_5_20_same_complex_doc_split_results/{model}/accuracy_{acc_type}.csv")
                acc = load_and_eval(path, topk=k)
                results.append({
                    'Model': model,
                    'Accuracy': round(acc, 2) if acc is not None else "NA",
                    'Top_k': k,
                    'Accuracy_type': acc_type
                })

    summary_dir = os.path.join(base_data_dir, "75_5_20_same_complex_doc_split_results/accuracy_summary")
    try:
        os.makedirs(summary_dir, exist_ok=True)
    except Exception as e:
        print(f"Failed to create summary directory {summary_dir}: {e}", file=sys.stderr)
        sys.exit(1)
    df_results = pd.DataFrame(results)
    output_path = os.path.join(summary_dir, "all_topk_accuracy.csv")
    try:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        df_results.to_csv(output_path, index=False)
        print(f"Top-k accuracy summary written to: {output_path}")
    except Exception as e:
        print(f"Failed to write output CSV: {e}", file=sys.stderr)
        sys.exit(2)

if __name__ == "__main__":
    main()

# 
# python accuracy_metrics/cal_top_k_accuracy.py --data pistachio