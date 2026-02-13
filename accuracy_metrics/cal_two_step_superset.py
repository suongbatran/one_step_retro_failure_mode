import argparse
import pandas as pd
import os
import sys
from tqdm import tqdm

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from accuracy_metrics.accuracy_utils import two_step_superset


def main():
    parser = argparse.ArgumentParser(description="Calculate two-step superset accuracy.")
    parser.add_argument("--superset_accuracy_file", type=str, required=True)
    parser.add_argument("--second_run_assigned_predictions_path", type=str, required=True)
    parser.add_argument("--n_best_two_step", type=int, default=5)
    args = parser.parse_args()

    superset_accuracy = pd.read_csv(args.superset_accuracy_file)
    n_best = args.n_best_two_step

    second_pass_dict = {}
    for i in range(1, n_best + 1):
        df = pd.read_csv(
            os.path.join(args.second_run_assigned_predictions_path,
                         f"cand_{i}_assigned_predictions.csv")
        ).set_index("id") # make sure id index is unique
        second_pass_dict[i - 1] = df

    # output
    acc_cols = [f"cand_precursor_{i+1}" for i in range(n_best)]
    out_path = os.path.join(os.path.dirname(args.superset_accuracy_file),
                            "accuracy_two_step_superset.csv")

    with open(out_path, "w") as f_out:
        f_out.write(",".join(["id", "rxn_smiles"] + acc_cols) + "\n")


    for _, row in tqdm(superset_accuracy.iterrows(),
                       total=superset_accuracy.shape[0],
                       desc="Calculating two-step superset accuracy"):

        row_id = row["id"]

        row_candidates = {
            i: second_pass_dict[i].loc[row_id]
            for i in range(n_best)
        }

        acc = two_step_superset(
            [row[f"cand_precursor_{i+1}"] for i in range(n_best)],
            row_candidates,
            n_best_two_step=n_best
        )

        with open(out_path, "a") as f_out:
            f_out.write(",".join([str(row["id"]), row["rxn_smiles"]] +
                                 [str(a) for a in acc]) + "\n")

    print(f"Saved to {out_path}")


if __name__ == "__main__":
    main()
