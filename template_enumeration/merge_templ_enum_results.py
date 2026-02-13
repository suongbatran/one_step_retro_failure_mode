#!/usr/bin/env python3
"""
Merge template enumeration results back into template-relevance predictions.

Two modes of operation:
1. Standard mode: Merge with template-relevance predictions
   - Template-relevance predictions with accuracy (predictions_with_accuracy.csv)
   - Default-set template enumeration results (default_set_reactions_enumerated.csv)
   - Optionally: filtered-set enumeration results (filtered_set_reactions_enumerated.csv)

2. Standalone mode (--standalone): Use enumeration results only
   - Default-set template enumeration results (input_enumerated.csv or default_set_reactions_enumerated.csv)
   - Optionally: filtered-set enumeration results (filtered_set_reactions_enumerated.csv)

Output is a single summary CSV with exact_match, exact_set_match, and superset_match.
"""

import argparse
import logging
import os
import pandas as pd

from utils.file_utils import setup_logging


def load_standalone_input(input_path: str) -> pd.DataFrame:
    """
    Load enumeration results for standalone mode.
    
    The enumeration results already contain exact_match, exact_set_match, superset_match.
    """
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Enumeration results not found: {input_path}")
    
    df = pd.read_csv(input_path, on_bad_lines="warn")
    
    # Ensure required columns exist
    required_cols = ["id", "rxn_smiles", "exact_match", "exact_set_match", "superset_match"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in enumeration results: {missing}")
    
    # Standardize column names (enumeration output may have 'product'/'reactants' instead of 'prod_smi'/'precursor')
    if "prod_smi" not in df.columns and "product" in df.columns:
        df["prod_smi"] = df["product"]
    if "precursor" not in df.columns and "reactants" in df.columns:
        df["precursor"] = df["reactants"]
    
    return df


def load_predictions_with_accuracy(templ_rel_results_path: str) -> pd.DataFrame:
    """
    Load template-relevance predictions with accuracy and derive exact_match from top-k.

    Uses the last cand_precursor_* column (e.g., cand_precursor_50) as the model
    exact_match flag (1 = correct within top-k).
    """
    path = os.path.join(templ_rel_results_path, "predictions_with_accuracy.csv")
    if not os.path.exists(path):
        raise FileNotFoundError(f"Predictions with accuracy not found under {templ_rel_results_path}")

    df = pd.read_csv(path)
    cand_cols = [c for c in df.columns if c.startswith("cand_precursor_")]
    if not cand_cols:
        raise ValueError("No cand_precursor_* columns found in predictions file")
    last_cand_col = sorted(cand_cols, key=lambda x: int(x.split("_")[-1]))[-1]
    df["exact_match"] = df[last_cand_col].astype(bool)

    cols_to_keep = ["id", "rxn_smiles", "prod_smi", "precursor", "exact_match"]
    return df[cols_to_keep]


def load_enumeration_results(path: str) -> pd.DataFrame:
    """Load a template enumeration output CSV (default or filtered set)."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Enumeration results not found: {path}")
    return pd.read_csv(path, on_bad_lines="warn")


def _merge_enumeration_into_model(
    model_df: pd.DataFrame,
    enum_df: pd.DataFrame,
    suffix: str = "_templ",
    merge_cols: list = None,
) -> pd.DataFrame:
    """
    Merge enumeration match columns into the model dataframe on rxn_smiles.

    For rows where the model already has exact_match=True, keep True.
    Otherwise use the enumeration result. Fills NaN with False for match cols.
    """
    if merge_cols is None:
        merge_cols = [
            "exact_match",
            "exact_set_match",
            "superset_match",
        ]
    rename_map = {c: c + suffix for c in merge_cols if c in enum_df.columns}
    enum_subset = enum_df.rename(columns=rename_map)
    cols_to_merge = ["rxn_smiles"] + list(rename_map.values())
    merged = model_df.merge(
        enum_subset[cols_to_merge],
        on="rxn_smiles",
        how="left",
    )

    bool_cols = [c for c in rename_map.values() if "match" in c]
    for c in bool_cols:
        if c in merged.columns:
            merged[c] = merged[c].fillna(False)

    # Combined: model correct => True; else use templ
    if "exact_match" + suffix in merged.columns:
        merged["exact_match"] = merged["exact_match"] | merged["exact_match" + suffix]
    if "exact_set_match" + suffix in merged.columns:
        merged["exact_set_match"] = merged.apply(
            lambda row: row["exact_match"] or row["exact_set_match" + suffix],
            axis=1,
        )
    if "superset_match" + suffix in merged.columns:
        merged["superset_match"] = merged.apply(
            lambda row: row["exact_match"] or row["superset_match" + suffix],
            axis=1,
        )

    return merged


def merge_default_set_enumeration(
    model_df: pd.DataFrame,
    default_enum_path: str,
) -> pd.DataFrame:
    """
    Merge default-set template enumeration results into model predictions.

    Drops temporary columns and keeps final id, rxn_smiles, prod_smi, precursor,
    exact_match, exact_set_match, superset_match.
    """
    enum_df = load_enumeration_results(default_enum_path)
    merged = _merge_enumeration_into_model(
        model_df,
        enum_df,
        suffix="_templ",
        merge_cols=["exact_match", "exact_set_match", "superset_match"],
    )
    final_cols = [
        "id", "rxn_smiles", "prod_smi", "precursor",
        "exact_match", "exact_set_match", "superset_match",
    ]
    return merged[[c for c in final_cols if c in merged.columns]]


def merge_filtered_set_enumeration(
    merged_df: pd.DataFrame,
    filtered_enum_path: str,
) -> pd.DataFrame:
    """
    Merge filtered-set (top-inf all) enumeration results into an already-merged dataframe.

    Updates exact_match, exact_set_match, superset_match from filtered results where
    the model was still wrong.
    """
    enum_df = load_enumeration_results(filtered_enum_path)
    suffix = "_templ_all"
    rename_map = {
        "exact_match": "exact_match" + suffix,
        "exact_set_match": "exact_set_match" + suffix,
        "superset_match": "superset_match" + suffix,
    }
    enum_subset = enum_df.rename(columns=rename_map)
    cols_to_merge = ["rxn_smiles"] + list(rename_map.values())
    merged = merged_df.merge(
        enum_subset[cols_to_merge],
        on="rxn_smiles",
        how="left",
    )

    for c in ["exact_match" + suffix, "exact_set_match" + suffix, "superset_match" + suffix]:
        if c in merged.columns:
            merged[c] = merged[c].fillna(False)

    merged["exact_match"] = merged["exact_match"] | merged["exact_match" + suffix]
    merged["exact_set_match"] = merged.apply(
        lambda row: row["exact_match"] or row["exact_set_match" + suffix],
        axis=1,
    )
    merged["superset_match"] = merged.apply(
        lambda row: row["exact_match"] or row["superset_match" + suffix],
        axis=1,
    )

    final_cols = [
        "id", "rxn_smiles", "prod_smi", "precursor",
        "exact_match", "exact_set_match", "superset_match",
    ]
    return merged[[c for c in final_cols if c in merged.columns]]


def merge_standalone(
    default_enum_path: str,
    filtered_enum_path: str = None,
    include_filtered: bool = False,
) -> pd.DataFrame:
    """
    Merge enumeration results in standalone mode (no template-relevance predictions).
    
    Args:
        default_enum_path: Path to default-set enumeration results
        filtered_enum_path: Path to filtered-set enumeration results (optional)
        include_filtered: Whether to include filtered-set results
        
    Returns:
        DataFrame with combined enumeration results
    """
    logging.info("Running in standalone mode (no template-relevance predictions)")
    
    # Load default enumeration as base
    logging.info("Loading default-set enumeration from %s", default_enum_path)
    base_df = load_standalone_input(default_enum_path)
    logging.info("Loaded %d rows from default enumeration", len(base_df))
    print_summary(base_df, "Default set only")
    
    # Merge filtered results if requested
    if include_filtered and filtered_enum_path and os.path.exists(filtered_enum_path):
        logging.info("Merging filtered-set enumeration from %s", filtered_enum_path)
        base_df = merge_filtered_set_enumeration(base_df, filtered_enum_path)
        logging.info("After merging filtered-set enumeration:")
        print_summary(base_df, "With filtered set")
    elif include_filtered:
        logging.warning("Filtered enumeration file not found, skipping: %s", filtered_enum_path)
    
    return base_df


def check_rxn_smiles_consistency(model_df: pd.DataFrame, enum_df: pd.DataFrame) -> None:
    """Log a warning if the same id has different rxn_smiles between model and enumeration."""
    common_ids = set(model_df["id"]).intersection(enum_df["id"])
    merged = model_df[["id", "rxn_smiles"]].merge(
        enum_df[["id", "rxn_smiles"]],
        on="id",
        suffixes=("_model", "_templ"),
    )
    mismatches = merged[merged["rxn_smiles_model"] != merged["rxn_smiles_templ"]]
    if len(mismatches) > 0:
        logging.warning("Mismatched rxn_smiles for same id: %d rows", len(mismatches))
    else:
        logging.info("All rxn_smiles match for common ids (n=%d)", len(common_ids))


def print_summary(df: pd.DataFrame, label: str = "") -> None:
    """Print match counts and percentages."""
    n = len(df)
    if n == 0:
        logging.warning("Empty dataframe, no summary")
        return
    exact = df["exact_match"].sum()
    exact_set = df["exact_set_match"].sum()
    superset = df["superset_match"].sum()
    logging.info("%sRows: %d", label + " " if label else "", n)
    logging.info("Exact matches: %d (%.2f%%)", exact, 100.0 * exact / n)
    logging.info("Exact set matches: %d (%.2f%%)", exact_set, 100.0 * exact_set / n)
    logging.info("Superset matches: %d (%.2f%%)", superset, 100.0 * superset / n)


def main():
    parser = argparse.ArgumentParser(
        description="Merge template enumeration results into template-relevance predictions",
    )
    parser.add_argument(
        "--templ_rel_results_path",
        type=str,
        default=None,
        help="Path to template-relevance results dir (contains predictions_with_accuracy.csv). "
             "Not required in standalone mode.",
    )
    parser.add_argument(
        "--results_dir",
        type=str,
        default=None,
        help="Path to template enumeration results dir (contains default_set_reactions_enumerated.csv)",
    )
    parser.add_argument(
        "--default_enum_file",
        type=str,
        default=None,
        help="Path to default-set enumeration results CSV. If not provided, uses "
             "<results_dir>/default_set_reactions_enumerated.csv",
    )
    parser.add_argument(
        "--filtered_enum_file",
        type=str,
        default=None,
        help="Path to filtered-set enumeration results CSV. If not provided, uses "
             "<results_dir>/filtered_set_reactions_enumerated.csv",
    )
    parser.add_argument(
        "--output_file",
        type=str,
        required=True,
        help="Output summary CSV path",
    )
    parser.add_argument(
        "--include_filtered",
        action="store_true",
        help="Also merge filtered-set enumeration for top-inf all",
    )
    parser.add_argument(
        "--standalone",
        action="store_true",
        help="Run in standalone mode without template-relevance predictions. "
             "Uses enumeration results directly.",
    )
    parser.add_argument(
        "--log_file",
        type=str,
        default=None,
        help="Log file path (default: stdout)",
    )
    args = parser.parse_args()

    setup_logging(log_file=args.log_file)

    # Determine paths for enumeration results
    if args.default_enum_file:
        default_enum_path = args.default_enum_file
    elif args.results_dir:
        default_enum_path = os.path.join(
            args.results_dir,
            "default_set_reactions_enumerated.csv",
        )
    else:
        parser.error("Either --default_enum_file or --results_dir must be provided")
    
    if args.filtered_enum_file:
        filtered_enum_path = args.filtered_enum_file
    elif args.results_dir:
        filtered_enum_path = os.path.join(
            args.results_dir,
            "filtered_set_reactions_enumerated.csv",
        )
    else:
        filtered_enum_path = None

    if args.standalone:
        # Standalone mode: use enumeration results directly
        merged = merge_standalone(
            default_enum_path,
            filtered_enum_path,
            include_filtered=args.include_filtered,
        )
    else:
        # Standard mode: merge with template-relevance predictions
        if not args.templ_rel_results_path:
            parser.error("--templ_rel_results_path is required when not using --standalone")
        
        logging.info("Loading predictions with accuracy from %s", args.templ_rel_results_path)
        model_df = load_predictions_with_accuracy(args.templ_rel_results_path)
        logging.info("Total predictions: %d", len(model_df))
        logging.info("Model exact matches (top-k): %d", model_df["exact_match"].sum())

        logging.info("Loading default-set enumeration from %s", default_enum_path)
        default_enum_df = load_enumeration_results(default_enum_path)
        check_rxn_smiles_consistency(model_df, default_enum_df)
        logging.info(
            "Default-set enum: exact=%d, exact_set=%d, superset=%d",
            default_enum_df["exact_match"].sum(),
            default_enum_df["exact_set_match"].sum(),
            default_enum_df["superset_match"].sum(),
        )

        merged = merge_default_set_enumeration(model_df, default_enum_path)
        logging.info("After merging default-set enumeration:")
        print_summary(merged, "Default set only")

        if args.include_filtered and filtered_enum_path and os.path.exists(filtered_enum_path):
            logging.info("Merging filtered-set enumeration from %s", filtered_enum_path)
            merged = merge_filtered_set_enumeration(merged, filtered_enum_path)
            logging.info("After merging filtered-set enumeration (top-inf all):")
            print_summary(merged, "With filtered set")
        elif args.include_filtered:
            logging.warning("Filtered enumeration file not found, skipping: %s", filtered_enum_path)

    os.makedirs(os.path.dirname(args.output_file) or ".", exist_ok=True)
    merged.to_csv(args.output_file, index=False)
    logging.info("Saved summary to %s", args.output_file)


if __name__ == "__main__":
    main()
