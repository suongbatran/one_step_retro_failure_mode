
import json
import logging
import pandas as pd
import concurrent.futures
from tqdm import tqdm
from pathlib import Path

from rdkit import Chem, RDLogger
from rdchiral.main import rdchiralReaction, rdchiralReactants, rdchiralRun
from typing import List, Tuple, Optional, Set, Dict, Any, Union

from pebble import ProcessPool, ThreadPool
from concurrent.futures import TimeoutError, ThreadPoolExecutor
from utils.chem_utils import canonicalize_smiles
from utils.file_utils import load_templates_as_list, setup_logging
from utils.templ_enum_utils import has_passed_filter, apply_template_to_product

# Disable RDKit warning and error messages
RDLogger.DisableLog('rdApp.*')
default_result = {
    'exact_template': "",
    'exact_match': False,
    'exact_set_template': "",
    'exact_set_match': False,
    'superset_template': "",
    'superset_match': False,
}

match_types = ["exact", "exact_set", "superset"]

def add_small_frag_to_prec(prec, small_frag):

    if small_frag:
        return canonicalize_smiles(prec + "." + small_frag)
    else:
        return prec

def check_exact_match(ground_truth: str, canonical_precs: List[str]) -> Tuple[bool, str]:
    """
    Check if there's an exact match between ground truth and any predicted precursor.
    
    Args:
        ground_truth: The ground truth SMILES string
        canonical_precs: List of canonicalized predicted precursors
    Returns:
        match_found
    """
    if ground_truth in canonical_precs:
        return True
    return False

def check_exact_set_match(ground_truth_set: Set[str], canonical_precs: List[str]) -> Tuple[bool, str]:
    """
    Check if there's an exact set match between ground truth and any predicted precursor.
    
    Args:
        ground_truth_set: Set of ground truth SMILES strings
        canonical_precs: List of canonicalized predicted precursors
        
    Returns:
        match_found
    """
    for pred in canonical_precs:
        predicted_set = set(pred.split('.'))
        if predicted_set == ground_truth_set:
            return True
    return False

def check_superset_match(ground_truth_set: Set[str], canonical_precs: List[str]) -> Tuple[bool, str]:
    """
    Check if there's a superset match between ground truth and any predicted precursor.
    
    Args:
        ground_truth_set: Set of ground truth SMILES strings
        canonical_precs: List of canonicalized predicted precursors
        
    Returns:
        match_found
    """
    for pred in canonical_precs:
        predicted_set = set(pred.split('.'))
        if ground_truth_set.issubset(predicted_set):
            return True
    return False


def match_templates(ground_truth, ground_truth_set, canonical_precs, template):
    """
    Checks for exact, exact set, and superset matches.
    Returns a dict with match results and the matching template(s).
    """
    result = {
        'exact_template': "",
        'exact_match': False,
        'exact_set_template': "",
        'exact_set_match': False,
        'superset_template': "",
        'superset_match': False,
    }

    # Check exact match
    if check_exact_match(ground_truth, canonical_precs):
        for match_type in match_types:
            result[f'{match_type}_template'] = template["reaction_smarts"]
            result[f'{match_type}_match'] = True
        return result

    # Check exact set match
    if check_exact_set_match(ground_truth_set, canonical_precs):
        for match_type in match_types[1:]:
            result[f'{match_type}_template'] = template["reaction_smarts"]
            result[f'{match_type}_match'] = True

    # Check superset match
    if check_superset_match(ground_truth_set, canonical_precs):
        result['superset_template'] = template["reaction_smarts"]
        result['superset_match'] = True

    return result


def process_row(args):
    # set num_cores to 0 to not use pool
    row, templates, num_cores = args
    result = {
        **row,
        **default_result,
    }
    try:
        product = row['product']
        ground_truth = row['reactants']
        candidate_small_frag = ""
        value = row.get('candidate_small_frag', None)
        if value is not None and not pd.isna(value):
            candidate_small_frag = value

        superset_template = ""
        ground_truth_set = set(ground_truth.split('.'))

        
        if num_cores > 0:
            # Parallelize over templates using ProcessPool
            with ProcessPool(max_workers=num_cores) as pool:
                futures = []

                for template in templates:
                    future = pool.schedule(
                        apply_template_to_product,
                        args=(template, product, None, candidate_small_frag),
                        timeout=10
                    )
                    futures.append((template, future))

                for template, future in tqdm(futures, desc=f"Applying templates to row {row['id']}"):
                    try:
                        canonical_precs = future.result()
                    except TimeoutError:
                        logging.warning(f"Template {template} timed out after 10 seconds")
                        continue
                    except Exception as e:
                        logging.error(f"Error applying template {template}: {e}")
                        continue

                    match_result = match_templates(
                        ground_truth, ground_truth_set, canonical_precs, template
                    )
                    for key in match_result:
                        if match_result[key]:
                            result[key] = match_result[key]
                    if match_result['exact_match']:
                        break
        else:
            # Serial version (no parallelization)
            prod = rdchiralReactants(product)
            for template in templates:
                canonical_precs = apply_template_to_product(
                    template, product, prod=prod, candidate_small_frag=candidate_small_frag
                )
                match_result = match_templates(
                    ground_truth, ground_truth_set, canonical_precs, template
                )
                for key in match_result:
                    if match_result[key]:
                        result[key] = match_result[key]
                if match_result['exact_match']:
                    break

    except Exception as e:
        logging.error(f"Error processing row {row['id']}: {e}")
    return result

class TemplateEnumerator:
    @staticmethod
    def load_templates(template_file: str) -> List[str]:
        """Load templates from a jsonl file and format them for rdchiral."""
        logging.info(f"Loading templates from {template_file}...")
        templates, _ = load_templates_as_list(template_file)
        
        for template in templates:
            template['reaction_smarts'] = "("+template['reaction_smarts'].replace(">>", ")>>(")+")"

        logging.info(f"Loaded {len(templates)} templates")
        return templates

    def __init__(self, template_file: str):
        """Initialize the TemplateEnumerator with a template file."""
        self.templates = self.load_templates(template_file)

    def check_columns(self, df: pd.DataFrame) -> bool:
        """Check if the input file has the required columns."""
        required_columns = ['id', 'product', 'reactants']
        is_updated = False
        if (
            'product' not in df.columns 
            or 'reactants' not in df.columns 
        ):
            if 'rxn_smiles' not in df.columns:
                raise ValueError("Either product or reactants or rxn_smiles column must be present")
            else:
                logging.info(f"Converting rxn_smiles to product and reactants")
                tqdm.pandas()
                df['product'] = df['rxn_smiles'].progress_apply(
                    lambda x: canonicalize_smiles(x.split('>')[-1])
                )
                df['reactants'] = df['rxn_smiles'].progress_apply(
                    lambda x: canonicalize_smiles(x.split('>')[0])
                )
                is_updated = True
        
        if len(df) != len(set(df['id'])):
            raise ValueError("Duplicate ids found in input file")

        return df, is_updated

    def load_data_and_check_columns(self, input_file: str) -> pd.DataFrame:
        df = pd.read_csv(input_file)
        logging.info(f"Loaded {len(df)} rows with columns {df.columns}")
        df, is_updated = self.check_columns(df)
        logging.info(f"Checking columns passed, final columns: {df.columns}")
        if is_updated:
            df.to_csv(input_file, index=False)
            logging.info(f"Updated input file with {len(df)} rows")
        return df
        
    def process_products(self, 
                        input_file: str,
                        num_cores: int = None,
                        output_file: str = None) -> None:
        """Process product-reactant pairs in parallel using ThreadPoolExecutor and per-template ProcessPool."""
        df = self.load_data_and_check_columns(input_file)
        
        if num_cores is None:
            num_cores = concurrent.futures.cpu_count()
        
        columns = list(df.columns) + list(default_result.keys())
        if output_file:
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            processed_ids = set()
            if output_path.exists():
                logging.info(f"Output file exists, loading processed IDs")
                try:
                    result_df = pd.read_csv(output_file, on_bad_lines='warn')
                    processed_ids = set(result_df['id'].values)
                    logging.info(f"Found {len(processed_ids)} already processed rows in output file")
                except Exception as e:
                    logging.error(f"Error reading output file: {e}")
                    logging.error("Attempting to read file with different parameters...")
                    result_df = pd.read_csv(output_file, on_bad_lines='skip')
                    processed_ids = set(result_df['id'].values)
                    logging.info(f"Found {len(processed_ids)} processed rows after skipping bad lines")
                if len(processed_ids) > 0:
                    result_df.to_csv(output_file, index=False)
                    logging.info(f"Rewrote output file with {len(processed_ids)} processed rows")
            else:
                header_df = pd.DataFrame(columns=columns)
                header_df.to_csv(output_file, index=False)
                logging.info("Created new output file with header")
        
        remaining_df = df[~df['id'].isin(processed_ids)]
        logging.info(f"Total rows: {len(df)}")
        logging.info(f"Remaining rows: {len(remaining_df)}")
        timed_out_ids = []

        if len(remaining_df):
            logging.info(f"First row in remaining_df: {remaining_df['id'].iloc[0]}")
        else:
            logging.info(f"No remaining rows for first pass")

        logging.info(f"Processing {len(remaining_df)} remaining rows")


        process_args = [(row, self.templates, 0) for _, row in remaining_df.iterrows()]
        success_ids = []
        
        # --- pebble.ProcessPool with 600s timeout ---
        # first pass parallelizes across rows (for efficiency)
        logging.info(f"Starting first pass with {num_cores} cores")
        with ProcessPool(max_workers=num_cores) as pool:
            futures = []
            for args in process_args:
                future = pool.schedule(process_row, args=(args,), timeout=600)
                futures.append((args, future))
            for args, future in tqdm(futures, total=len(futures), desc="Processing rows (first pass)"):
                row_id = args[0]['id']
                try:
                    result = future.result()
                    if isinstance(result, Exception):
                        raise result
                    if output_file:
                        result_df = pd.DataFrame([result], columns=columns)
                        result_df.to_csv(output_file, mode='a', header=False, index=False)
                    success_ids.append(row_id)
                except TimeoutError:
                    logging.warning(f"Row {row_id} timed out in process pool")
                    timed_out_ids.append(row_id)
                except Exception as e:
                    logging.error(f"Error processing row {row_id}: {e}")
                    timed_out_ids.append(row_id)
        logging.info(f"Completed process pool pass. {len(success_ids)} rows succeeded, {len(timed_out_ids)} timed out or errored.")

        # --- ThreadPool for timed-out rows ---
        # second pass parallelizes across templates for timed-out rows
        logging.info(f"Starting second pass with {num_cores} cores, {len(timed_out_ids)} timed-out rows")
        if timed_out_ids:
            logging.info(f"Retrying {len(timed_out_ids)} timed-out rows with parallelization across templates")
            retry_df = df[df['id'].isin(timed_out_ids)]
            process_args_retry = [(row, self.templates, num_cores) for _, row in retry_df.iterrows()]
            for args in process_args_retry:
                try:
                    result = process_row(args)
                    if output_file:
                        result_df = pd.DataFrame([result], columns=columns)
                        result_df.to_csv(output_file, mode='a', header=False, index=False)
                    success_ids.append(args[0]['id'])
                except Exception as e:
                    logging.error(f"Error processing row {args[0]['id']} in thread pool: {e}")

        logging.info(f"Completed processing all rows")
        logging.info(f"Successfully processed {len(success_ids)} new rows")
        logging.info(f"Failed to process {len(process_args) - len(success_ids)} rows")

    


if __name__ == "__main__":
    import argparse
    import os
    
    parser = argparse.ArgumentParser(description='Run template enumeration on a dataset')
    parser.add_argument('--template_file', type=str, required=True,
                      help='Path to the template jsonl file')
    parser.add_argument('--input_file', type=str, required=True,
                      help='Path to input CSV file with product and reactant columns')
    parser.add_argument('--output_file', type=str, required=True,
                      help='Path to save results')
    parser.add_argument('--num_cores', type=int, default=None,
                      help='Number of CPU cores to use (default: all available)')
    parser.add_argument('--log_file', type=str, default=None,
                      help='Path to log file (default: None, logs to console)')
    args = parser.parse_args()

    # Setup logging
    setup_logging(log_file=args.log_file)

    # Initialize enumerator
    enumerator = TemplateEnumerator(args.template_file)
    
    # Process products and save to output file
    logging.info("Processing products...")
    enumerator.process_products(
        input_file=args.input_file,
        num_cores=args.num_cores,
        output_file=args.output_file
    )
    
    # Load results from output file
    logging.info(f"Loading results from {args.output_file}...")
    try:
        results_df = pd.read_csv(args.output_file, on_bad_lines='warn')
        logging.info(f"Loaded {len(results_df)} results from output file")
    except Exception as e:
        logging.error(f"Error reading output file: {e}")
        logging.error("Attempting to read file with different parameters...")
        results_df = pd.read_csv(args.output_file, on_bad_lines='skip')
        logging.info(f"Loaded {len(results_df)} results from output file after skipping bad lines")
    
    # Print statistics
    total = len(results_df)
    exact_matches = results_df['exact_match'].sum()
    exact_set_matches = results_df['exact_set_match'].sum()
    superset_matches = results_df['superset_match'].sum()
    
    logging.info("\nResults Summary:")
    logging.info(f"Total reactions processed: {total}")
    logging.info(f"Exact matches found: {exact_matches} ({exact_matches/total*100:.1f}%)")
    logging.info(f"Exact set matches found: {exact_set_matches} ({exact_set_matches/total*100:.1f}%)")
    logging.info(f"Superset matches found: {superset_matches} ({superset_matches/total*100:.1f}%)")
