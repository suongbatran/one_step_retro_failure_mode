import os
import json
import random
import logging
import concurrent.futures
import pandas as pd

from tqdm import tqdm
from pathlib import Path
from rdkit import Chem, RDLogger
from rdchiral.main import rdchiralReaction, rdchiralReactants, rdchiralRun

from typing import List, Tuple, Optional, Set, Dict, Any, Union

from pebble import ProcessPool, ThreadPool
from concurrent.futures import TimeoutError, ThreadPoolExecutor
from utils.templ_enum_utils import has_passed_filter, apply_template_to_product
from utils.file_utils import setup_logging, load_templates_as_list
from utils.chem_utils import canonicalize_smiles

RDLogger.DisableLog('rdApp.*')
default_result = {
    # default template set (min freq = 5)
    'num_applied_default': 0,
    'num_all_precs_default': 0,
    'num_unique_precs_default': 0,
    # all templates (all)
    'num_applied_all': 0,
    'num_all_precs_all': 0,
    'num_unique_precs_all': 0,
}


def process_row(args):
    # set num_cores to 0 to not use pool
    # Support both formats (row, templates, num_cores) and (row, default, filtered_templates, num_cores)
    if len(args) == 4:
        row, default, filtered_templates, num_cores = args
        default = default if default is not None else []
        filtered_templates = filtered_templates if filtered_templates is not None else []
    else:
        row, templates, num_cores = args
        default = templates
        filtered_templates = []
    
    result = {
        **row,
        **default_result,
    }
    try:
        product = row['product']
        all_precs = []
        num_applied_templates = 0
        
        # Initialize prod for serial version (will be reused for both template sets)
        prod = None
        if num_cores == 0:
            prod = rdchiralReactants(product)
    
       
        if num_cores > 0:
            # Parallelize over templates using ProcessPool
            with ProcessPool(max_workers=num_cores) as pool:

                # Apply default first
                futures = []
                for template in default:
                    future = pool.schedule(
                        apply_template_to_product,
                        args=(template, product, None),
                        timeout=10
                    )
                    futures.append((template, future))

                for template, future in tqdm(futures, desc=f"Applying templates to row {row['id']}"):
                    try:
                        canonical_precs = future.result()
                        if canonical_precs:
                            all_precs.extend(canonical_precs)
                            num_applied_templates += 1
                    except TimeoutError:
                        logging.warning(f"Template {template} timed out after 10 seconds")
                        continue
                    except Exception as e:
                        logging.error(f"Error applying template {template}: {e}")
                        continue
                    
                unique_precs = list(set(all_precs))
                result['num_unique_precs_default'] = len(unique_precs)
                result['num_all_precs_default'] = len(all_precs)
                result['num_applied_default'] = num_applied_templates
                unique_precs.clear()
        
                futures = []
                for template in filtered_templates:
                    future = pool.schedule(
                        apply_template_to_product,
                        args=(template, product, None),
                        timeout=10
                    )
                    futures.append((template, future))
                for template, future in tqdm(futures, desc=f"Applying filtered_templates to row {row['id']}"):
                    try:
                        canonical_precs = future.result()
                        if canonical_precs:
                            all_precs.extend(canonical_precs)
                            num_applied_templates += 1
                    except TimeoutError:
                        logging.warning(f"Template {template} timed out after 10 seconds")
                        continue
                    except Exception as e:
                        logging.error(f"Error applying template {template}: {e}")
                        continue
                
                unique_precs = list(set(all_precs))
                result['num_unique_precs_all'] = len(unique_precs)
                result['num_all_precs_all'] = len(all_precs)
                result['num_applied_all'] = num_applied_templates
        else:
            # Serial version (no parallelization)
            for template in default:
                canonical_precs = apply_template_to_product(
                    template, product, prod=prod
                )
                if canonical_precs:
                    all_precs.extend(canonical_precs)
                    num_applied_templates += 1

            unique_precs = list(set(all_precs))
            result['num_unique_precs_default'] = len(unique_precs)
            result['num_all_precs_default'] = len(all_precs)
            result['num_applied_default'] = num_applied_templates
            unique_precs.clear()

            for template in filtered_templates:
                canonical_precs = apply_template_to_product(
                    template, product, prod=prod
                )
                if canonical_precs:
                    all_precs.extend(canonical_precs)
                    num_applied_templates += 1
            
            unique_precs = list(set(all_precs))
            result['num_unique_precs_all'] = len(unique_precs)
            result['num_all_precs_all'] = len(all_precs)
            result['num_applied_all'] = num_applied_templates
    except Exception as e:
        logging.error(f"Error processing row {row['id']}: {e}")
    return result

class TemplateEnumeratorNumPrecs:
    @staticmethod
    def load_templates(template_path: str) -> List[str]:
        """Load templates from a jsonl file and format them for rdchiral."""
        logging.info(f"Loading templates from {template_path}...")
        
        templates, _ = load_templates_as_list(template_path)
        for template in templates:
            template['reaction_smarts'] = "("+template['reaction_smarts'].replace(">>", ")>>(")+")"

        logging.info(f"Loaded {len(templates)} templates")
        return templates

    def __init__(
        self, 
        template_path: str, 
        filtered_template_path: Optional[str] = None, 
        random_seed: int = 42,
        num_cores: int = None
        ):
        """Initialize the TemplateEnumerator.
        
        Args:
            template_path: Path to first template file
            filtered_template_path: Optional path to second template file
            random_seed: Random seed for sampling
            num_cores: Number of cores to use
        """
        self.templates = self.load_templates(template_path)
        self.filtered_templates = self.load_templates(filtered_template_path) if filtered_template_path else []
        self.random_seed = random_seed
        self.num_cores = num_cores


    def load_data(self, input_file: str, sampled_products_file: str, num_samples: int) -> pd.DataFrame:

        if os.path.exists(sampled_products_file):
            df_sampled = pd.read_csv(sampled_products_file)
            logging.info(f"Loaded {len(df_sampled)} sampled products from {sampled_products_file}")
        else:
            df_sampled = self.random_sample_products(input_file, num_samples)
            df_sampled.to_csv(sampled_products_file, index=False)
            logging.info(f"Saved {len(df_sampled)} sampled products to {sampled_products_file}")
        return df_sampled
    

    def random_sample_products(self, input_file: str, num_samples: int) -> pd.DataFrame:
        # Read CSV
        df = pd.read_csv(input_file)
        
        # Generate an 'id' column if not present
        if 'id' not in df.columns:
            df['id'] = range(len(df))
        
        # Determine product column based on file format
        # processed_test.csv has '1' column with products
        # raw_test.csv has 'rxn_smiles' column (reactants>>product)
        if '1' in df.columns:
            product_col = '1'
        elif 'rxn_smiles' in df.columns:
            # Extract product from rxn_smiles (format: reactants>>product)
            df['product'] = df['rxn_smiles'].apply(lambda x: canonicalize_smiles(x.split('>')[-1]))
            product_col = 'product'
        else:
            raise ValueError("Input file must have either '1' column (processed_test.csv) or 'rxn_smiles' column (raw_test.csv)")
        
        products = df[product_col].dropna().unique()
        
        # Randomly sample products
        random.seed(self.random_seed)
        sampled_products = random.sample(list(products), min(num_samples, len(products)))
        
        # Keep first row for each sampled product and save to CSV
        df_sampled = df[df[product_col].isin(sampled_products)].drop_duplicates(subset=product_col, keep='first')
        df_sampled = df_sampled[['id', product_col]].rename(columns={product_col: 'product'})
        
        return df_sampled
        
    def process_products(self, input_file: str, sampled_products_file: str, output_file: str, num_samples: int) -> None:
        """Process product-reactant pairs in parallel using ThreadPoolExecutor and per-template ProcessPool."""
        df = self.load_data(input_file, sampled_products_file, num_samples)
        
        if self.num_cores is None:
            self.num_cores = concurrent.futures.cpu_count()
        
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

        process_args = [(row, self.templates, self.filtered_templates, 0) for _, row in remaining_df.iterrows()]
        success_ids = []
        #
        # --- pebble.ProcessPool with 600s timeout ---
        # first pass parallelizes across rows
        logging.info(f"Starting first pass with {self.num_cores} cores")
        with ProcessPool(max_workers=self.num_cores) as pool:
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
        logging.info(f"Starting second pass with {self.num_cores} cores, {len(timed_out_ids)} timed-out rows")
        if timed_out_ids:
            logging.info(f"Retrying {len(timed_out_ids)} timed-out rows with parallelization across templates")
            retry_df = df[df['id'].isin(timed_out_ids)]
            process_args_retry = [(row, self.templates, self.filtered_templates, self.num_cores) for _, row in retry_df.iterrows()]
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
    import pandas as pd
    
    parser = argparse.ArgumentParser(description='Run template enumeration on a dataset')
    parser.add_argument('--template_path', type=str, required=True,
                      help='Path to the template jsonl file')
    parser.add_argument('--filtered_template_path', type=str, required=False, default=None,
                      help='Path to the second template jsonl file')
    parser.add_argument('--input_file', type=str, required=True,
                      help='Path to input CSV file with product and reactant columns')
    parser.add_argument('--sampled_products_file', type=str, required=True,
                      help='Path to sampled products file')
    parser.add_argument('--output_file', type=str, required=True,
                      help='Path to save results')
    parser.add_argument('--num_samples', type=int, required=False, default=1000,
                      help='Number of samples to take from the input file')
    parser.add_argument('--num_cores', type=int, default=None,
                      help='Number of CPU cores to use (default: all available)')
    parser.add_argument('--log_file', type=str, default=None,
                      help='Path to log file (default: None, logs to console)')
    args = parser.parse_args()

    # Setup logging
    setup_logging(log_file=args.log_file)

    # Initialize enumerator
    enumerator = TemplateEnumeratorNumPrecs(
        template_path=args.template_path,
        filtered_template_path=args.filtered_template_path,
        num_cores=args.num_cores
    )
    
    # Process products and save to output file
    logging.info("Processing products...")
    enumerator.process_products(
        input_file=args.input_file,
        sampled_products_file=args.sampled_products_file,
        num_samples=args.num_samples,
        output_file=args.output_file
    )
    

