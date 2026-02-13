#!/usr/bin/env python3
"""
Adaptive dataset splitting with KS-distance analysis.
Splits dataset using document-based splitting and finds the split with the lowest maximum KS-statistic
for all specified complexity features.
"""

import os
import sys
import json
sys.path.append("../")

import logging
logging.getLogger('rdkit').setLevel(logging.ERROR)

import random
import tempfile
import argparse
import pandas as pd
import numpy as np
from datetime import datetime
from scipy.stats import ks_2samp
from utils.document_split import document_split
from utils.file_utils import setup_logging


def perform_ks_analysis(train_data, test_data, cols_to_plot, logger=None):
    """
    Perform KS-test between train and test splits for specified columns.
    
    Args:
        train_data: Training data DataFrame
        test_data: Test data DataFrame
        cols_to_plot: List of column names to analyze
        logger: Logger instance to use for logging
    
    Returns:
        dict: Dictionary with KS statistics for each column
    """
    results = {}
    
    for col in cols_to_plot:
        if col in train_data.columns and col in test_data.columns:
            data1 = train_data[col].dropna().values
            data2 = test_data[col].dropna().values
            
            if len(data1) > 0 and len(data2) > 0:
                res = ks_2samp(data1, data2, alternative="two-sided")
                results[col] = {
                    'statistic': res.statistic
                }
            else:
                results[col] = {
                    'statistic': np.nan
                }
        else:
            if logger:
                logger.warning(f"Column {col} not found in data")
            results[col] = {
                'statistic': np.nan
            }
    
    return results


def adaptive_split_with_ks_test(input_file, features_file, output_dir, 
                               test_val_split=0.25, test_split=0.8, 
                               max_iterations=100,
                               cols_to_plot=None, random_seed=None, 
                               log_dir=None, log_level=logging.INFO):
    """
    Perform adaptive dataset splitting with KS-distance analysis.
    
    Args:
        input_file: Path to input JSON file with reactions
        features_file: Path to CSV file with complexity features
        output_dir: Directory to save split CSV files
        test_val_split: Fraction of reactions for validation + test
        test_split: Fraction of validation + test set for test set
        max_iterations: Maximum number of splitting attempts
        cols_to_plot: List of columns to analyze (if None, uses default)
        random_seed: Random seed for reproducibility
    """
    
    # Setup logging
    log_dir_to_use = log_dir if log_dir else output_dir
    log_file, logger = setup_logging(log_dir_to_use, "document_split_with_ks_dist", log_level)
    logger.info(f"Log file created: {log_file}")
    
    if cols_to_plot is None:
        cols_to_plot = ['num_react_atom_forward', 'num_react_atom_retro', 
                       'num_change_ring']
    
    # Load features data
    logger.info(f"Loading features from {features_file}")
    features_df = pd.read_csv(features_file)
    
    
    best_split = None
    best_ks_statistic = None
    best_iteration = 0
    
    logger.info(f"Starting adaptive splitting")
    logger.info(f"Analyzing columns: {cols_to_plot}")
    
    for iteration in range(max_iterations):
        logger.info("--------------------------------")
        logger.info(f"--- Iteration {iteration + 1} ---")
        
        # Use different random seed for each iteration
        current_seed = random_seed + iteration if random_seed is not None else iteration
        
        # Create temporary directory for this iteration
        with tempfile.TemporaryDirectory() as temp_dir:
            # Perform document split with return_dataframes=True
            logger.info(f"Current seed: {current_seed}")
            result = document_split(
                input_file, temp_dir, test_val_split, test_split, return_dataframes=True, random_state=current_seed
            )
            if result is None:
                logger.error("Document split failed")
                continue
            train_data, val_data, test_data = result
        
        # Merge with features
        train_with_features = train_data.merge(features_df, on=["id", "rxn_smiles"], how="left")
        test_with_features = test_data.merge(features_df, on=["id", "rxn_smiles"], how="left")
        
        # Perform KS analysis
        ks_results = perform_ks_analysis(train_with_features, test_with_features, cols_to_plot, logger)
        
        min_ks_statistic = float('inf')
        max_ks_statistic = float('-inf')
        
        logger.info("KS-test results:")
        for col, result in ks_results.items():
            statistic = result['statistic']
            logger.info(f"  {col}: KS-stat={statistic:.6f}")
            
            if not np.isnan(statistic):
                min_ks_statistic = min(min_ks_statistic, statistic)
                max_ks_statistic = max(max_ks_statistic, statistic)
        
        logger.info(f"Minimum KS-statistic: {min_ks_statistic:.6f}")
        logger.info(f"Maximum KS-statistic: {max_ks_statistic:.6f}")
        # Update best split if this is better
        if best_ks_statistic is None or max_ks_statistic < best_ks_statistic:
            best_split = (train_data, val_data, test_data)
            best_ks_statistic = max_ks_statistic
            best_iteration = iteration + 1
            logger.info(f"New best split found! Max KS-statistic: {max_ks_statistic:.4f}")

    if best_split is None:
        logger.error("Failed to find any valid split")
        return False
    
    # Save the best split
    train_data, val_data, test_data = best_split
    
    # Create final dataframes with only required columns
    df_train = train_data[["id", "rxn_smiles"]]
    df_val = val_data[["id", "rxn_smiles"]]
    df_test = test_data[["id", "rxn_smiles"]]
    
    # Save to CSV
    os.makedirs(output_dir, exist_ok=True)
    
    df_train.to_csv(os.path.join(output_dir, "raw_train.csv"), index=False)
    df_val.to_csv(os.path.join(output_dir, "raw_val.csv"), index=False)
    df_test.to_csv(os.path.join(output_dir, "raw_test.csv"), index=False)
    
    logger.info("--------------------------------")
    logger.info("Final Results:")
    logger.info(f"Best split found at iteration {best_iteration}")
    logger.info(f"Best (lowest) maximum KS-statistic: {best_ks_statistic:.4f}")
    logger.info(f"Files saved to: {output_dir}")
    
    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform adaptive dataset splitting with KS-distance analysis.")
    parser.add_argument('--input', type=str, required=True, help='Path to input JSON file with reactions')
    parser.add_argument('--features', type=str, required=True, help='Path to CSV file with complexity features')
    parser.add_argument('--output', type=str, required=True, help='Output directory to save split CSV files')
    parser.add_argument('--test_val_split', type=float, default=0.25, help='Fraction of reactions for validation + test')
    parser.add_argument('--test_split', type=float, default=0.8, help='Fraction of validation + test set for test set')
    parser.add_argument('--max_iterations', type=int, default=100, help='Maximum number of splitting attempts')
    parser.add_argument('--random_seed', type=int, default=42, help='Random seed for reproducibility')
    parser.add_argument('--log_level', type=str, default='INFO', 
                       choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], 
                       help='Logging level')
    parser.add_argument('--log_dir', type=str, default=None, 
                       help='Directory to save log files (default: output directory)')
    
    args = parser.parse_args()
    
    # Define columns to analyze
    cols_to_plot = ['num_react_atom_forward', 'num_react_atom_retro', 'num_change_ring']
    
    success = adaptive_split_with_ks_test(
        input_file=args.input,
        features_file=args.features,
        output_dir=args.output,
        test_val_split=args.test_val_split,
        test_split=args.test_split,
        max_iterations=args.max_iterations,
        cols_to_plot=cols_to_plot,
        random_seed=args.random_seed,
        log_dir=args.log_dir,
        log_level=getattr(logging, args.log_level)
    )
    
    if success:
        print("Script completed successfully!")
    else:
        print("Script failed!")
        sys.exit(1) 