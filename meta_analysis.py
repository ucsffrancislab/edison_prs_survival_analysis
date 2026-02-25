#!/usr/bin/env python3
"""
Meta-analysis of survival results across datasets
Author: Edison Scientific
Description: Fixed-effects inverse-variance weighted meta-analysis
"""

import pandas as pd
import numpy as np
import argparse
from pathlib import Path
from typing import List
import sys


def inverse_variance_meta(coef_list, se_list):
    """
    Perform fixed-effects inverse variance weighted meta-analysis
    
    Parameters:
    -----------
    coef_list : array-like
        Effect sizes (beta coefficients) from each study
    se_list : array-like
        Standard errors from each study
    
    Returns:
    --------
    dict with meta-analysis results
    """
    coef_array = np.array(coef_list)
    se_array = np.array(se_list)
    
    # Calculate weights (inverse variance)
    weights = 1 / (se_array ** 2)
    
    # Weighted mean effect
    meta_coef = np.sum(weights * coef_array) / np.sum(weights)
    
    # Standard error of meta effect
    meta_se = np.sqrt(1 / np.sum(weights))
    
    # Z-score and p-value
    meta_z = meta_coef / meta_se
    meta_p = 2 * (1 - stats.norm.cdf(abs(meta_z)))
    
    # Confidence intervals
    meta_lower = meta_coef - 1.96 * meta_se
    meta_upper = meta_coef + 1.96 * meta_se
    
    # Heterogeneity (Cochran's Q)
    Q = np.sum(weights * (coef_array - meta_coef) ** 2)
    df = len(coef_array) - 1
    
    # I-squared
    I2 = max(0, (Q - df) / Q) if Q > 0 else 0
    
    return {
        'meta_coef': meta_coef,
        'meta_se': meta_se,
        'meta_z': meta_z,
        'meta_p': meta_p,
        'meta_lower_95': meta_lower,
        'meta_upper_95': meta_upper,
        'meta_exp_coef': np.exp(meta_coef),
        'meta_exp_lower_95': np.exp(meta_lower),
        'meta_exp_upper_95': np.exp(meta_upper),
        'Q': Q,
        'I2': I2,
        'n_studies': len(coef_array)
    }


def meta_analyze_results(result_files: List[str], output_file: str) -> pd.DataFrame:
    """
    Meta-analyze results from multiple datasets
    
    Parameters:
    -----------
    result_files : list
        List of paths to dataset-specific result files
    output_file : str
        Path to save meta-analysis results
    
    Returns:
    --------
    DataFrame with meta-analysis results
    """
    # Load all results
    all_results = []
    for f in result_files:
        df = pd.read_csv(f, sep='\t')
        all_results.append(df)
    
    # Combine
    combined = pd.concat(all_results, ignore_index=True)
    
    print(f"Loaded {len(all_results)} datasets with {len(combined)} total results")
    
    # Group by model and subtype
    grouped = combined.groupby(['model', 'subtype'])
    
    meta_results = []
    
    for (model, subtype), group in grouped:
        # Need at least 2 studies for meta-analysis
        if len(group) < 2:
            continue
        
        # Get coefficients and standard errors
        coefs = group['coef'].values
        ses = group['se_coef'].values
        
        # Perform meta-analysis
        meta = inverse_variance_meta(coefs, ses)
        
        # Add metadata
        meta['model'] = model
        meta['subtype'] = subtype
        meta['total_samples'] = group['n_samples'].sum()
        meta['total_events'] = group['n_events'].sum()
        
        # Add individual study results
        for i, row in group.iterrows():
            meta[f"{row['dataset']}_coef"] = row['coef']
            meta[f"{row['dataset']}_se"] = row['se_coef']
            meta[f"{row['dataset']}_p"] = row['p']
            meta[f"{row['dataset']}_n"] = row['n_samples']
        
        meta_results.append(meta)
    
    # Convert to DataFrame
    meta_df = pd.DataFrame(meta_results)
    
    # Sort by p-value
    meta_df = meta_df.sort_values('meta_p')
    
    # Save
    meta_df.to_csv(output_file, sep='\t', index=False)
    print(f"Meta-analysis complete: {len(meta_df)} results saved to {output_file}")
    
    return meta_df


def main():
    parser = argparse.ArgumentParser(description='Meta-analysis of PGS survival results')
    parser.add_argument('--input', required=True, nargs='+', help='Input result files from each dataset')
    parser.add_argument('--output', required=True, help='Output file for meta-analysis results')
    
    args = parser.parse_args()
    
    # Import scipy here to avoid import issues
    global stats
    from scipy import stats
    
    print(f"Meta-analyzing {len(args.input)} datasets...")
    meta_analyze_results(args.input, args.output)


if __name__ == '__main__':
    main()
