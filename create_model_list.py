#!/usr/bin/env python3
"""
Create model list for survival analysis
Author: Edison Scientific
Description: Extract model names from scores file and create subsets
"""

import pandas as pd
import argparse
from pathlib import Path


def create_model_list(scores_file: str, output_file: str, subset_pct: float = None):
    """
    Create a list of model names from scores file
    
    Parameters:
    -----------
    scores_file : str
        Path to scores file
    output_file : str
        Path to output model list file
    subset_pct : float, optional
        If provided, select a subset of models (e.g., 0.03 for 3%)
    """
    # Load scores
    scores = pd.read_csv(scores_file, sep=',', nrows=1)
    
    # Get model names (exclude sample column)
    sample_col = 'sample' if 'sample' in scores.columns else scores.columns[0]
    all_models = [col for col in scores.columns if col != sample_col]
    
    print(f"Total models in {scores_file}: {len(all_models)}")
    
    # Separate custom and PGS models
    custom_models = [m for m in all_models if not m.startswith('PGS')]
    pgs_models = [m for m in all_models if m.startswith('PGS')]
    
    print(f"  Custom models: {len(custom_models)}")
    print(f"  PGS models: {len(pgs_models)}")
    
    # Select subset if requested
    if subset_pct is not None:
        # Always include custom models
        # Sample PGS models
        n_pgs_subset = int(len(pgs_models) * subset_pct)
        
        # Sample evenly across the range
        step = max(1, len(pgs_models) // n_pgs_subset)
        pgs_subset = pgs_models[::step][:n_pgs_subset]
        
        selected_models = custom_models + pgs_subset
        print(f"\nSelected {len(selected_models)} models ({len(selected_models)/len(all_models)*100:.1f}%):")
        print(f"  Custom: {len(custom_models)}")
        print(f"  PGS subset: {len(pgs_subset)}")
    else:
        selected_models = all_models
        print(f"\nUsing all {len(selected_models)} models")
    
    # Write to file
    with open(output_file, 'w') as f:
        for model in selected_models:
            f.write(f"{model}\n")
    
    print(f"\nModel list saved to: {output_file}")
    
    return selected_models


def main():
    parser = argparse.ArgumentParser(description='Create model list for survival analysis')
    parser.add_argument('--scores', required=True, help='Path to scores file')
    parser.add_argument('--output', required=True, help='Output file for model list')
    parser.add_argument('--subset', type=float, help='Subset percentage (e.g., 0.03 for 3%)')
    
    args = parser.parse_args()
    
    create_model_list(args.scores, args.output, args.subset)


if __name__ == '__main__':
    main()
