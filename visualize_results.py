#!/usr/bin/env python3
"""
Visualization of meta-analysis results
Author: Edison Scientific
Description: Create volcano plots, forest plots, and Kaplan-Meier curves
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from pathlib import Path
import sys

try:
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import multivariate_logrank_test
except ImportError:
    print("ERROR: lifelines package not found. Install with: pip install lifelines")
    sys.exit(1)


def create_volcano_plot(meta_df: pd.DataFrame, subtype: str, output_dir: str):
    """
    Create volcano plot for a specific subtype
    
    Parameters:
    -----------
    meta_df : DataFrame
        Meta-analysis results
    subtype : str
        Subtype to plot
    output_dir : str
        Directory to save plots
    """
    # Filter to subtype
    sub_df = meta_df[meta_df['subtype'] == subtype].copy()
    
    if len(sub_df) == 0:
        print(f"No results for subtype: {subtype}")
        return
    
    # Calculate -log10(p)
    sub_df['neg_log_p'] = -np.log10(sub_df['meta_p'])
    
    # Determine significance
    sub_df['significant'] = sub_df['meta_p'] < 0.05
    
    # Create plot
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Plot non-significant points
    non_sig = sub_df[~sub_df['significant']]
    ax.scatter(non_sig['meta_coef'], non_sig['neg_log_p'], 
               c='gray', alpha=0.5, s=20, label='Non-significant')
    
    # Plot significant points
    sig = sub_df[sub_df['significant']]
    ax.scatter(sig['meta_coef'], sig['neg_log_p'], 
               c='red', alpha=0.7, s=30, label='Significant (p<0.05)')
    
    # Add significance line
    ax.axhline(-np.log10(0.05), color='blue', linestyle='--', linewidth=1, alpha=0.5, label='p=0.05')
    ax.axvline(0, color='black', linestyle='-', linewidth=0.5, alpha=0.3)
    
    # Labels
    ax.set_xlabel('Effect Size (log HR)', fontsize=12)
    ax.set_ylabel('-log10(p-value)', fontsize=12)
    ax.set_title(f'Volcano Plot: {subtype.replace("_", " ").title()}', fontsize=14, fontweight='bold')
    ax.legend()
    
    # Save
    output_file = f"{output_dir}/volcano_{subtype}.png"
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {output_file}")


def create_forest_plot(meta_df: pd.DataFrame, model: str, subtype: str, 
                       dataset_results: dict, output_dir: str):
    """
    Create forest plot for a specific model and subtype showing all datasets
    
    Parameters:
    -----------
    meta_df : DataFrame
        Meta-analysis results (single row for this model+subtype)
    model : str
        PGS model name
    subtype : str
        Subtype name
    dataset_results : dict
        Individual dataset results
    output_dir : str
        Directory to save plots
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Get meta-analysis result
    meta_row = meta_df[(meta_df['model'] == model) & (meta_df['subtype'] == subtype)]
    
    if len(meta_row) == 0:
        print(f"No meta-analysis result for {model} - {subtype}")
        return
    
    meta_row = meta_row.iloc[0]
    
    # Collect dataset results
    datasets = []
    hrs = []
    lowers = []
    uppers = []
    
    for ds in ['cidr', 'i370', 'onco', 'tcga']:
        coef_col = f"{ds}_coef"
        se_col = f"{ds}_se"
        n_col = f"{ds}_n"
        
        if coef_col in meta_row.index and pd.notna(meta_row[coef_col]):
            datasets.append(ds.upper())
            
            coef = meta_row[coef_col]
            se = meta_row[se_col]
            n = meta_row[n_col] if pd.notna(meta_row[n_col]) else 0
            
            hr = np.exp(coef)
            lower = np.exp(coef - 1.96 * se)
            upper = np.exp(coef + 1.96 * se)
            
            hrs.append(hr)
            lowers.append(lower)
            uppers.append(upper)
    
    # Add meta-analysis result
    datasets.append('META-ANALYSIS')
    hrs.append(meta_row['meta_exp_coef'])
    lowers.append(meta_row['meta_exp_lower_95'])
    uppers.append(meta_row['meta_exp_upper_95'])
    
    # Create forest plot
    y_pos = np.arange(len(datasets))
    
    # Plot individual studies
    for i in range(len(datasets) - 1):
        ax.plot([lowers[i], uppers[i]], [y_pos[i], y_pos[i]], 'k-', linewidth=1.5)
        ax.plot(hrs[i], y_pos[i], 'ks', markersize=8)
    
    # Plot meta-analysis (diamond shape)
    i = len(datasets) - 1
    ax.plot([lowers[i], uppers[i]], [y_pos[i], y_pos[i]], 'b-', linewidth=2.5)
    ax.plot(hrs[i], y_pos[i], 'bD', markersize=10)
    
    # Reference line at HR=1
    ax.axvline(1, color='red', linestyle='--', linewidth=1, alpha=0.5)
    
    # Labels
    ax.set_yticks(y_pos)
    ax.set_yticklabels(datasets)
    ax.set_xlabel('Hazard Ratio (95% CI)', fontsize=12)
    ax.set_title(f'{model}\n{subtype.replace("_", " ").title()}', fontsize=12, fontweight='bold')
    ax.set_xlim([min(lowers) * 0.8, max(uppers) * 1.2])
    
    # Add text annotations
    for i, (hr, lower, upper) in enumerate(zip(hrs, lowers, uppers)):
        text = f"{hr:.2f} ({lower:.2f}-{upper:.2f})"
        ax.text(max(uppers) * 1.15, y_pos[i], text, va='center', fontsize=9)
    
    plt.tight_layout()
    
    # Save
    safe_model = model.replace('/', '_').replace(' ', '_')
    output_file = f"{output_dir}/forest_{safe_model}_{subtype}.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {output_file}")


def create_km_curves(model: str, subtype: str, dataset_name: str,
                     scores_file: str, covariates_file: str,
                     filters: dict, output_dir: str):
    """
    Create Kaplan-Meier curves comparing high vs low PGS
    
    Parameters:
    -----------
    model : str
        PGS model name
    subtype : str
        Subtype name
    dataset_name : str
        Dataset name
    scores_file : str
        Path to scores file
    covariates_file : str
        Path to covariates file
    filters : dict
        Filters for subtype
    output_dir : str
        Directory to save plots
    """
    # Load data
    scores = pd.read_csv(scores_file, sep=',')
    covariates = pd.read_csv(covariates_file, sep=',')
    
    # Rename for merging
    sample_col = 'sample' if 'sample' in scores.columns else scores.columns[0]
    scores.rename(columns={sample_col: 'IID'}, inplace=True)
    
    # Filter to subtype
    df = covariates.copy()
    for key, value in filters.items():
        if key in df.columns:
            df = df[df[key] == value]
    
    # Merge with scores
    df = df.merge(scores[['IID', model]], on='IID', how='inner')
    
    # Remove missing data
    df = df.dropna(subset=['survdays', 'vstatus', model])
    
    if len(df) < 20:
        print(f"Insufficient samples for KM plot: {dataset_name} - {model} - {subtype}")
        return
    
    # Split into high/low based on median
    median_score = df[model].median()
    df['risk_group'] = (df[model] > median_score).map({True: 'High PRS', False: 'Low PRS'})
    
    # Create KM plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    kmf = KaplanMeierFitter()
    
    for group in ['Low PRS', 'High PRS']:
        mask = df['risk_group'] == group
        kmf.fit(df.loc[mask, 'survdays'], 
                df.loc[mask, 'vstatus'],
                label=f'{group} (n={mask.sum()})')
        kmf.plot_survival_function(ax=ax, ci_show=True)
    
    # Log-rank test
    try:
        result = multivariate_logrank_test(
            df['survdays'],
            df['risk_group'],
            df['vstatus']
        )
        p_value = result.p_value
    except:
        p_value = np.nan
    
    # Labels
    ax.set_xlabel('Time (days)', fontsize=12)
    ax.set_ylabel('Survival Probability', fontsize=12)
    title = f'{model} - {subtype.replace("_", " ").title()}\n{dataset_name.upper()}'
    if not np.isnan(p_value):
        title += f' (log-rank p={p_value:.3e})'
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.legend(loc='best')
    
    plt.tight_layout()
    
    # Save
    safe_model = model.replace('/', '_').replace(' ', '_')
    output_file = f"{output_dir}/km_{dataset_name}_{safe_model}_{subtype}.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {output_file}")


def main():
    parser = argparse.ArgumentParser(description='Visualize PGS survival results')
    parser.add_argument('--meta', required=True, help='Meta-analysis results file')
    parser.add_argument('--output-dir', required=True, help='Output directory for plots')
    parser.add_argument('--top-n', type=int, default=20, help='Number of top hits to plot')
    parser.add_argument('--data-dir', help='Directory with original data files (for KM plots)')
    
    args = parser.parse_args()
    
    # Create output directory
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    
    # Load meta-analysis results
    meta_df = pd.read_csv(args.meta, sep='\t')
    print(f"Loaded {len(meta_df)} meta-analysis results")
    
    # Get unique subtypes
    subtypes = meta_df['subtype'].unique()
    print(f"\nCreating volcano plots for {len(subtypes)} subtypes...")
    
    for subtype in subtypes:
        create_volcano_plot(meta_df, subtype, args.output_dir)
    
    # Get top N hits overall
    top_hits = meta_df.nsmallest(args.top_n, 'meta_p')
    print(f"\nCreating forest plots for top {len(top_hits)} hits...")
    
    for i, row in top_hits.iterrows():
        create_forest_plot(meta_df, row['model'], row['subtype'], {}, args.output_dir)
    
    # KM plots if data directory provided
    if args.data_dir:
        print(f"\nCreating KM curves for top {len(top_hits)} hits...")
        
        # Define subtype filters
        SUBTYPES = {
            'all_cases': {'case': 1},
            'idh_wildtype': {'case': 1, 'idh': 0},
            'idh_mutant': {'case': 1, 'idh': 1},
            'hgg_idh_wildtype': {'case': 1, 'grade': 'HGG', 'idh': 0},
            'hgg_idh_mutant': {'case': 1, 'grade': 'HGG', 'idh': 1},
            'lgg_idh_wildtype': {'case': 1, 'grade': 'LGG', 'idh': 0},
            'lgg_idh_mutant': {'case': 1, 'grade': 'LGG', 'idh': 1},
            'lgg_idh_mutant_pq_intact': {'case': 1, 'grade': 'LGG', 'idh': 1, 'pq': 0},
            'lgg_idh_mutant_pq_codel': {'case': 1, 'grade': 'LGG', 'idh': 1, 'pq': 1}
        }
        
        for i, row in top_hits.iterrows():
            model = row['model']
            subtype = row['subtype']
            
            # Create KM for each dataset
            for ds in ['cidr', 'i370', 'onco', 'tcga']:
                coef_col = f"{ds}_coef"
                if coef_col in row.index and pd.notna(row[coef_col]):
                    scores_file = f"{args.data_dir}/{ds}.scores.z-scores.txt.gz"
                    cov_file = f"{args.data_dir}/{ds}-covariates.tsv"
                    
                    if Path(scores_file).exists() and Path(cov_file).exists():
                        try:
                            create_km_curves(model, subtype, ds, scores_file, cov_file,
                                           SUBTYPES[subtype], args.output_dir)
                        except Exception as e:
                            print(f"Error creating KM for {ds}-{model}-{subtype}: {e}")
    
    print(f"\nAll visualizations saved to {args.output_dir}")


if __name__ == '__main__':
    main()
