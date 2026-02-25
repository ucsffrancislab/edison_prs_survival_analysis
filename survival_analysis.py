#!/usr/bin/env python3
"""
PRS Survival Analysis for Glioma
Author: Edison Scientific
Description: Cox proportional hazards models for PGS associations with glioma survival
"""

import pandas as pd
import numpy as np
import argparse
import json
import sys
from pathlib import Path
from typing import List, Dict, Tuple
import warnings
warnings.filterwarnings('ignore')

try:
    from lifelines import CoxPHFitter
    from lifelines.exceptions import ConvergenceError
except ImportError:
    print("ERROR: lifelines package not found. Install with: pip install lifelines")
    sys.exit(1)


class SurvivalAnalysis:
    """Class to perform Cox PH regression for PGS models"""
    
    # Define subtypes as class variable
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
    
    COVARIATES = ['source', 'age', 'sex', 'grade', 'rad', 'chemo', 
                  'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8']
    
    def __init__(self, dataset_name: str):
        """Initialize with dataset name"""
        self.dataset_name = dataset_name
        self.scores = None
        self.covariates = None
        
    def load_data(self, scores_file: str, covariates_file: str) -> None:
        """Load scores and covariates"""
        print(f"Loading data for {self.dataset_name}...")
        self.scores = pd.read_csv(scores_file, sep=',')
        self.covariates = pd.read_csv(covariates_file, sep=',')
        
        # Merge on sample ID
        sample_col = 'sample' if 'sample' in self.scores.columns else self.scores.columns[0]
        id_col = 'IID' if 'IID' in self.covariates.columns else self.covariates.columns[0]
        
        self.scores.rename(columns={sample_col: 'IID'}, inplace=True)
        
        print(f"  Scores: {self.scores.shape}")
        print(f"  Covariates: {self.covariates.shape}")
        
    def filter_subtype(self, filters: Dict) -> pd.DataFrame:
        """Filter data based on subtype criteria"""
        df = self.covariates.copy()
        
        for key, value in filters.items():
            if key in df.columns:
                df = df[df[key] == value]
        
        return df
    
    def run_cox_model(self, model_name: str, subtype_name: str, filters: Dict) -> Dict:
        """
        Run Cox PH model for a single PGS model and subtype
        
        Returns dict with results or None if failed
        """
        # Filter to subtype
        subset = self.filter_subtype(filters)
        
        if len(subset) < 10:  # Need minimum samples
            return None
        
        # Merge with scores
        subset = subset.merge(self.scores[['IID', model_name]], on='IID', how='inner')
        
        # Remove rows with missing survival data
        subset = subset.dropna(subset=['survdays', 'vstatus', model_name])
        
        if len(subset) < 10:
            return None
        
        # Prepare data for Cox model
        # Need: survdays, vstatus, PGS score, and covariates
        cox_data = subset[['survdays', 'vstatus', model_name]].copy()
        
        # Add covariates that are available and have data
        available_covariates = []
        for cov in self.COVARIATES:
            if cov in subset.columns:
                # Only add covariate if it has at least some non-missing data
                if subset[cov].notna().sum() > 0:
                    cox_data[cov] = subset[cov]
                    available_covariates.append(cov)
        
        # Handle categorical variables
        if 'sex' in cox_data.columns:
            cox_data['sex'] = (cox_data['sex'] == 'M').astype(int)
        if 'grade' in cox_data.columns:
            cox_data['grade'] = (cox_data['grade'] == 'HGG').astype(int)
        if 'source' in cox_data.columns:
            # One-hot encode source
            source_dummies = pd.get_dummies(cox_data['source'], prefix='source', drop_first=True)
            cox_data = pd.concat([cox_data.drop('source', axis=1), source_dummies], axis=1)
            available_covariates.remove('source')
            available_covariates.extend(source_dummies.columns.tolist())
        
        # Drop rows with any missing data
        # This will remove rows where any of the available covariates are missing
        cox_data = cox_data.dropna()
        
        if len(cox_data) < 10:
            return None
        
        try:
            # Fit Cox PH model
            cph = CoxPHFitter(penalizer=0.0)
            cph.fit(cox_data, duration_col='survdays', event_col='vstatus')
            
            # Extract results for the PGS model
            summary = cph.summary
            
            if model_name in summary.index:
                result = {
                    'dataset': self.dataset_name,
                    'model': model_name,
                    'subtype': subtype_name,
                    'n_samples': len(cox_data),
                    'n_events': int(cox_data['vstatus'].sum()),
                    'coef': summary.loc[model_name, 'coef'],
                    'exp_coef': summary.loc[model_name, 'exp(coef)'],
                    'se_coef': summary.loc[model_name, 'se(coef)'],
                    'z': summary.loc[model_name, 'z'],
                    'p': summary.loc[model_name, 'p'],
                    'lower_95': summary.loc[model_name, 'exp(coef) lower 95%'],
                    'upper_95': summary.loc[model_name, 'exp(coef) upper 95%']
                }
                return result
            else:
                return None
                
        except (ConvergenceError, Exception) as e:
            # Model didn't converge or other error
            return None
    
    def analyze_models(self, models: List[str], output_file: str) -> None:
        """
        Run Cox models for all specified models and subtypes
        Write results incrementally to output file
        """
        results = []
        total_analyses = len(models) * len(self.SUBTYPES)
        counter = 0
        
        print(f"\nRunning {total_analyses} analyses ({len(models)} models × {len(self.SUBTYPES)} subtypes)...")
        
        for model_name in models:
            if model_name not in self.scores.columns:
                print(f"Warning: {model_name} not found in scores")
                continue
            
            for subtype_name, filters in self.SUBTYPES.items():
                counter += 1
                if counter % 50 == 0:
                    print(f"  Progress: {counter}/{total_analyses} ({counter/total_analyses*100:.1f}%)")
                
                result = self.run_cox_model(model_name, subtype_name, filters)
                if result is not None:
                    results.append(result)
        
        # Save results
        if results:
            results_df = pd.DataFrame(results)
            results_df.to_csv(output_file, sep='\t', index=False)
            print(f"\nCompleted: {len(results)} successful analyses saved to {output_file}")
        else:
            print("No results to save")


def main():
    parser = argparse.ArgumentParser(description='PGS Survival Analysis')
    parser.add_argument('--dataset', required=True, help='Dataset name (e.g., cidr, i370, onco, tcga)')
    parser.add_argument('--scores', required=True, help='Path to z-scores file')
    parser.add_argument('--covariates', required=True, help='Path to covariates file')
    parser.add_argument('--models', required=True, help='Comma-separated list of model names OR path to file with model names')
    parser.add_argument('--output', required=True, help='Output file path')
    
    args = parser.parse_args()
    
    # Parse models
    if Path(args.models).exists():
        # Load from file
        with open(args.models, 'r') as f:
            models = [line.strip() for line in f if line.strip()]
    else:
        # Parse from command line
        models = [m.strip() for m in args.models.split(',')]
    
    print(f"Dataset: {args.dataset}")
    print(f"Models to analyze: {len(models)}")
    
    # Run analysis
    analysis = SurvivalAnalysis(args.dataset)
    analysis.load_data(args.scores, args.covariates)
    analysis.analyze_models(models, args.output)


if __name__ == '__main__':
    main()
