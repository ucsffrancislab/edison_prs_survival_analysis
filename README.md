# PGS Survival Analysis Pipeline

Comprehensive Cox proportional hazards survival analysis for polygenic risk scores (PGS) in glioma datasets.

## Overview

This pipeline performs survival analysis on 4 glioma datasets (cidr, i370, onco, tcga) across 9 molecular subtypes using 5110+ PGS models. It includes:

- Cox proportional hazards regression for each PGS model
- Adjustment for clinical covariates and principal components
- Meta-analysis across datasets using inverse variance weighting
- Visualization: volcano plots, forest plots, and Kaplan-Meier curves

## Requirements

### Python Packages

```bash
pip install pandas numpy scipy matplotlib seaborn lifelines
```

Or if using conda:

```bash
conda install pandas numpy scipy matplotlib seaborn
pip install lifelines
```

### System Requirements

- Python 3.7+
- For full analysis: 32-64 GB RAM per dataset
- For parallel analysis: 64 CPU cores, 128 GB RAM recommended

## Files

### Analysis Scripts

- `survival_analysis.py` - Main survival analysis script
- `meta_analysis.py` - Meta-analysis across datasets
- `visualize_results.py` - Generate plots
- `create_model_list.py` - Helper to create model lists

### Execution Scripts

- `test_pipeline.sh` - Test pipeline on subset (RUN THIS FIRST)
- `archive/run_survival_analysis.sh` - SLURM array job script
- `run_pipeline.sh` - SLURM single parallel job script

## Quick Start

### 1. Test the Pipeline

Before running the full analysis, test on a small subset:

```bash
# Make scripts executable
chmod +x *.sh *.py

# Run test (analyzes ~3% of models, takes ~10-30 minutes)
./test_pipeline.sh
```

This will:
- Create a test model list with ~162 models (3% of 5110)
- Run survival analysis on all 4 datasets
- Perform meta-analysis
- Generate visualization plots
- Save results to `test_results/`

### 2. Review Test Results

```bash
# View top results
head -n 20 test_results/meta_analysis_results.txt

# Check plots
ls test_results/plots/
```

### 3. Run Full Analysis

Once satisfied with test results, create full model list and submit to SLURM:

```bash
# Create full model list (all 5110 models)
python3 create_model_list.py \
    --scores cidr.scores.z-scores.txt.gz \
    --output model_list.txt

# Option A: Array job (4 separate jobs, one per dataset)
sbatch archve/run_survival_analysis.sh

# Option B: Single parallel job (all datasets in one job)
sbatch run_pipeline.sh
```

## Analysis Details

### Subtypes Analyzed

1. **all_cases** - All glioma cases
2. **idh_wildtype** - IDH wildtype cases
3. **idh_mutant** - IDH mutant cases
4. **hgg_idh_wildtype** - High-grade glioma, IDH wildtype
5. **hgg_idh_mutant** - High-grade glioma, IDH mutant
6. **lgg_idh_wildtype** - Low-grade glioma, IDH wildtype
7. **lgg_idh_mutant** - Low-grade glioma, IDH mutant
8. **lgg_idh_mutant_pq_intact** - LGG, IDH mutant, 1p/19q intact
9. **lgg_idh_mutant_pq_codel** - LGG, IDH mutant, 1p/19q codeleted

### Covariates

All models are adjusted for:
- `source` - Data source (one-hot encoded)
- `age` - Age at diagnosis
- `sex` - Sex (M/F)
- `grade` - Tumor grade (HGG/LGG)
- `rad` - Radiation therapy
- `chemo` - Chemotherapy
- `PC1-PC8` - First 8 genetic principal components

### Statistical Methods

- **Cox Proportional Hazards**: Evaluates association between PGS and survival time
- **Meta-analysis**: Fixed-effects inverse variance weighted method
- **Multiple Testing**: Not corrected in meta-analysis (report all p-values)

## Output Files

### Per-Dataset Results

`results/{dataset}_survival_results.txt` - Tab-delimited file with columns:
- `dataset` - Dataset name
- `model` - PGS model name
- `subtype` - Molecular subtype
- `n_samples` - Number of samples
- `n_events` - Number of death events
- `coef` - Log hazard ratio
- `exp_coef` - Hazard ratio
- `se_coef` - Standard error
- `z` - Z-statistic
- `p` - P-value
- `lower_95`, `upper_95` - 95% confidence intervals

### Meta-Analysis Results

`results/meta_analysis_results.txt` - Combined results across datasets:
- All columns from individual datasets (prefixed with dataset name)
- `meta_coef` - Meta-analyzed log hazard ratio
- `meta_exp_coef` - Meta-analyzed hazard ratio
- `meta_se` - Meta-analysis standard error
- `meta_z` - Meta-analysis Z-statistic
- `meta_p` - Meta-analysis P-value
- `Q` - Cochran's Q statistic for heterogeneity
- `I2` - I-squared heterogeneity metric
- `n_studies` - Number of datasets contributing

### Visualizations

`results/plots/` directory contains:

1. **Volcano plots** (`volcano_{subtype}.png`)
   - One per subtype
   - X-axis: Effect size (log HR)
   - Y-axis: -log10(p-value)
   - Highlights significant associations (p<0.05)

2. **Forest plots** (`forest_{model}_{subtype}.png`)
   - Top 20 associations
   - Shows effect estimates from each dataset
   - Meta-analysis estimate (diamond)
   - 95% confidence intervals

3. **Kaplan-Meier curves** (`km_{dataset}_{model}_{subtype}.png`)
   - Top 20 associations, all datasets
   - Compares high vs low PGS (median split)
   - Includes log-rank test p-value

## Customization

### Analyze Different Model Subsets

```bash
# Analyze only custom models
grep "scoring_system" test_model_list.txt > custom_models.txt

# Analyze specific PGS catalog range
sed -n '/PGS001000/,/PGS002000/p' model_list.txt > pgs_subset.txt

# Run with custom list
python3 survival_analysis.py \
    --dataset cidr \
    --scores cidr.scores.z-scores.txt.gz \
    --covariates cidr-covariates.tsv \
    --models custom_models.txt \
    --output results/custom_results.txt
```

### Adjust Resource Allocation

Edit SLURM scripts to modify:
```bash
#SBATCH --cpus-per-task=32    # Fewer CPUs
#SBATCH --mem=64G             # Less memory
#SBATCH --time=12:00:00       # Shorter time
```

### Change Top N Visualizations

```bash
python3 visualize_results.py \
    --meta results/meta_analysis_results.txt \
    --output-dir results/plots \
    --top-n 50 \              # Plot top 50 instead of 20
    --data-dir .
```

## Troubleshooting

### Import Errors

```bash
# Install missing packages
pip install lifelines pandas numpy scipy matplotlib seaborn
```

### Memory Issues

Reduce models analyzed:
```bash
# Analyze 1% of models
python3 create_model_list.py --scores cidr.scores.z-scores.txt.gz --output model_list_1pct.txt --subset 0.01
```

### Convergence Warnings

Some Cox models may not converge (small sample sizes, collinearity). These are automatically skipped and logged as `None` in results.

### Missing Data

The pipeline automatically handles:
- Samples with missing survival data (excluded)
- Samples with missing covariates (excluded)
- Models with insufficient samples (<10, skipped)

## Performance Estimates

Based on 4 datasets × 5110 models × 9 subtypes = ~184,000 Cox models:

- **Test (3% models)**: 10-30 minutes per dataset
- **Full analysis per dataset**: 6-12 hours
- **Array job (parallel)**: 6-12 hours total
- **Single job (sequential)**: 24-48 hours total
- **Meta-analysis**: 5-10 minutes
- **Visualizations**: 30-60 minutes

## Expected Results

### Sample Sizes by Dataset/Subtype

Approximate sample sizes:
- **all_cases**: 446-2165 per dataset
- **idh_wildtype**: 215-786 per dataset
- **idh_mutant**: 111-655 per dataset
- **Specific subtypes**: 10-400 per dataset

Some rare subtypes may have insufficient samples in some datasets.

### Statistical Power

- Large effects (HR>2): Well powered in most subtypes
- Moderate effects (HR=1.5): Requires larger subtypes
- Small effects (HR<1.3): May be underpowered in rare subtypes

## Citation

If you use this pipeline, please cite:
- **lifelines**: Davidson-Pilon, C. (2019). lifelines: survival analysis in Python. Journal of Open Source Software, 4(40), 1317.
- **PGS Catalog**: Polygenic Score Catalog (https://www.pgscatalog.org/)

## Support

For questions or issues:
1. Check test results first
2. Review error logs in `logs/` directory
3. Verify data file formats match expected structure

## License

This pipeline is provided as-is for research use.

https://platform.edisonscientific.com/trajectories/beee0ba7-2d74-46b3-9aaf-b4af19e7fa30


