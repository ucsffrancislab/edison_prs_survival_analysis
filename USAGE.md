# PGS Survival Analysis - Usage Guide

## Quick Start Commands

### 1. Installation
```bash
chmod +x *.sh *.py
./install_requirements.sh
```

### 2. Test Pipeline (IMPORTANT - Run this first!)
```bash
./test_pipeline.sh
```
This will analyze ~3% of models (~15-30 minutes) and verify everything works.

### 3. Full Analysis

After successful test, choose one option:

**Option A: SLURM Array Job** (Recommended for most HPC clusters)
```bash
# Create full model list
python3 create_model_list.py \
    --scores cidr.scores.z-scores.txt.gz \
    --output model_list.txt

# Submit array job (4 jobs, one per dataset)
sbatch run_survival_analysis.sh
```

**Option B: Single Parallel Job** (For large multi-CPU nodes)
```bash
# Create full model list
python3 create_model_list.py \
    --scores cidr.scores.z-scores.txt.gz \
    --output model_list.txt

# Submit single job that runs all datasets in parallel
sbatch run_parallel_survival.sh
```

### 4. Monitor Progress

```bash
# Check job status
squeue -u $USER

# View logs
tail -f logs/survival_*.out

# Check output
ls -lh results/
```

### 5. View Results

```bash
# Top meta-analysis results
head -n 20 results/meta_analysis_results.txt | column -t

# View plots
ls results/plots/

# Copy plots to local machine for viewing
scp -r your_server:path/to/results/plots/ ./
```

## File Structure

```
.
├── survival_analysis.py         # Main analysis script
├── meta_analysis.py              # Meta-analysis script
├── visualize_results.py          # Plotting script
├── create_model_list.py          # Helper to create model lists
│
├── run_survival_analysis.sh      # SLURM array job script
├── run_parallel_survival.sh      # SLURM parallel job script
├── test_pipeline.sh              # Test script (run first!)
├── install_requirements.sh       # Package installation
│
├── README.md                     # Detailed documentation
├── USAGE.md                      # This file
│
├── Data files:
│   ├── cidr.scores.z-scores.txt.gz
│   ├── cidr-covariates.tsv
│   ├── i370.scores.z-scores.txt.gz
│   ├── i370-covariates.tsv
│   ├── onco.scores.z-scores.txt.gz
│   ├── onco-covariates.tsv
│   ├── tcga.scores.z-scores.txt.gz
│   └── tcga-covariates.tsv
│
└── Output:
    ├── model_list.txt            # List of models to analyze
    ├── results/
    │   ├── cidr_survival_results.txt
    │   ├── i370_survival_results.txt
    │   ├── onco_survival_results.txt
    │   ├── tcga_survival_results.txt
    │   ├── meta_analysis_results.txt
    │   └── plots/
    │       ├── volcano_*.png
    │       ├── forest_*.png
    │       └── km_*.png
    └── logs/
        └── survival_*.{out,err}
```

## Common Tasks

### Analyze Specific Model Subset

```bash
# Create custom model list
echo "PGS000001" > custom_models.txt
echo "PGS000002" >> custom_models.txt
echo "idhmut_scoring_system" >> custom_models.txt

# Run analysis
python3 survival_analysis.py \
    --dataset cidr \
    --scores cidr.scores.z-scores.txt.gz \
    --covariates cidr-covariates.tsv \
    --models custom_models.txt \
    --output custom_results.txt
```

### Analyze Different Percentage of Models

```bash
# 1% of models (~51 models)
python3 create_model_list.py \
    --scores cidr.scores.z-scores.txt.gz \
    --output model_list_1pct.txt \
    --subset 0.01

# 5% of models (~256 models)
python3 create_model_list.py \
    --scores cidr.scores.z-scores.txt.gz \
    --output model_list_5pct.txt \
    --subset 0.05

# Then edit run_survival_analysis.sh to use this file
# Change: MODEL_LIST="model_list_1pct.txt"
```

### Run Single Dataset

```bash
python3 survival_analysis.py \
    --dataset cidr \
    --scores cidr.scores.z-scores.txt.gz \
    --covariates cidr-covariates.tsv \
    --models model_list.txt \
    --output results/cidr_survival_results.txt
```

### Meta-analyze Manually

```bash
python3 meta_analysis.py \
    --input results/cidr_survival_results.txt \
            results/i370_survival_results.txt \
            results/onco_survival_results.txt \
            results/tcga_survival_results.txt \
    --output results/meta_analysis_results.txt
```

### Generate Visualizations

```bash
# Top 20 hits (default)
python3 visualize_results.py \
    --meta results/meta_analysis_results.txt \
    --output-dir results/plots \
    --top-n 20 \
    --data-dir .

# Top 50 hits
python3 visualize_results.py \
    --meta results/meta_analysis_results.txt \
    --output-dir results/plots_top50 \
    --top-n 50 \
    --data-dir .

# Without KM curves (faster)
python3 visualize_results.py \
    --meta results/meta_analysis_results.txt \
    --output-dir results/plots_no_km \
    --top-n 20
```

## Troubleshooting

### Test Fails
```bash
# Check logs
cat logs/survival_*.log

# Check intermediate files
ls -lh test_results/
head test_results/*_survival_results.txt
```

### Out of Memory
```bash
# Reduce model subset
python3 create_model_list.py --scores cidr.scores.z-scores.txt.gz --output model_list_small.txt --subset 0.01

# OR increase memory in SLURM script
# Edit run_survival_analysis.sh:
#SBATCH --mem=64G  # Increase from 32G
```

### Jobs Not Running
```bash
# Check queue
squeue -u $USER

# Check SLURM limits
sacctmgr show qos
sinfo

# Reduce resources if needed
# Edit SLURM script to use fewer CPUs/memory/time
```

### Missing Plots
```bash
# Re-run visualization only
python3 visualize_results.py \
    --meta results/meta_analysis_results.txt \
    --output-dir results/plots \
    --top-n 20 \
    --data-dir .
```

## Expected Runtime

### Test Pipeline (3% models)
- Per dataset: 5-10 minutes
- Total: 15-30 minutes

### Full Analysis (all 5110 models)
- Per dataset: 6-12 hours
- Array job (parallel): 6-12 hours total
- Sequential: 24-48 hours total
- Meta-analysis: 5-10 minutes
- Visualizations: 30-60 minutes

## Resource Requirements

### Test
- Memory: 8-16 GB per dataset
- CPUs: 1-4 per dataset
- Disk: ~100 MB

### Full Analysis
- Memory: 32-64 GB per dataset
- CPUs: 4-16 per dataset (doesn't parallelize much)
- Disk: ~1-2 GB for results

### Parallel Job
- Memory: 128 GB total (32 GB × 4 datasets)
- CPUs: 64 (16 per dataset)
- Disk: ~1-2 GB

## Getting Help

1. **Check README.md** for detailed documentation
2. **Review test results** in test_results/ directory
3. **Check error logs** in logs/ directory
4. **Verify data files** are in correct format

## Citation

If using this pipeline, cite:
- lifelines: Davidson-Pilon, C. (2019). lifelines: survival analysis in Python. JOSS, 4(40), 1317.
- PGS Catalog: https://www.pgscatalog.org/
