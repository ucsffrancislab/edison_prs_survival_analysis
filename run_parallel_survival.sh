#!/bin/bash
#SBATCH --job-name=pgs_survival_parallel
#SBATCH --output=logs/survival_parallel_%j.out
#SBATCH --error=logs/survival_parallel_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=48:00:00

# PGS Survival Analysis - Single Large Parallel Job
# This script runs all datasets in parallel using GNU parallel

# Activate your conda environment if needed
# source activate your_env

# Load modules if needed
# module load python/3.9
# module load parallel

# Set variables
DATA_DIR="."
OUTPUT_DIR="results"
MODEL_LIST="model_list.txt"
N_JOBS=4  # One per dataset

# Create output directory
mkdir -p ${OUTPUT_DIR}
mkdir -p logs

echo "=================================================="
echo "PGS Survival Analysis (Parallel)"
echo "=================================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Running on host: $(hostname)"
echo "Starting at: $(date)"
echo "=================================================="

# Function to run analysis for one dataset
run_dataset() {
    DATASET=$1
    echo "Starting analysis for ${DATASET}..."
    
    python3 survival_analysis.py \
        --dataset ${DATASET} \
        --scores ${DATA_DIR}/${DATASET}.scores.z-scores.txt.gz \
        --covariates ${DATA_DIR}/${DATASET}-covariates.tsv \
        --models ${MODEL_LIST} \
        --output ${OUTPUT_DIR}/${DATASET}_survival_results.txt \
        > logs/${DATASET}_survival.log 2>&1
    
    echo "Completed analysis for ${DATASET}"
}

export -f run_dataset
export DATA_DIR OUTPUT_DIR MODEL_LIST

# Run all datasets in parallel
parallel -j ${N_JOBS} run_dataset ::: cidr i370 onco tcga

echo "=================================================="
echo "All datasets completed"
echo "=================================================="

# Run meta-analysis
echo "Starting meta-analysis..."
python3 meta_analysis.py \
    --input ${OUTPUT_DIR}/cidr_survival_results.txt \
            ${OUTPUT_DIR}/i370_survival_results.txt \
            ${OUTPUT_DIR}/onco_survival_results.txt \
            ${OUTPUT_DIR}/tcga_survival_results.txt \
    --output ${OUTPUT_DIR}/meta_analysis_results.txt

echo "Meta-analysis completed"
echo "=================================================="

# Generate visualizations
echo "Generating visualizations..."
python3 visualize_results.py \
    --meta ${OUTPUT_DIR}/meta_analysis_results.txt \
    --output-dir ${OUTPUT_DIR}/plots \
    --top-n 20 \
    --data-dir ${DATA_DIR}

echo "Visualizations completed"
echo "=================================================="
echo "Completed at: $(date)"
echo "=================================================="
