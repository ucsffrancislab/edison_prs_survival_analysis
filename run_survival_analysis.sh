#!/bin/bash
#SBATCH --job-name=pgs_survival
#SBATCH --output=logs/survival_%A_%a.out
#SBATCH --error=logs/survival_%A_%a.err
#SBATCH --array=0-3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00

# PGS Survival Analysis - SLURM Array Job
# This script runs survival analysis for one dataset per array task

# Activate your conda environment if needed
# source activate your_env

# Load modules if needed
# module load python/3.9

# Set variables
DATA_DIR="."
OUTPUT_DIR="results"
MODEL_LIST="model_list.txt"

# Create output directory
mkdir -p ${OUTPUT_DIR}
mkdir -p logs

# Dataset array
DATASETS=(cidr i370 onco tcga)
DATASET=${DATASETS[$SLURM_ARRAY_TASK_ID]}

echo "=================================================="
echo "PGS Survival Analysis"
echo "=================================================="
echo "Dataset: ${DATASET}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Running on host: $(hostname)"
echo "Starting at: $(date)"
echo "=================================================="

# Run survival analysis
python3 survival_analysis.py \
    --dataset ${DATASET} \
    --scores ${DATA_DIR}/${DATASET}.scores.z-scores.txt.gz \
    --covariates ${DATA_DIR}/${DATASET}-covariates.tsv \
    --models ${MODEL_LIST} \
    --output ${OUTPUT_DIR}/${DATASET}_survival_results.txt

echo "=================================================="
echo "Completed at: $(date)"
echo "=================================================="
