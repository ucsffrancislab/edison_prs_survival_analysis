#!/bin/bash
#SBATCH --job-name=pgs_survival_parallel
#SBATCH --output=survival_parallel_%j.out
#SBATCH --error=survival_parallel_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=490G
#SBATCH --time=14-0

# PGS Survival Analysis - Single Large Parallel Job
# This script runs all datasets in parallel using GNU parallel

# ── Locate the directory this script lives in ────────────────────────────────
# Under SLURM, sbatch copies the script to a spool directory so dirname "$0"
# points to the wrong place. Use scontrol to recover the original path.
# Outside SLURM (interactive/local), dirname "$0" works fine.
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    PIPELINE_DIR=$(dirname "$(scontrol show job "$SLURM_JOB_ID" \
        | awk '/Command=/{sub(/.*Command=/, ""); print $1}')")
else
    PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"
fi

# Activate your conda environment if needed
# source activate your_env

# Load modules if needed
# module load python/3.9
# module load parallel

# Set variables
MODEL_LIST="model_list.txt"
N_JOBS=4  # One per dataset
DATA_DIR="input"
OUTPUT_DIR="results"
#PYTHON_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --datadir)
            DATA_DIR="$2"
            shift 2
            ;;
        --outdir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --models)
            MODEL_LIST="$2"
            shift 2
            ;;
#        *)
#            PYTHON_ARGS+=("$1")
#            shift
#            ;;
    esac
done

if [[ -z "$OUTPUT_DIR" ]]; then
    echo "ERROR: --outdir is required" >&2
    echo "Usage: sbatch run_pipeline.sh --outdir <directory> [pipeline options]" >&2
    exit 1
fi


# ── Create output directory ───────────────────────────────────────────────────
mkdir -p ${OUTPUT_DIR}
mkdir -p logs

echo "=================================================="
echo "PGS Survival Analysis (Parallel)"
echo "=================================================="
echo "Job ID:       ${SLURM_JOB_ID:-local}"
echo "CPUs:         ${SLURM_CPUS_PER_TASK:-$(nproc)}"
echo "Running on:   $(hostname)"
echo "Pipeline dir: $PIPELINE_DIR"
echo "Starting at:  $(date)"
echo "=================================================="

# Function to run analysis for one dataset
run_dataset() {
    DATASET=$1
    echo "Starting analysis for ${DATASET}..."
    
    python3 "$PIPELINE_DIR/survival_analysis.py" \
        --dataset ${DATASET} \
        --scores ${DATA_DIR}/${DATASET}.scores.z-scores.txt.gz \
        --covariates ${DATA_DIR}/${DATASET}-covariates.csv \
        --models ${MODEL_LIST} \
        --output ${OUTPUT_DIR}/${DATASET}_survival_results.txt \
        > logs/${DATASET}_survival.log 2>&1
    
    echo "Completed analysis for ${DATASET}"
}

export -f run_dataset
export DATA_DIR OUTPUT_DIR MODEL_LIST PIPELINE_DIR

# Run all datasets in parallel
parallel -j ${N_JOBS} run_dataset ::: cidr i370 onco tcga

echo "=================================================="
echo "All datasets completed"
echo "=================================================="

# Run meta-analysis
echo "Starting meta-analysis..."
python3 "$PIPELINE_DIR/meta_analysis.py" \
    --input ${OUTPUT_DIR}/cidr_survival_results.txt \
            ${OUTPUT_DIR}/i370_survival_results.txt \
            ${OUTPUT_DIR}/onco_survival_results.txt \
            ${OUTPUT_DIR}/tcga_survival_results.txt \
    --output ${OUTPUT_DIR}/meta_analysis_results.txt

echo "Meta-analysis completed"
echo "=================================================="

# Generate visualizations
echo "Generating visualizations..."
python3 "$PIPELINE_DIR/visualize_results.py" \
    --meta ${OUTPUT_DIR}/meta_analysis_results.txt \
    --output-dir ${OUTPUT_DIR}/plots \
    --top-n 100 \
    --data-dir ${DATA_DIR}

echo "Visualizations completed"
echo "=================================================="
echo "Completed at: $(date)"
echo "=================================================="
