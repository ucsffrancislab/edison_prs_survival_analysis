#!/bin/bash

# Test Pipeline for PGS Survival Analysis
# This script tests the analysis on a small subset before running the full analysis

set -e  # Exit on error

# ── Locate the directory this script lives in ────────────────────────────────
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    PIPELINE_DIR=$(dirname "$(scontrol show job "$SLURM_JOB_ID" \
        | awk '/Command=/{sub(/.*Command=/, ""); print $1}')")
else
    PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"
fi

echo "=================================================="
echo "PGS Survival Analysis - Test Pipeline"
echo "=================================================="
echo "Pipeline dir: $PIPELINE_DIR"
echo "Starting at: $(date)"
echo ""

# Create directories
echo "Creating directories..."
mkdir -p test_results
mkdir -p test_results/plots
mkdir -p logs

# Step 1: Create test model list (3% of models)
echo ""
echo "Step 1: Creating test model list (3% of models)..."
python3 "$PIPELINE_DIR/create_model_list.py" \
    --scores cidr.scores.z-scores.txt.gz \
    --output test_model_list.txt \
    --subset 0.03

# Count models
N_MODELS=$(wc -l < test_model_list.txt)
echo "Test will analyze ${N_MODELS} models"

# Step 2: Run survival analysis for each dataset
echo ""
echo "Step 2: Running survival analysis for each dataset..."

for DATASET in cidr i370 onco tcga; do
    echo ""
    echo "  Analyzing ${DATASET}..."
    
    python3 "$PIPELINE_DIR/survival_analysis.py" \
        --dataset ${DATASET} \
        --scores ${DATASET}.scores.z-scores.txt.gz \
        --covariates ${DATASET}-covariates.csv \
        --models test_model_list.txt \
        --output test_results/${DATASET}_survival_results.txt
    
    # Check if output was created
    if [ -f test_results/${DATASET}_survival_results.txt ]; then
        N_RESULTS=$(wc -l < test_results/${DATASET}_survival_results.txt)
        echo "    ✓ ${DATASET} complete: ${N_RESULTS} results"
    else
        echo "    ✗ ${DATASET} failed: no output file"
        exit 1
    fi
done

# Step 3: Run meta-analysis
echo ""
echo "Step 3: Running meta-analysis..."

python3 "$PIPELINE_DIR/meta_analysis.py" \
    --input test_results/cidr_survival_results.txt \
            test_results/i370_survival_results.txt \
            test_results/onco_survival_results.txt \
            test_results/tcga_survival_results.txt \
    --output test_results/meta_analysis_results.txt

if [ -f test_results/meta_analysis_results.txt ]; then
    N_META=$(wc -l < test_results/meta_analysis_results.txt)
    echo "  ✓ Meta-analysis complete: ${N_META} results"
else
    echo "  ✗ Meta-analysis failed"
    exit 1
fi

# Step 4: Generate visualizations
echo ""
echo "Step 4: Generating visualizations..."

python3 "$PIPELINE_DIR/visualize_results.py" \
    --meta test_results/meta_analysis_results.txt \
    --output-dir test_results/plots \
    --top-n 10 \
    --data-dir .

# Count plots
N_PLOTS=$(find test_results/plots -name "*.png" | wc -l)
echo "  ✓ Visualizations complete: ${N_PLOTS} plots created"

# Step 5: Summary
echo ""
echo "=================================================="
echo "Test Pipeline Summary"
echo "=================================================="
echo "Models analyzed: ${N_MODELS}"
echo "Datasets: 4 (cidr, i370, onco, tcga)"
echo "Meta-analysis results: ${N_META}"
echo "Plots created: ${N_PLOTS}"
echo ""
echo "Test results saved to: test_results/"
echo ""

# Show top 5 results
echo "Top 5 meta-analysis results:"
echo "----------------------------"
head -n 6 test_results/meta_analysis_results.txt | cut -f1-6

echo ""
echo "=================================================="
echo "Test completed successfully at: $(date)"
echo "=================================================="
echo ""
echo "Next steps:"
echo "1. Review test results in test_results/ directory"
echo "2. If satisfied, create full model list:"
echo "   python3 $PIPELINE_DIR/create_model_list.py --scores cidr.scores.z-scores.txt.gz --output model_list.txt"
echo ""
echo "3. Submit full analysis to SLURM:"
echo "   sbatch $PIPELINE_DIR/run_survival_analysis.sh  (for array job)"
echo "   OR"
echo "   sbatch $PIPELINE_DIR/run_parallel_survival.sh  (for single parallel job)"
echo ""
