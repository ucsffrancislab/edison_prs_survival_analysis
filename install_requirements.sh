#!/bin/bash

# Installation Script for PGS Survival Analysis
# Run this script before using the pipeline

set -e

echo "=================================================="
echo "Installing Required Packages"
echo "=================================================="

# Check Python version
PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}')
echo "Python version: $PYTHON_VERSION"

# Install packages
echo ""
echo "Installing Python packages..."
pip install --upgrade pip
pip install pandas numpy scipy matplotlib seaborn lifelines

echo ""
echo "=================================================="
echo "Installation Complete"
echo "=================================================="
echo ""
echo "Packages installed:"
pip list | grep -E "pandas|numpy|scipy|matplotlib|seaborn|lifelines"

echo ""
echo "=================================================="
echo "Next Steps:"
echo "=================================================="
echo "1. Run test pipeline:"
echo "   ./test_pipeline.sh"
echo ""
echo "2. If test succeeds, run full analysis:"
echo "   sbatch archive/run_survival_analysis.sh  (array job)"
echo "   OR"
echo "   sbatch run_pipeline.sh  (parallel job)"
echo ""
