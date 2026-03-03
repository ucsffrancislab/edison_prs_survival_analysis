
# SURVIVAL ANALYSIS BUG FIX

## Problem Identified

The survival analysis was only producing results for 3 subtypes:
- all_cases
- idh_wildtype  
- idh_mutant

The other 6 molecular subtypes were silently failing, resulting in ~66% of analyses not running.

## Root Cause

When filtering data to specific molecular subtypes, covariates that are part of the filtering criteria become constant (zero variance):

**Examples:**
- In `hgg_idh_wildtype` subtype: ALL samples have `grade='HGG'` and `idh=0`
- In `lgg_idh_mutant_pq_codel` subtype: ALL samples have `grade='LGG'`, `idh=1`, and `pq=1`

After categorical encoding (e.g., `grade='HGG'` → `grade=1`), these become constant covariates.

Cox proportional hazards models **cannot** estimate coefficients for covariates with zero variance, causing ConvergenceError and analysis failure.

## Solution Implemented

Modified `survival_analysis.py` to detect and remove zero-variance covariates after categorical encoding and before fitting Cox models.

**Changes in `run_cox_model()` method (lines 122-147):**

```python
# Remove covariates with zero variance (constant after filtering)
# These cause convergence issues in Cox models
zero_variance_cols = []
for col in cox_data.columns:
    if col not in ['survdays', 'vstatus', model_name]:
        if cox_data[col].nunique() == 1:
            zero_variance_cols.append(col)

if zero_variance_cols:
    cox_data = cox_data.drop(columns=zero_variance_cols)
```

## Testing Results

**BEFORE FIX:**
- TCGA with 11 models: 33 results (3 subtypes only)
- Missing: ~66 analyses

**AFTER FIX:**
- TCGA with 5 models: 45 results (all 9 subtypes)
- CIDR with 5 models: 35 results (7 viable subtypes*)
- i370 with 5 models: 40 results (8 viable subtypes*)
- ONCO with 5 models: 44 results (8-9 viable subtypes*)

*Some subtypes have insufficient samples (N<10) in certain datasets, which is expected.

## Impact

✓ All molecular subtypes now run successfully (when sample size ≥10)
✓ No more silent failures due to zero-variance covariates
✓ ~3x increase in successful analyses
✓ Results are statistically valid (only relevant covariates included)

## Notes

- Covariates with zero variance provide no information and should not be in the model
- Removing them is the correct statistical approach
- Sample size cutoff (N≥10) remains in place as a safeguard
