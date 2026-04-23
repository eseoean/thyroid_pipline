# Thyroid Step6/Step7 BRCA-Rule Validation - 2026-04-23

## Goal

Re-run the final thyroid model output with the same separated flow used for
PAAD and BRCA:

1. Step 6 external cohort validation
2. Step 7 ADMET reranking with the exact BRCA rules

## Input

- Final thyroid model basis: `mixed_groupcv_pan_lincs`
- Source Top30:
  - `results/mixed_groupcv/crossattention_residualmlp_lightgbm_groupcv/mixed_groupcv_ensemble_top30_drugs.csv`
- External validation source reused from the existing thyroid run:
  - `external_validation/mixed_groupcv_pan_lincs/thyroid_target_expression.csv`
  - `external_validation/mixed_groupcv_pan_lincs/thyroid_survival_validation.csv`

## Step 6

New output folder:

- `external_validation/mixed_groupcv_pan_lincs_step6_external_cohort/`

Primary rule:

- `external_cohort_validated = target_expressed in TCGA-THCA`

Key counts:

- External target-expression support: `19/30`
- Survival significant: `0/30`
- Known thyroid drugs in Top30: `2/30`

Top15 for Step 7 now comes from Step 6 reranking rather than the previous
model-only Top15.

## Step 7

New output folder:

- `admet/mixed_groupcv_pan_lincs_step7_admet_brca_rules/`

BRCA rules copied unchanged:

- same 22 ADMET assays
- same toxicity penalties
- same known-approved override set
- `Approved` if known-approved, else `Candidate` if safety score >= 4, else `Caution`

Category counts:

- `Approved = 5`
- `Candidate = 8`
- `Caution = 2`

## Main files

- Step 6 report:
  - `reports/mixed_groupcv_pan_lincs_step6_external_cohort/THYROID_STEP6_EXTERNAL_COHORT_VALIDATION_20260423.md`
- Step 7 report:
  - `reports/mixed_groupcv_pan_lincs_step7_admet_brca_rules/THYROID_STEP7_ADMET_BRCA_RULES_20260423.md`
- Combined run summary:
  - `reports/thyroid_step6_step7_brca_run_summary_20260423.json`
- Implementation:
  - `scripts/21_run_thyroid_step6_external_step7_brca_rules.py`
