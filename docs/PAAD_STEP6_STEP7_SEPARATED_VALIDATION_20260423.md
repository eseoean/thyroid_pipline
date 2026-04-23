# PAAD Step6/Step7 Separated Validation

Date: 2026-04-23

## BRCA reference flow checked

The teammate BRCA flow separates the steps:

1. Step 6 METABRIC external validation
   - Input: model Step 5 Top30 drugs.
   - Method A: target expression in METABRIC BRCA patients.
   - Method B: survival stratification by target expression.
   - Method C: known BRCA drug precision at K.
   - Output: `models/metabric_results*/top15_validated.csv` and `step6_metabric_results.json`.

2. Step 7 ADMET gate
   - Input: Step 6 `top15_validated.csv`.
   - Looks up 22 ADMET assays after external validation.
   - Output: `models/admet_results*/final_drug_candidates.csv` and `step7_admet_results.json`.

This means ADMET is not part of the METABRIC target-expression count such as `29/30`; it is a downstream safety annotation/filtering step.

## PAAD implementation

Implemented the same separation in:

- `scripts/19_run_paad_step6_external_step7_admet.py`

Step 6 uses the existing PAAD Top50 external source table from:

- `external_validation/paad/groupcv4_drug_v2_sources/top50_external_validation_v2.csv`

but excludes ADMET/SIDER from Step 6 counting.

## PAAD Step 6 outputs

Local:

- `external_validation/paad/groupcv4_drug_step6_external_cohort/all_30_scores.csv`
- `external_validation/paad/groupcv4_drug_step6_external_cohort/top15_validated.csv`
- `external_validation/paad/groupcv4_drug_step6_external_cohort/step6_external_cohort_results.json`
- `reports/paad/groupcv4_drug_step6_external_cohort/PAAD_STEP6_EXTERNAL_COHORT_VALIDATION_20260423.md`

S3:

- `s3://say2-4team/PAAD_outputs/external_validation/paad/groupcv4_drug_step6_external_cohort/`
- `s3://say2-4team/PAAD_outputs/reports/paad/groupcv4_drug_step6_external_cohort/`

Key Step 6 counts:

| Metric | Count |
|---|---:|
| External PAAD cohort/protein validated, GEO or CPTAC | 27/30 |
| GEO external validation, GSE62452 or GSE71729 | 27/30 |
| Both GEO cohorts validated | 17/30 |
| Any disease reference including TCGA-GTEx | 28/30 |
| TCGA-GTEx target up | 11/30 |
| GSE62452 target up | 25/30 |
| GSE71729 target up | 19/30 |
| CPTAC protein up | 19/30 |
| Survival significant, p < 0.05 | 2/30 |
| PRISM evidence | 27/30 |

## PAAD Step 7 outputs

Local:

- `admet/paad/groupcv4_drug_step7_admet_after_external/final_drug_candidates.csv`
- `admet/paad/groupcv4_drug_step7_admet_after_external/step7_admet_ranked_candidates.csv`
- `admet/paad/groupcv4_drug_step7_admet_after_external/step7_admet_results.json`
- `reports/paad/groupcv4_drug_step7_admet_after_external/PAAD_STEP7_ADMET_AFTER_EXTERNAL_20260423.md`

S3:

- `s3://say2-4team/PAAD_outputs/admet/paad/groupcv4_drug_step7_admet_after_external/`
- `s3://say2-4team/PAAD_outputs/reports/paad/groupcv4_drug_step7_admet_after_external/`

Step 7 category counts:

| ADMET category | Count |
|---|---:|
| Candidate | 10 |
| Caution | 5 |

Top Step 7 ranked candidates:

| Rank | Drug | ADMET category |
|---:|---|---|
| 1 | Topotecan | Candidate |
| 2 | AZD5438 | Candidate |
| 3 | BI-2536 | Candidate |
| 4 | SN-38 | Candidate |
| 5 | Camptothecin | Candidate |
| 6 | Vinorelbine | Candidate |
| 7 | MK-1775 | Candidate |
| 8 | MG-132 | Candidate |
| 9 | OSI-027 | Candidate |
| 10 | Staurosporine | Candidate |

