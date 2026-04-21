# Thyroid Cancer Final Handoff - 2026-04-21

## 1. Recovery and Upload Scope

This handoff summarizes the recovered thyroid cancer work after the previous Codex context was restored. I rechecked the local repository, recovered session traces, Git remote, S3 layout, model output files, external validation files, ADMET summaries, and final candidate tables before preparing upload.

- Git repository: `git@github.com:eseoean/thyroid_pipline.git`
- GitHub URL: `https://github.com/eseoean/thyroid_pipline`
- Current branch: `codex/thyroid-hybrid-pipeline`
- Raw S3 source: `s3://say2-4team/thyroid_raw/`
- Result S3 target: `s3://say2-4team/20260409_eseo/20260421_thyroid/`

The existing older `phase5_final_results/FINAL_REPORT.md` and README historical summary describe the earlier random sample 3-fold run. The final basis for this upload is the newer `mixed_groupcv_pan_lincs` run, because it uses drug-level GroupCV and the final mixed model ensemble.

## 2. Final Model Basis

- Final variant: `mixed_groupcv_pan_lincs`
- Input set: `numeric + strong context + smiles + pan-cancer LINCS`
- Cross-validation: GroupCV by `canonical_drug_id`
- Ensemble members: CrossAttention, ResidualMLP, LightGBM
- Main result path: `results/mixed_groupcv/crossattention_residualmlp_lightgbm_groupcv/`
- Final candidate path: `phase5_final_results/mixed_groupcv_pan_lincs/final_comprehensive_candidates.csv`

## 3. Input Set Composition

The final input set is `numeric_strong_context_smiles_pan_lincs`.

| item | value |
| --- | ---: |
| rows | 3,387 |
| total features | 6,655 |
| numeric / pan-LINCS feature dimension | 6,559 |
| SMILES SVD features | 64 |
| strong context one-hot features | 32 |
| replaced LINCS features | 1,027 |
| pan-LINCS covered drugs | 116 |
| row-level pan-LINCS coverage | 0.4830 |
| dtype | float32 |

Earlier source-screening context:

- Raw screened response: 4,037 rows, 16 cell lines, 295 screened drugs
- Final SMILES-valid modeling pool: 3,387 rows, 16 cell lines, 243 drugs
- Removed before final modeling: 52 drugs with missing SMILES

The final array and metadata are intended for S3 rather than Git because the main model input array is large:

- `data/X_numeric_strong_context_smiles_pan_lincs.npy`
- `data/y_train.npy`
- `data/row_metadata.parquet`
- `data/thyroid_response_pairs.parquet`

## 4. Model Performance

Final GroupCV ensemble metrics:

| model / ensemble | n | Spearman | Pearson | RMSE | MAE | R2 | NDCG@20 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Spearman-weighted ensemble | 3,387 | 0.6331 | 0.6731 | 2.1392 | 1.5924 | 0.4520 | 0.8638 |
| Equal-weight ensemble | 3,387 | 0.6329 | 0.6729 | 2.1396 | 1.5925 | 0.4518 | 0.8638 |

Individual GroupCV model summary:

| model | Spearman | RMSE |
| --- | ---: | ---: |
| CrossAttention | 0.6118 | 2.2119 |
| LightGBM | 0.5928 | 2.2467 |
| ResidualMLP | 0.5846 | 2.2713 |

Interpretation: this is lower than the earlier random sample 3-fold OOF score, as expected, because drug-level GroupCV evaluates generalization to held-out drugs and is the more conservative final estimate.

## 5. External Validation

External validation was rerun on the `mixed_groupcv_pan_lincs` Top30 and Top15 candidates.

- External cohort: TCGA-THCA expression and clinical/survival files
- Expression genes available: 296
- Patient count: 572
- Mean target gene match rate: 0.9126
- Top15 target expressed: 9 / 15
- Known thyroid positive-control precision: P@5 = 0.0, P@10 = 0.0, P@20 = 0.0

Top15 model-only candidates:

1. Dactinomycin
2. Docetaxel
3. Topotecan
4. Vinblastine
5. Epirubicin
6. Mitoxantrone
7. Camptothecin
8. Lestaurtinib
9. Staurosporine
10. Piperlongumine
11. MG-132
12. Paclitaxel
13. BX795
14. Teniposide
15. Vinorelbine

Notable external-validation note: Vinorelbine had target expression percentage 0.9965 and survival p-value 0.0680 with `high_expression_longer_survival` direction.

## 6. ADMET Review

ADMET was run on the final Top15 candidates using 22 assay tables.

| item | value |
| --- | ---: |
| candidate count | 15 |
| assay count | 22 |
| Approved | 1 |
| Candidate | 8 |
| Caution | 6 |
| toxicity flag count | 6 |
| low-confidence toxic signal count | 10 |

ADMET category by drug:

- Approved: Vinblastine
- Candidate: Topotecan, Camptothecin, Lestaurtinib, Staurosporine, Piperlongumine, MG-132, BX795, Vinorelbine
- Caution: Dactinomycin, Docetaxel, Epirubicin, Mitoxantrone, Paclitaxel, Teniposide

Implementation note: `thyroid_pipeline/core.py` was adjusted so no-match AMES/DILI/hERG positives are counted as `low_confidence_toxic_signals` instead of regular toxicity flags, and `canonical_drug_id` is cast consistently before ADMET joins.

## 7. Knowledge Validation and Final Tiers

Final Tier summary:

| tier | count |
| --- | ---: |
| Tier 1 | 1 |
| Tier 2 | 3 |
| Tier 3 | 5 |
| Excluded | 6 |

Final candidate categories:

- True repurposing or exploratory candidate: 9
- Thyroid indication expansion candidate: 6

Final ranked candidates:

| final group | drugs |
| --- | --- |
| Tier 1 | Vinorelbine |
| Tier 2 | Lestaurtinib, Staurosporine, Vinblastine |
| Tier 3 | Topotecan, Camptothecin, MG-132, BX795, Piperlongumine |
| Excluded | Dactinomycin, Epirubicin, Mitoxantrone, Teniposide, Docetaxel, Paclitaxel |

The final table is `phase5_final_results/mixed_groupcv_pan_lincs/final_comprehensive_candidates.csv`.

## 8. Upload Manifest

Git upload should include:

- Pipeline code: `thyroid_pipeline/`
- Experiment scripts: `scripts/04_*.py` through `scripts/10_*.py`
- Final documentation: `docs/`
- Small result artifacts: `results/`, `reports/`, `external_validation/`, `admet/`, `knowledge_validation/`, `phase5_final_results/`

S3 upload should include the Git-uploaded artifacts plus large model-input files:

- `model_inputs/X_numeric_strong_context_smiles_pan_lincs.npy`
- `model_inputs/y_train.npy`
- `model_inputs/row_metadata.parquet`
- `model_inputs/thyroid_response_pairs.parquet`
- final feature-name JSON files where present

Excluded from S3 result bundle unless explicitly needed:

- `data/source_staging/` raw/staging cache
- bulky raw data already represented under `s3://say2-4team/thyroid_raw/`
