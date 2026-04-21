# Thyroid Final Summary - mixed_groupcv_pan_lincs

This directory contains the final thyroid cancer candidate outputs for the 2026-04-21 conservative run.

## Final Basis

- Variant: `mixed_groupcv_pan_lincs`
- Input: `numeric + strong context + smiles + pan-cancer LINCS`
- CV: GroupCV by `canonical_drug_id`
- Ensemble: CrossAttention + ResidualMLP + LightGBM
- Final table: `final_comprehensive_candidates.csv`

## Performance

| ensemble | n | Spearman | Pearson | RMSE | MAE | R2 | NDCG@20 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Spearman-weighted | 3,387 | 0.6331 | 0.6731 | 2.1392 | 1.5924 | 0.4520 | 0.8638 |
| Equal-weight | 3,387 | 0.6329 | 0.6729 | 2.1396 | 1.5925 | 0.4518 | 0.8638 |

## Input Composition

- Rows: 3,387
- Features: 6,655
- Numeric / pan-LINCS features: 6,559
- SMILES SVD features: 64
- Strong context features: 32
- Replaced LINCS features: 1,027
- Pan-LINCS drugs: 116
- Row-level pan-LINCS coverage: 0.4830

## Validation Summary

- TCGA-THCA patient count: 572
- Target gene match rate mean: 0.9126
- Top15 target expressed: 9 / 15
- ADMET categories: Approved 1, Candidate 8, Caution 6
- Final tiers: Tier 1 = 1, Tier 2 = 3, Tier 3 = 5, Excluded = 6

## Final Tiers

| tier | drugs |
| --- | --- |
| Tier 1 | Vinorelbine |
| Tier 2 | Lestaurtinib, Staurosporine, Vinblastine |
| Tier 3 | Topotecan, Camptothecin, MG-132, BX795, Piperlongumine |
| Excluded | Dactinomycin, Epirubicin, Mitoxantrone, Teniposide, Docetaxel, Paclitaxel |
