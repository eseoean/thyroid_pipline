# Thyroid Mixed GroupCV Ensemble - 2026-04-21

## 목적

`CrossAttention + ResidualMLP + LightGBM` 조합을 drug-level GroupCV 기준으로 평가했다. CrossAttention/ResidualMLP는 기존 추가 DL GroupCV OOF를 재사용했고, LightGBM은 같은 `canonical_drug_id` GroupKFold split으로 새로 학습했다.

## 설정

- Input set: `numeric_strong_context_smiles_pan_lincs`
- CV: GroupKFold 3-fold by `canonical_drug_id`
- Members: `CrossAttention`, `ResidualMLP`, `LightGBM`
- Ensemble: Spearman-weighted average와 equal-weight average를 함께 산출

## 개별 모델 GroupCV 성능

| model | spearman | pearson | rmse | mae | r2 | ndcg_at_20 | prediction_variance |
| --- | --- | --- | --- | --- | --- | --- | --- |
| CrossAttention | 0.6118 | 0.6507 | 2.2119 | 1.6545 | 0.4142 | 0.8728 | 4.3311 |
| LightGBM | 0.5928 | 0.6313 | 2.2467 | 1.6807 | 0.3956 | 0.8386 | 3.7106 |
| ResidualMLP | 0.5846 | 0.6328 | 2.2713 | 1.6846 | 0.3823 | 0.8696 | 4.9111 |

## 앙상블 GroupCV 성능

| ensemble | spearman | pearson | rmse | mae | r2 | ndcg_at_20 | prediction_variance |
| --- | --- | --- | --- | --- | --- | --- | --- |
| spearman_weighted | 0.6331 | 0.6731 | 2.1392 | 1.5924 | 0.4520 | 0.8638 | 3.8702 |
| equal_weight | 0.6329 | 0.6729 | 2.1396 | 1.5925 | 0.4518 | 0.8638 | 3.8725 |

## Diversity

| model_a | model_b | prediction_spearman_corr | residual_pearson_corr | mean_abs_prediction_gap |
| --- | --- | --- | --- | --- |
| CrossAttention | ResidualMLP | 0.8365 | 0.8920 | 0.8364 |
| CrossAttention | LightGBM | 0.7991 | 0.8527 | 0.9171 |
| ResidualMLP | LightGBM | 0.7906 | 0.8565 | 0.8974 |

## Spearman-weighted Top15 Drugs

| rank | drug_name | mixed_groupcv_score | mean_pred_ln_ic50 | target_genes | PATHWAY_NAME_NORMALIZED | classification |
| --- | --- | --- | --- | --- | --- | --- |
| 1.0000 | Dactinomycin | 4.3256 | -4.3256 | TOP1;TOP2A | Other | screened_candidate |
| 2.0000 | Docetaxel | 3.3417 | -3.3417 | BCL2;NR1I2;TUBB1 | Mitosis | indication_expansion |
| 3.0000 | Topotecan | 2.8960 | -2.8960 | TOP1;TOP1MT | DNA replication | screened_candidate |
| 4.0000 | Vinblastine | 2.3322 | -2.3322 | JUN;TUBA1A;TUBB;TUBB2A;TUBD1;TUBE1;TUBG1 | Mitosis | screened_candidate |
| 5.0000 | Epirubicin | 1.8480 | -1.8480 | ANTHRACYCLINE;TOP2A;TOP2B | DNA replication | indication_expansion |
| 6.0000 | Mitoxantrone | 1.6478 | -1.6478 | TOP2;TOP2A;TOP2B | DNA replication | indication_expansion |
| 7.0000 | Camptothecin | 0.8951 | -0.8951 | TOP1 | DNA replication | indication_expansion |
| 8.0000 | Lestaurtinib | 0.5204 | -0.5204 | FLT3;JAK2;NTRK1;NTRK2;NTRK3 | Other kinases | screened_candidate |
| 9.0000 | Staurosporine | 0.4478 | -0.4478 | CDK2;CHRM1;CSK;GSK3B;IKBKB;ITK;LCK;MAPKAPK2;PDPK1;PIK3CG;PIM1;PRKCA;PRKCQ;SYK;ZAP70 | RTK signaling | screened_candidate |
| 10.0000 | Piperlongumine | 0.2498 | -0.2498 |  | Other | screened_candidate |
| 11.0000 | MG-132 | 0.1759 | -0.1759 | CAPN1;PROTEASOME | Protein stability and degradation | screened_candidate |
| 12.0000 | Paclitaxel | 0.1459 | -0.1459 | BCL2;NR1I2;TUBB1 | Mitosis | indication_expansion |
| 13.0000 | BX795 | 0.0037 | -0.0037 | AURKB;AURKC;IKK;TBK1 | Other kinases | screened_candidate |
| 14.0000 | Teniposide | -0.0667 | 0.0667 | TOP2A;TOP2B | DNA replication | screened_candidate |
| 15.0000 | Vinorelbine | -0.1815 | 0.1815 | BCL2;MAP4;MAPK1;TUBB | Mitosis | indication_expansion |

## Weights

```json
{
  "spearman_weighted": {
    "CrossAttention": 0.3419292286029331,
    "LightGBM": 0.3313320535456792,
    "ResidualMLP": 0.3267387178513877
  },
  "equal_weight": {
    "CrossAttention": 0.3333333333333333,
    "ResidualMLP": 0.3333333333333333,
    "LightGBM": 0.3333333333333333
  }
}
```

## 산출물

- `results/mixed_groupcv/crossattention_residualmlp_lightgbm_groupcv/individual_groupcv_metrics.csv`
- `results/mixed_groupcv/crossattention_residualmlp_lightgbm_groupcv/ensemble_groupcv_metrics.csv`
- `results/mixed_groupcv/crossattention_residualmlp_lightgbm_groupcv/ensemble_diversity.csv`
- `results/mixed_groupcv/crossattention_residualmlp_lightgbm_groupcv/mixed_groupcv_ensemble_top30_drugs.csv`
- `results/mixed_groupcv/crossattention_residualmlp_lightgbm_groupcv/mixed_groupcv_ensemble_results.json`
- `reports/qc_mixed_groupcv_crossattention_residualmlp_lightgbm_20260421.json`
