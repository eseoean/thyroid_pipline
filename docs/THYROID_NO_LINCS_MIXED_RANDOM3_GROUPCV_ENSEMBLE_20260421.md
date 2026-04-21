# Thyroid No-LINCS Mixed Ensemble - 2026-04-21

## 목적

LINCS feature를 제거한 입력셋에서 `CrossAttention + ResidualMLP + LightGBM` 조합을 random sample 3-fold와 drug-level GroupCV 기준으로 각각 평가했다.

## 설정

- Input set: `numeric_strong_context_smiles_no_lincs`
- Removed feature family: LINCS
- Members: `CrossAttention`, `ResidualMLP`, `LightGBM`
- CV 1: random sample 3-fold OOF
- CV 2: GroupKFold 3-fold by `canonical_drug_id`
- Ensemble: Spearman-weighted average와 equal-weight average

## Random Sample 3-fold 개별 모델 성능

| cv | model | spearman | pearson | rmse | mae | r2 | ndcg_at_20 | prediction_variance |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| random3 | CrossAttention | 0.8866 | 0.9230 | 1.1201 | 0.8482 | 0.8498 | 0.9553 | 7.8462 |
| random3 | LightGBM | 0.8848 | 0.9189 | 1.1440 | 0.8685 | 0.8433 | 0.9558 | 6.5584 |
| random3 | ResidualMLP | 0.8740 | 0.9106 | 1.2051 | 0.9208 | 0.8261 | 0.9216 | 7.7506 |

## Random Sample 3-fold 앙상블 성능

| cv | ensemble | spearman | pearson | rmse | mae | r2 | ndcg_at_20 | prediction_variance |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| random3 | spearman_weighted | 0.8924 | 0.9262 | 1.0902 | 0.8299 | 0.8577 | 0.9541 | 7.2347 |
| random3 | equal_weight | 0.8923 | 0.9261 | 1.0904 | 0.8301 | 0.8576 | 0.9506 | 7.2357 |

## Random Sample 3-fold Diversity

| cv | model_a | model_b | prediction_spearman_corr | residual_pearson_corr | mean_abs_prediction_gap |
| --- | --- | --- | --- | --- | --- |
| random3 | CrossAttention | ResidualMLP | 0.9585 | 0.8311 | 0.5143 |
| random3 | CrossAttention | LightGBM | 0.9759 | 0.8788 | 0.3939 |
| random3 | ResidualMLP | LightGBM | 0.9536 | 0.7941 | 0.5612 |

## GroupCV 개별 모델 성능

| cv | model | spearman | pearson | rmse | mae | r2 | ndcg_at_20 | prediction_variance |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| groupcv | LightGBM | 0.5550 | 0.5986 | 2.3172 | 1.7596 | 0.3571 | 0.8560 | 3.3605 |
| groupcv | CrossAttention | 0.5190 | 0.5899 | 2.3489 | 1.7841 | 0.3393 | 0.8830 | 3.8941 |
| groupcv | ResidualMLP | 0.4507 | 0.4995 | 2.5577 | 1.9059 | 0.2167 | 0.9058 | 3.3935 |

## GroupCV 앙상블 성능

| cv | ensemble | spearman | pearson | rmse | mae | r2 | ndcg_at_20 | prediction_variance |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| groupcv | spearman_weighted | 0.5561 | 0.6114 | 2.2888 | 1.7239 | 0.3727 | 0.8911 | 3.0493 |
| groupcv | equal_weight | 0.5524 | 0.6079 | 2.2972 | 1.7288 | 0.3681 | 0.8934 | 3.0437 |

## GroupCV Diversity

| cv | model_a | model_b | prediction_spearman_corr | residual_pearson_corr | mean_abs_prediction_gap |
| --- | --- | --- | --- | --- | --- |
| groupcv | CrossAttention | ResidualMLP | 0.7872 | 0.8945 | 0.9160 |
| groupcv | CrossAttention | LightGBM | 0.7409 | 0.8763 | 0.8762 |
| groupcv | ResidualMLP | LightGBM | 0.6578 | 0.8414 | 1.0864 |

## GroupCV Spearman-weighted Top15 Drugs

| rank | drug_name | mixed_score | mean_pred_ln_ic50 | target_genes | PATHWAY_NAME_NORMALIZED | classification |
| --- | --- | --- | --- | --- | --- | --- |
| 1.0000 | Docetaxel | 3.1651 | -3.1651 | BCL2;NR1I2;TUBB1 | Mitosis | indication_expansion |
| 2.0000 | Dactinomycin | 2.9674 | -2.9674 | TOP1;TOP2A | Other | screened_candidate |
| 3.0000 | Vinblastine | 1.7630 | -1.7630 | JUN;TUBA1A;TUBB;TUBB2A;TUBD1;TUBE1;TUBG1 | Mitosis | screened_candidate |
| 4.0000 | Topotecan | 1.3173 | -1.3173 | TOP1;TOP1MT | DNA replication | screened_candidate |
| 5.0000 | Vinorelbine | 1.0925 | -1.0925 | BCL2;MAP4;MAPK1;TUBB | Mitosis | indication_expansion |
| 6.0000 | Paclitaxel | 1.0442 | -1.0442 | BCL2;NR1I2;TUBB1 | Mitosis | indication_expansion |
| 7.0000 | Lestaurtinib | 0.7708 | -0.7708 | FLT3;JAK2;NTRK1;NTRK2;NTRK3 | Other kinases | screened_candidate |
| 8.0000 | Teniposide | 0.4586 | -0.4586 | TOP2A;TOP2B | DNA replication | screened_candidate |
| 9.0000 | Temsirolimus | 0.2816 | -0.2816 | MTOR | PI3K/MTOR signaling | indication_expansion |
| 10.0000 | Rapamycin | -0.3974 | 0.3974 | MTOR | PI3K/MTOR signaling | indication_expansion |
| 11.0000 | Irinotecan | -0.5055 | 0.5055 | TOP1;TOP1MT | DNA replication | indication_expansion |
| 12.0000 | Tanespimycin | -0.6476 | 0.6476 | HSP90;HSP90AA1;HSP90AB1 | Protein stability and degradation | indication_expansion |
| 13.0000 | SN-38 | -0.7205 | 0.7205 | TOP1 | DNA replication | screened_candidate |
| 14.0000 | Epirubicin | -0.9144 | 0.9144 | ANTHRACYCLINE;TOP2A;TOP2B | DNA replication | indication_expansion |
| 15.0000 | Refametinib | -1.4375 | 1.4375 | MAP2K1;MAP2K2 | ERK MAPK signaling | screened_candidate |

## Weights

```json
{
  "random3": {
    "spearman_weighted": {
      "CrossAttention": 0.33513361914294987,
      "LightGBM": 0.3344798566538829,
      "ResidualMLP": 0.33038652420316716
    },
    "equal_weight": {
      "CrossAttention": 0.3333333333333333,
      "ResidualMLP": 0.3333333333333333,
      "LightGBM": 0.3333333333333333
    }
  },
  "groupcv": {
    "spearman_weighted": {
      "LightGBM": 0.3640043258500276,
      "CrossAttention": 0.34041146190230853,
      "ResidualMLP": 0.29558421224766385
    },
    "equal_weight": {
      "CrossAttention": 0.3333333333333333,
      "ResidualMLP": 0.3333333333333333,
      "LightGBM": 0.3333333333333333
    }
  }
}
```

## 산출물

- `/Users/skku_aws2_18/team4_project/pre_project/thyroid_pipline/results/mixed_no_lincs/crossattention_residualmlp_lightgbm_no_lincs/random3/individual_metrics.csv`
- `/Users/skku_aws2_18/team4_project/pre_project/thyroid_pipline/results/mixed_no_lincs/crossattention_residualmlp_lightgbm_no_lincs/random3/ensemble_metrics.csv`
- `/Users/skku_aws2_18/team4_project/pre_project/thyroid_pipline/results/mixed_no_lincs/crossattention_residualmlp_lightgbm_no_lincs/random3/ensemble_diversity.csv`
- `/Users/skku_aws2_18/team4_project/pre_project/thyroid_pipline/results/mixed_no_lincs/crossattention_residualmlp_lightgbm_no_lincs/random3/mixed_ensemble_top30_drugs.csv`
- `/Users/skku_aws2_18/team4_project/pre_project/thyroid_pipline/results/mixed_no_lincs/crossattention_residualmlp_lightgbm_no_lincs/groupcv/individual_metrics.csv`
- `/Users/skku_aws2_18/team4_project/pre_project/thyroid_pipline/results/mixed_no_lincs/crossattention_residualmlp_lightgbm_no_lincs/groupcv/ensemble_metrics.csv`
- `/Users/skku_aws2_18/team4_project/pre_project/thyroid_pipline/results/mixed_no_lincs/crossattention_residualmlp_lightgbm_no_lincs/groupcv/ensemble_diversity.csv`
- `/Users/skku_aws2_18/team4_project/pre_project/thyroid_pipline/results/mixed_no_lincs/crossattention_residualmlp_lightgbm_no_lincs/groupcv/mixed_ensemble_top30_drugs.csv`
- `reports/qc_no_lincs_mixed_random3_groupcv_ensemble_20260421.json`
