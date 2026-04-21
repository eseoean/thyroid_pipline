# Thyroid Additional DL Models - 2026-04-21

## 목적

기존 thyroid 파이프라인은 ML tree 계열과 sklearn `FlatMLP` 중심이었다. BRCA 파이프라인에서 사용했던 DL 계열을 thyroid primary input set에 추가 적용해, DL 단독 성능과 DL-only ensemble/diversity를 확인했다.

## 실행 설정

- Input set: `numeric_strong_context_smiles_no_lincs`
- CV: random sample 3-fold OOF
- 추가 학습 모델: `ResidualMLP`, `WideDeep`, `CrossAttention`
- 비교용 포함: 기존 pan-cancer LINCS `FlatMLP_existing` OOF
- Device: PyTorch auto device

## 추가 DL 모델 성능

| input_set | model | spearman | pearson | rmse | mae | r2 | ndcg_at_20 | train_oof_spearman_gap | elapsed_sec |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| numeric_strong_context_smiles_no_lincs | CrossAttention | 0.8866 | 0.9230 | 1.1201 | 0.8482 | 0.8498 | 0.9553 | 0.0450 | 10.0059 |
| numeric_strong_context_smiles_no_lincs | ResidualMLP | 0.8740 | 0.9106 | 1.2051 | 0.9208 | 0.8261 | 0.9216 | 0.0881 | 24.2079 |

## DL-only Weighted Ensemble

- Ensemble Spearman: `0.8892`
- Ensemble RMSE: `1.1121`
- Ensemble R2: `0.8519`
- Weights: `{"CrossAttention": 0.503566454740134, "ResidualMLP": 0.4964335452598661}`

## DL Ensemble Top 30 중 상위 15

| rank | drug_name | dl_ensemble_score | mean_pred_ln_ic50 | target_genes | PATHWAY_NAME_NORMALIZED | classification |
| --- | --- | --- | --- | --- | --- | --- |
| 1.0000 | Romidepsin | 5.6158 | -5.6158 | HDAC1;HDAC2;HDAC3;HDAC4;HDAC6;HDAC8 | Chromatin histone acetylation | indication_expansion |
| 2.0000 | Bortezomib | 5.3360 | -5.3360 | PROTEASOME;PRSS1;PSMB1;PSMB5 | Protein stability and degradation | indication_expansion |
| 3.0000 | Sepantronium bromide | 5.0775 | -5.0775 | BIRC5 | Apoptosis regulation | screened_candidate |
| 4.0000 | Dactinomycin | 4.0338 | -4.0338 | TOP1;TOP2A | Other | screened_candidate |
| 5.0000 | Docetaxel | 3.7833 | -3.7833 | BCL2;NR1I2;TUBB1 | Mitosis | indication_expansion |
| 6.0000 | Vinorelbine | 3.3968 | -3.3968 | BCL2;MAP4;MAPK1;TUBB | Mitosis | indication_expansion |
| 7.0000 | SN-38 | 3.1021 | -3.1021 | TOP1 | DNA replication | screened_candidate |
| 8.0000 | Daporinad | 3.0904 | -3.0904 | NAMPT | Metabolism | screened_candidate |
| 9.0000 | Staurosporine | 3.0897 | -3.0897 | CDK2;CHRM1;CSK;GSK3B;IKBKB;ITK;LCK;MAPKAPK2;PDPK1;PIK3CG;PIM1;PRKCA;PRKCQ;SYK;ZAP70 | RTK signaling | screened_candidate |
| 10.0000 | Vinblastine | 2.7853 | -2.7853 | JUN;TUBA1A;TUBB;TUBB2A;TUBD1;TUBE1;TUBG1 | Mitosis | screened_candidate |
| 11.0000 | Paclitaxel | 2.3480 | -2.3480 | BCL2;NR1I2;TUBB1 | Mitosis | indication_expansion |
| 12.0000 | Camptothecin | 2.2512 | -2.2512 | TOP1 | DNA replication | indication_expansion |
| 13.0000 | Dinaciclib | 2.1458 | -2.1458 | CDK1;CDK2;CDK5;CDK9 | Cell cycle | screened_candidate |
| 14.0000 | Luminespib | 1.9525 | -1.9525 | HSP90 | Protein stability and degradation | screened_candidate |
| 15.0000 | MG-132 | 1.4842 | -1.4842 | CAPN1;PROTEASOME | Protein stability and degradation | screened_candidate |

## DL Diversity

| model_a | model_b | prediction_spearman_corr | residual_pearson_corr | mean_abs_prediction_gap |
| --- | --- | --- | --- | --- |
| CrossAttention | ResidualMLP | 0.9585 | 0.8311 | 0.5143 |

## 산출물

- `results/additional_dl/numeric_strong_context_smiles_no_lincs/additional_dl_metrics_summary.csv`
- `results/additional_dl/numeric_strong_context_smiles_no_lincs/additional_dl_fold_metrics.csv`
- `results/additional_dl/numeric_strong_context_smiles_no_lincs/oof/`
- `results/additional_dl/numeric_strong_context_smiles_no_lincs/ensemble/additional_dl_ensemble_results.json`
- `results/additional_dl/numeric_strong_context_smiles_no_lincs/ensemble/additional_dl_ensemble_top30_drugs.csv`
- `results/additional_dl/numeric_strong_context_smiles_no_lincs/ensemble/additional_dl_ensemble_diversity.csv`
