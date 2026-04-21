# Thyroid No-LINCS Ablation Experiment - 2026-04-21

## 목적

갑상선암 파이프라인에서 BRCA/MCF7 기반 LINCS feature가 질환 특이적이지 않다는 리스크가 있어, LINCS 관련 컬럼을 모두 제거한 입력셋으로 random sample 3-fold 학습과 GroupCV stress test를 다시 수행했다.

## 제거 기준

feature name에 `lincs`가 포함된 모든 컬럼을 제거했다. 여기에는 `drug__lincs__*` signature feature와 `drug__has_lincs_signature`, `drug__lincs_signature_source_direct`, `drug__lincs_signature_source_recovered` availability/source flag가 포함된다.

## 입력셋 차원 변화

| input_set | source_input_set | rows | source_features | removed_lincs_features | remaining_features | array_path | feature_names_path |
| --- | --- | --- | --- | --- | --- | --- | --- |
| numeric_no_lincs | numeric | 3387.0000 | 6559.0000 | 1027.0000 | 5532.0000 | /Users/skku_aws2_18/team4_project/pre_project/thyroid_pipline/data/X_numeric_no_lincs.npy | /Users/skku_aws2_18/team4_project/pre_project/thyroid_pipline/data/numeric_no_lincs_feature_names.json |
| numeric_smiles_no_lincs | numeric_smiles | 3387.0000 | 6623.0000 | 1027.0000 | 5596.0000 | /Users/skku_aws2_18/team4_project/pre_project/thyroid_pipline/data/X_numeric_smiles_no_lincs.npy | /Users/skku_aws2_18/team4_project/pre_project/thyroid_pipline/data/numeric_smiles_no_lincs_feature_names.json |
| numeric_strong_context_smiles_no_lincs | numeric_strong_context_smiles | 3387.0000 | 6655.0000 | 1027.0000 | 5628.0000 | /Users/skku_aws2_18/team4_project/pre_project/thyroid_pipline/data/X_numeric_strong_context_smiles_no_lincs.npy | /Users/skku_aws2_18/team4_project/pre_project/thyroid_pipline/data/numeric_strong_context_smiles_no_lincs_feature_names.json |

## Random Sample 3-fold Best Model 비교

| input_set | variant | model | spearman | rmse | r2 | ndcg_at_20 |
| --- | --- | --- | --- | --- | --- | --- |
| numeric | baseline_with_lincs | LightGBM | 0.8776 | 1.1587 | 0.8392 | 0.9692 |
| numeric | no_lincs | LightGBM | 0.8762 | 1.1906 | 0.8303 | 0.9339 |
| numeric_smiles | baseline_with_lincs | LightGBM | 0.8833 | 1.1353 | 0.8457 | 0.9654 |
| numeric_smiles | no_lincs | LightGBM | 0.8845 | 1.1493 | 0.8418 | 0.9490 |
| numeric_strong_context_smiles | baseline_with_lincs | LightGBM | 0.8835 | 1.1264 | 0.8481 | 0.9689 |
| numeric_strong_context_smiles | no_lincs | LightGBM | 0.8848 | 1.1440 | 0.8433 | 0.9558 |

## GroupCV Stress Test 비교

| variant | input_set | model | n | spearman | pearson | rmse | mae | r2 | ndcg_at_20 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| baseline_with_lincs | numeric_strong_context_smiles | ExtraTrees | 3387.0000 | 0.5112 | 0.4903 | 2.6134 | 1.9301 | 0.1822 | 0.7222 |
| no_lincs | numeric_strong_context_smiles_no_lincs | ExtraTrees | 3387.0000 | 0.4180 | 0.4391 | 2.7657 | 2.0235 | 0.0841 | 0.8615 |

## No-LINCS Ensemble

- Primary input set: `numeric_strong_context_smiles_no_lincs`
- Best single model: `LightGBM`
- Ensemble Spearman: `0.8828`
- Ensemble RMSE: `1.1956`
- Ensemble R2: `0.8288`

## 산출물

- `data/X_numeric_no_lincs.npy`
- `data/X_numeric_smiles_no_lincs.npy`
- `data/X_numeric_strong_context_smiles_no_lincs.npy`
- `results/no_lincs/random3/`
- `results/no_lincs/groupcv_stress_test/`
- `results/no_lincs/ensemble/`
- `reports/qc_no_lincs_training_20260421.json`
- `reports/qc_no_lincs_ensemble_20260421.json`
