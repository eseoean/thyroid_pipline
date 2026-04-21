# Thyroid Pan-Cancer LINCS Experiment - 2026-04-21

## 목적

MCF7-only LINCS feature의 유방암 세포주 편향을 줄이기 위해, LINCS GSE92742의 tumor sample-type cancer cell line 전체를 사용한 pan-cancer drug perturbation signature로 LINCS block을 교체했다.

## Pan-Cancer LINCS 정의

- Source: `GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz`
- Cell filter: `sample_type == tumor`
- Perturbation filter: `pert_type == trt_cp`
- Drug matching: thyroid screened drug의 normalized name 우선, canonical SMILES 보조
- Feature replacement: 기존 `drug__lincs__*` 1024개 gene feature와 LINCS availability flag 3개를 pan-cancer 값으로 교체

## Signature QC

- Mapped drugs: `116`
- Mapped cells: `53`
- Mapped unique signatures: `8430`
- Matched GCTX gene features: `1024`

## 입력셋

| input_set | source_input_set | rows | features | replaced_lincs_features | pan_lincs_drugs | row_level_has_pan_lincs_ratio | array_path | feature_names_path |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| numeric_pan_lincs | numeric | 3387.0000 | 6559.0000 | 1027.0000 | 116.0000 | 0.4830 | /Users/skku_aws2_18/team4_project/pre_project/thyroid_pipline/data/X_numeric_pan_lincs.npy | /Users/skku_aws2_18/team4_project/pre_project/thyroid_pipline/data/numeric_pan_lincs_feature_names.json |
| numeric_smiles_pan_lincs | numeric_smiles | 3387.0000 | 6623.0000 | 1027.0000 | 116.0000 | 0.4830 | /Users/skku_aws2_18/team4_project/pre_project/thyroid_pipline/data/X_numeric_smiles_pan_lincs.npy | /Users/skku_aws2_18/team4_project/pre_project/thyroid_pipline/data/numeric_smiles_pan_lincs_feature_names.json |
| numeric_strong_context_smiles_pan_lincs | numeric_strong_context_smiles | 3387.0000 | 6655.0000 | 1027.0000 | 116.0000 | 0.4830 | /Users/skku_aws2_18/team4_project/pre_project/thyroid_pipline/data/X_numeric_strong_context_smiles_pan_lincs.npy | /Users/skku_aws2_18/team4_project/pre_project/thyroid_pipline/data/numeric_strong_context_smiles_pan_lincs_feature_names.json |

## Random Sample 3-fold Best Model 비교

| source_input_set | variant | model | spearman | rmse | r2 | ndcg_at_20 |
| --- | --- | --- | --- | --- | --- | --- |
| numeric | baseline_mcf7_lincs | LightGBM | 0.8776 | 1.1587 | 0.8392 | 0.9692 |
| numeric | no_lincs | LightGBM | 0.8762 | 1.1906 | 0.8303 | 0.9339 |
| numeric | pancancer_lincs | LightGBM | 0.8761 | 1.1769 | 0.8342 | 0.9345 |
| numeric_smiles | baseline_mcf7_lincs | LightGBM | 0.8833 | 1.1353 | 0.8457 | 0.9654 |
| numeric_smiles | no_lincs | LightGBM | 0.8845 | 1.1493 | 0.8418 | 0.9490 |
| numeric_smiles | pancancer_lincs | LightGBM | 0.8803 | 1.1600 | 0.8389 | 0.9350 |
| numeric_strong_context_smiles | baseline_mcf7_lincs | LightGBM | 0.8835 | 1.1264 | 0.8481 | 0.9689 |
| numeric_strong_context_smiles | no_lincs | LightGBM | 0.8848 | 1.1440 | 0.8433 | 0.9558 |
| numeric_strong_context_smiles | pancancer_lincs | LightGBM | 0.8801 | 1.1584 | 0.8393 | 0.9317 |

## GroupCV Stress Test 비교

| variant | input_set | model | n | spearman | pearson | rmse | mae | r2 | ndcg_at_20 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| baseline_mcf7_lincs | numeric_strong_context_smiles | ExtraTrees | 3387.0000 | 0.5112 | 0.4903 | 2.6134 | 1.9301 | 0.1822 | 0.7222 |
| no_lincs | numeric_strong_context_smiles_no_lincs | ExtraTrees | 3387.0000 | 0.4180 | 0.4391 | 2.7657 | 2.0235 | 0.0841 | 0.8615 |
| pancancer_lincs | numeric_strong_context_smiles_pan_lincs | ExtraTrees | 3387.0000 | 0.5036 | 0.5127 | 2.5695 | 1.9267 | 0.2094 | 0.8640 |

## Pan-Cancer LINCS Ensemble

- Primary input set: `numeric_strong_context_smiles_pan_lincs`
- Best single model: `LightGBM`
- Ensemble Spearman: `0.8787`
- Ensemble RMSE: `1.1962`
- Ensemble R2: `0.8287`

## 산출물

- `data/source_staging/lincs/lincs_drug_signature_pancancer_20260421.parquet`
- `data/X_numeric_pan_lincs.npy`
- `data/X_numeric_smiles_pan_lincs.npy`
- `data/X_numeric_strong_context_smiles_pan_lincs.npy`
- `results/pancancer_lincs/random3/`
- `results/pancancer_lincs/groupcv_stress_test/`
- `results/pancancer_lincs/ensemble/`
- `reports/qc_pancancer_lincs_signature_20260421.json`
- `reports/qc_pancancer_lincs_training_20260421.json`
