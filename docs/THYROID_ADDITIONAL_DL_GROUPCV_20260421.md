# Thyroid Additional DL GroupCV - 2026-04-21

## 목적

random sample 3-fold에서 좋게 나온 DL 모델들이 약물 ID 기준 GroupCV에서도 유지되는지 확인했다. GroupCV는 `canonical_drug_id`를 group으로 두기 때문에, validation fold의 약물은 train fold에서 빠진다. 따라서 random split보다 훨씬 보수적인 unseen-drug stress test이다.

## 실행 설정

- Input set: `numeric_strong_context_smiles_no_lincs`
- CV: GroupKFold 3-fold by `canonical_drug_id`
- Models: `ResidualMLP`, `WideDeep`, `CrossAttention`
- Device: PyTorch auto device

## DL GroupCV 성능

| input_set | model | cv | spearman | pearson | rmse | mae | r2 | ndcg_at_20 | train_oof_spearman_gap | elapsed_sec |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| numeric_strong_context_smiles_no_lincs | CrossAttention | groupcv_by_drug | 0.5190 | 0.5899 | 2.3489 | 1.7841 | 0.3393 | 0.8830 | 0.3983 | 11.4604 |
| numeric_strong_context_smiles_no_lincs | ResidualMLP | groupcv_by_drug | 0.4507 | 0.4995 | 2.5577 | 1.9059 | 0.2167 | 0.9058 | 0.3492 | 10.6261 |

## Fold별 성능

| model | fold | valid_drug_groups | best_epoch | valid_spearman | valid_rmse | valid_r2 |
| --- | --- | --- | --- | --- | --- | --- |
| CrossAttention | 1.0000 | 82.0000 | 3.0000 | 0.5670 | 2.1910 | 0.3942 |
| CrossAttention | 2.0000 | 80.0000 | 22.0000 | 0.4619 | 2.1672 | 0.2722 |
| CrossAttention | 3.0000 | 81.0000 | 13.0000 | 0.5136 | 2.6562 | 0.3262 |
| ResidualMLP | 1.0000 | 82.0000 | 10.0000 | 0.5154 | 2.3724 | 0.2897 |
| ResidualMLP | 2.0000 | 80.0000 | 7.0000 | 0.4239 | 2.2811 | 0.1937 |
| ResidualMLP | 3.0000 | 81.0000 | 12.0000 | 0.4180 | 2.9654 | 0.1602 |

## 기존 ML GroupCV 기준

- Baseline model: ``
- Baseline Spearman: ``
- Baseline RMSE: ``
- Baseline R2: ``

## 해석

GroupCV는 random sample 3-fold 성능보다 낮게 나오는 것이 정상이다. 여기서 성능이 유지되는 모델은 drug identity leakage에 덜 의존하고, unseen-drug 일반화에 상대적으로 강하다고 볼 수 있다.

## 산출물

- `results/additional_dl/numeric_strong_context_smiles_no_lincs/groupcv/additional_dl_groupcv_metrics_summary.csv`
- `results/additional_dl/numeric_strong_context_smiles_no_lincs/groupcv/additional_dl_groupcv_fold_metrics.csv`
- `results/additional_dl/numeric_strong_context_smiles_no_lincs/groupcv/oof/`
- `reports/qc_additional_dl_groupcv_20260421.json`
