# Thyroid Cancer Pipeline Actual Run Summary 2026-04-20

이 문서는 `say2-4team/thyroid_raw/`에 모은 실제 원천 데이터를 사용해 갑상선암 약물 재창출 파이프라인을 end-to-end로 실행한 결과를 요약한다. 현재 요약은 `v2_smiles_required` 기준이며, primary 추천 pool에서는 SMILES가 없는 약물을 제거했다.

## 결론

현재 데이터로 thyroid cancer screened-response 기반 추천 파이프라인은 실행 가능한 상태다.

- 원천 screened response는 GDSC2 기반 `4,037`개 drug-cell line response pair, `16`개 THCA cell line, `295`개 screened drug로 시작했다.
- Primary 모델 입력은 SMILES missing drug `52`개를 제거한 뒤 `3,387`개 row, `16`개 cell line, `243`개 SMILES-valid drug로 구성했다.
- 가장 좋은 random sample 3-fold OOF는 `Numeric + SMILES / LightGBM`이며, Spearman `0.8833`, RMSE `1.1251`, R2 `0.8484`다.
- 최종 추천 입력셋인 `Numeric + Strong Context + SMILES / LightGBM`도 Spearman `0.8831`, RMSE `1.1226`, R2 `0.8491`로 거의 동일한 수준이다.
- 외부검증, ADMET, KG/clinical validation까지 통과한 최종 우선 후보는 `Vinorelbine`, `Staurosporine`, `Vinblastine`, `Romidepsin`이다.
- 단, random sample 3-fold 성능은 cell-line/drug overlap 영향을 받을 수 있으므로, 실제 일반화 성능은 GroupCV stress test를 함께 봐야 한다.

## Source Coverage

원천 데이터는 `s3://say2-4team/thyroid_raw/` 기준으로 정리했다.

| 항목 | 값 |
|---|---:|
| S3 object count | `132` |
| S3 total size | `1,973,354,972` bytes |
| THCA response rows before SMILES filter | `4,037` |
| THCA response rows after SMILES filter | `3,387` |
| THCA cell lines | `16` |
| Screened drugs before SMILES filter | `295` |
| SMILES-valid screened drugs | `243` |
| Removed SMILES-missing drugs | `52` |
| TCGA-THCA expression samples | `572` |
| ClinicalTrials.gov thyroid cancer drug studies | `570` |
| ADMET assays | `22` |

## Model-Ready Input QC

`scripts/02_build_model_ready_from_thyroid_raw.py`가 source staging에서 실제 모델 입력 파일을 생성한다.

| 입력/검사 | 값 |
|---|---:|
| Response table | `3,387` rows |
| Sample feature table | `16 x 4,114` |
| Drug feature table | `243 x 1,560` |
| Drug annotations | `243 x 11` |
| Numeric input dimension | `5,332` |
| SMILES SVD dimension | `64` |
| Strong context dimension | `32` |
| Numeric + SMILES dimension | `5,396` |
| Numeric + Strong Context + SMILES dimension | `5,428` |
| SMILES present / parse OK | `243 / 243` |
| Target gene present after SMILES filter | `212 / 243` |
| LINCS direct signature drugs | `101 / 243` |
| CRISPR feature cell lines | `9 / 16` |
| TCGA-THCA target expression genes written | `296` |
| TCGA-THCA clinical samples | `572` |

SMILES filter 결과:

| 기준 | Rows | Cell lines | Drugs |
|---|---:|---:|---:|
| Filter 전 screened response | `4,037` | `16` | `295` |
| Filter 후 primary model input | `3,387` | `16` | `243` |
| 제거된 SMILES missing set | `650` rows | `0` cell-line loss | `52` drugs |

## Random Sample 3-Fold OOF 성능

주 성능 비교는 random sample 3-fold OOF 기준이다.

### Numeric

| Model | Spearman | RMSE | R2 | NDCG@20 |
|---|---:|---:|---:|---:|
| LightGBM | `0.8779` | `1.1630` | `0.8380` | `0.9595` |
| XGBoost | `0.8598` | `1.2821` | `0.8032` | `0.9622` |
| RandomForest | `0.8583` | `1.2592` | `0.8101` | `0.9534` |
| LightGBM_DART | `0.8529` | `1.5890` | `0.6977` | `0.9532` |
| ExtraTrees | `0.8353` | `1.3649` | `0.7769` | `0.9550` |
| FlatMLP | `0.7811` | `1.6166` | `0.6871` | `0.9103` |

### Numeric + SMILES

| Model | Spearman | RMSE | R2 | NDCG@20 |
|---|---:|---:|---:|---:|
| LightGBM | `0.8833` | `1.1251` | `0.8484` | `0.9605` |
| XGBoost | `0.8756` | `1.2139` | `0.8236` | `0.9644` |
| LightGBM_DART | `0.8653` | `1.5302` | `0.7196` | `0.9524` |
| RandomForest | `0.8648` | `1.2354` | `0.8172` | `0.9460` |
| ExtraTrees | `0.8461` | `1.3186` | `0.7918` | `0.9564` |
| FlatMLP | `0.7900` | `1.5747` | `0.7031` | `0.9299` |

### Numeric + Strong Context + SMILES

| Model | Spearman | RMSE | R2 | NDCG@20 |
|---|---:|---:|---:|---:|
| LightGBM | `0.8831` | `1.1226` | `0.8491` | `0.9595` |
| XGBoost | `0.8761` | `1.2051` | `0.8261` | `0.9649` |
| LightGBM_DART | `0.8666` | `1.5218` | `0.7227` | `0.9551` |
| RandomForest | `0.8658` | `1.2342` | `0.8176` | `0.9450` |
| ExtraTrees | `0.8445` | `1.3239` | `0.7901` | `0.9555` |
| FlatMLP | `0.8034` | `1.5280` | `0.7204` | `0.9111` |

## v1 대비 v2 변화

SMILES missing drug를 유지하던 v1과 제거한 v2를 비교하면 random sample 3-fold에서는 성능이 크게 개선됐다.

| 기준 | Rows | Drugs | Primary input best model | Spearman | RMSE | R2 |
|---|---:|---:|---|---:|---:|---:|
| v1 all screened drugs | `4,037` | `295` | Numeric + Strong Context + SMILES / LightGBM | `0.8370` | `1.2689` | `0.7947` |
| v2 SMILES-required | `3,387` | `243` | Numeric + Strong Context + SMILES / LightGBM | `0.8831` | `1.1226` | `0.8491` |

해석:

- SMILES 0-vector 또는 invalid structure noise가 제거되면서 random3 기준 rank/회귀 성능이 좋아졌다.
- Cell line 수는 `16`개 그대로라 sample 축 손실은 없다.
- 다만 drug 수가 `295`개에서 `243`개로 줄기 때문에, 제거된 52개는 별도 recovery 후보로만 관리한다.

## GroupCV Stress Test

GroupCV는 random sample보다 훨씬 어려운 검증이며, unseen group generalization을 보기 위한 보수적 기준이다.

| Input set | Model | OOF Spearman | RMSE | R2 | NDCG@20 |
|---|---|---:|---:|---:|---:|
| Numeric + Strong Context + SMILES | ExtraTrees | `0.4515` | `2.7152` | `0.1172` | `0.6264` |

해석:

- random sample 3-fold 성능은 모델 비교와 후보 ranking 용도로 사용 가능하다.
- GroupCV Spearman은 v1 `0.4481`에서 v2 `0.4515`로 거의 유지됐지만, RMSE/R2는 악화됐다.
- 즉, SMILES missing 제거는 random split의 noisy input 문제를 줄였지만, 갑상선암의 핵심 병목은 여전히 `16`개 cell line이라는 sample-axis 제한이다.
- 따라서 최종 후보는 모델 점수만이 아니라 TCGA-THCA 외부검증, ADMET, KG/clinical evidence를 함께 통과한 후보로 좁히는 것이 안전하다.

## Ensemble Top Candidates

앙상블은 `Numeric + Strong Context + SMILES` 입력셋 모델 OOF 성능을 기준으로 weighted ensemble을 만든 뒤, 같은 약물명이 여러 ID로 반복되는 경우 normalized drug name 기준으로 중복 제거했다.

| Rank | Drug | Score | Classification | Main targets |
|---:|---|---:|---|---|
| 1 | Bortezomib | `4.1600` | indication_expansion | PROTEASOME; PRSS1; PSMB1; PSMB5 |
| 2 | Dactinomycin | `4.1477` | screened_candidate | TOP1; TOP2A |
| 3 | Docetaxel | `3.3544` | indication_expansion | BCL2; NR1I2; TUBB1 |
| 4 | Romidepsin | `3.2340` | indication_expansion | HDAC1; HDAC2; HDAC3; HDAC4; HDAC6; HDAC8 |
| 5 | Sepantronium bromide | `3.0155` | screened_candidate | BIRC5 |
| 6 | Staurosporine | `2.7263` | screened_candidate | CDK2; GSK3B; PIK3CG; PRKCA; SYK and others |
| 7 | SN-38 | `2.6638` | screened_candidate | TOP1 |
| 8 | Vinblastine | `2.5744` | screened_candidate | JUN; TUBA1A; TUBB family |
| 9 | Camptothecin | `2.1505` | indication_expansion | TOP1 |
| 10 | Vinorelbine | `2.0955` | indication_expansion | BCL2; MAP4; MAPK1; TUBB |

## External Validation, ADMET, KG 결과

외부검증은 TCGA-THCA expression/clinical/survival을 사용했다.

| QC | 값 |
|---|---:|
| External target gene match rate mean | `0.9109` |
| Target expressed count in top15 | `10` |
| Clinical data available | `true` |
| Top-k precision p@5 / p@10 / p@20 | `0.00 / 0.00 / 0.05` |
| ADMET candidate count | `15` |
| ADMET assay count | `22` |
| ADMET Approved / Candidate / Caution | `1 / 2 / 12` |
| KG Tier 1 / Tier 2 / Tier 3 / Excluded | `1 / 2 / 1 / 11` |

최종 comprehensive validation 결과:

| Tier | Drug | Ensemble score | ADMET category | Final category |
|---|---|---:|---|---|
| Tier 1 | Vinorelbine | `2.0955` | Candidate | Thyroid indication expansion candidate |
| Tier 2 | Staurosporine | `2.7263` | Caution | True repurposing or exploratory candidate |
| Tier 2 | Vinblastine | `2.5744` | Approved | True repurposing or exploratory candidate |
| Tier 3 | Romidepsin | `3.2340` | Candidate | Thyroid indication expansion candidate |

## Missingness 최소화 계획

현재 정책은 hard filter와 soft indicator를 분리한다.

| 결측 축 | 정책 | 이유 |
|---|---|---|
| SMILES missing drug | primary pool에서 제거 | 구조 표현, ADMET, 최종 검증 모두에서 필수성이 큼 |
| LINCS missing drug | 제거하지 않고 availability flag 유지 | hard filter 시 `243`개 drug가 `101`개로 줄어 후보 pool이 과도하게 작아짐 |
| CRISPR missing cell line | 제거하지 않고 fallback feature 보강 | hard filter 시 `16`개 cell line이 `9`개로 줄어 sample bottleneck이 심해짐 |

Local source audit 결과, 제거된 52개 SMILES-missing drug는 현재 local DrugBank/ChEMBL/LINCS name match만으로는 복구되지 않았다. 대신 LINCS는 local metadata에서 `16`개 bridge 후보가 있고, CRISPR missing 7개 cell line은 DepMap model bridge가 있으므로 expression/CNV/mutation fallback feature 보강이 다음 우선순위다.

## 재현 명령

```bash
python3 scripts/01_acquire_datasets.py --config config/thyroid_pipeline_config.json --manifest config/data_manifest.thyroid_raw_20260420.json
make build-model-ready
make audit-missingness
python3 scripts/run_all.py --config config/thyroid_pipeline_config.json
```

빠른 재계산이 필요한 경우, 모델 재학습 없이 후단 검증만 다시 돌릴 수 있다.

```bash
python3 scripts/run_all.py --config config/thyroid_pipeline_config.json --steps ensemble external admet knowledge report
```

## 산출물 위치

| 산출물 | 경로 |
|---|---|
| Model-ready source QC | `reports/qc_source_to_model_ready_20260420.json` |
| Missingness audit summary | `reports/missingness/missingness_source_audit_summary.json` |
| Model input QC | `reports/qc_step5_model_inputs.json` |
| Random3 metrics | `results/*_metrics_summary.csv` |
| GroupCV stress test | `results/groupcv_stress_test/numeric_strong_context_smiles_ExtraTrees_groupcv.json` |
| Ensemble top30 | `results/ensemble/thyroid_ensemble_top30_drugs.csv` |
| External validation | `external_validation/top15_validated.csv` |
| ADMET candidates | `admet/final_drug_candidates.csv` |
| Final candidates | `phase5_final_results/final_comprehensive_candidates.csv` |
| Final report | `phase5_final_results/FINAL_REPORT.html` |
