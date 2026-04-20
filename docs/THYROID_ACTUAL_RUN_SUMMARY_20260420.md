# Thyroid Cancer Pipeline Actual Run Summary 2026-04-20

이 문서는 `say2-4team/thyroid_raw/`에 모은 실제 원천 데이터를 사용해 갑상선암 약물 재창출 파이프라인을 end-to-end로 실행한 결과를 요약한다.

## 결론

현재 데이터로 thyroid cancer screened-response 기반 추천 파이프라인은 실행 가능한 상태다.

- 추천의 주축은 GDSC2에서 실제 반응값이 존재하는 `295`개 screened drug이다.
- 최종 모델 입력은 `4,037`개 drug-cell line response pair와 `16`개 THCA cell line으로 구성됐다.
- 가장 좋은 random sample 3-fold OOF 모델은 `Numeric + Strong Context + SMILES / LightGBM`이다.
- 외부검증, ADMET, KG/clinical validation까지 통과한 최종 우선 후보는 `Vinorelbine`, `Staurosporine`, `Vinblastine`, `Romidepsin`이다.
- 단, random sample 3-fold 성능은 cell-line/drug overlap 영향을 받을 수 있으므로, 실제 일반화 성능은 GroupCV stress test를 함께 봐야 한다.

## Source Coverage

원천 데이터는 `s3://say2-4team/thyroid_raw/` 기준으로 정리했다.

| 항목 | 값 |
|---|---:|
| S3 object count | `132` |
| S3 total size | `1,973,354,972` bytes |
| THCA response rows | `4,037` |
| THCA cell lines | `16` |
| Screened drugs | `295` |
| TCGA-THCA expression samples | `572` |
| ClinicalTrials.gov thyroid cancer drug studies | `570` |
| ADMET assays | `22` |

## Model-Ready Input QC

`scripts/02_build_model_ready_from_thyroid_raw.py`가 source staging에서 실제 모델 입력 파일을 생성한다.

| 입력/검사 | 값 |
|---|---:|
| Response table | `4,037` rows |
| Sample feature table | `16 x 4,114` |
| Drug feature table | `295 x 1,560` |
| Drug annotations | `295 x 11` |
| Numeric input dimension | `5,334` |
| SMILES SVD dimension | `64` |
| Strong context dimension | `32` |
| Numeric + SMILES dimension | `5,398` |
| Numeric + Strong Context + SMILES dimension | `5,430` |
| SMILES present / parse OK | `243 / 243` |
| Target gene present | `238` |
| LINCS matched drugs | `101` |
| TCGA-THCA target expression genes written | `312` |
| TCGA-THCA clinical rows | `579` |

## Random Sample 3-Fold OOF 성능

주 성능 비교는 random sample 3-fold OOF 기준이다.

### Numeric

| Model | Spearman | RMSE | R2 | NDCG@20 |
|---|---:|---:|---:|---:|
| LightGBM | `0.8166` | `1.3474` | `0.7685` | `0.9647` |
| XGBoost | `0.7970` | `1.4460` | `0.7334` | `0.9540` |
| LightGBM_DART | `0.7896` | `1.7093` | `0.6274` | `0.9510` |
| RandomForest | `0.7861` | `1.4484` | `0.7325` | `0.9526` |
| ExtraTrees | `0.7643` | `1.5471` | `0.6948` | `0.9529` |
| FlatMLP | `0.7010` | `1.7924` | `0.5903` | `0.9265` |

### Numeric + SMILES

| Model | Spearman | RMSE | R2 | NDCG@20 |
|---|---:|---:|---:|---:|
| LightGBM | `0.8181` | `1.3343` | `0.7730` | `0.9681` |
| XGBoost | `0.8119` | `1.3931` | `0.7525` | `0.9618` |
| LightGBM_DART | `0.7988` | `1.6696` | `0.6445` | `0.9555` |
| RandomForest | `0.7914` | `1.4376` | `0.7365` | `0.9592` |
| ExtraTrees | `0.7717` | `1.5184` | `0.7060` | `0.9493` |
| FlatMLP | `0.7520` | `1.6211` | `0.6649` | `0.9436` |

### Numeric + Strong Context + SMILES

| Model | Spearman | RMSE | R2 | NDCG@20 |
|---|---:|---:|---:|---:|
| LightGBM | `0.8370` | `1.2689` | `0.7947` | `0.9697` |
| XGBoost | `0.8228` | `1.3549` | `0.7659` | `0.9616` |
| LightGBM_DART | `0.8096` | `1.6453` | `0.6548` | `0.9602` |
| RandomForest | `0.8046` | `1.3908` | `0.7534` | `0.9601` |
| ExtraTrees | `0.7820` | `1.4855` | `0.7186` | `0.9501` |
| FlatMLP | `0.7598` | `1.5791` | `0.6820` | `0.9399` |

## GroupCV Stress Test

GroupCV는 random sample보다 훨씬 어려운 검증이며, unseen group generalization을 보기 위한 보수적 기준이다.

| Input set | Model | OOF Spearman | RMSE | R2 | NDCG@20 |
|---|---|---:|---:|---:|---:|
| Numeric + Strong Context + SMILES | ExtraTrees | `0.4481` | `2.5048` | `0.2000` | `0.8018` |

해석:

- random sample 3-fold 성능은 모델 비교와 후보 ranking 용도로 사용 가능하다.
- GroupCV 성능은 실제 외삽 성능이 더 낮을 수 있음을 보여준다.
- 따라서 최종 후보는 모델 점수만이 아니라 TCGA-THCA 외부검증, ADMET, KG/clinical evidence를 함께 통과한 후보로 좁히는 것이 안전하다.

## Ensemble Top Candidates

앙상블은 `Numeric + Strong Context + SMILES` 입력셋 모델 OOF 성능을 기준으로 weighted ensemble을 만든 뒤, 같은 약물명이 여러 ID로 반복되는 경우 normalized drug name 기준으로 중복 제거했다.

| Rank | Drug | Score | Classification | Main targets |
|---:|---|---:|---|---|
| 1 | Bortezomib | `4.1811` | indication_expansion | PROTEASOME; PRSS1; PSMB1; PSMB5 |
| 2 | Dactinomycin | `4.1544` | screened_candidate | TOP1; TOP2A |
| 3 | Romidepsin | `3.9073` | indication_expansion | HDAC1; HDAC2; HDAC3; HDAC4; HDAC6; HDAC8 |
| 4 | Docetaxel | `3.5656` | indication_expansion | BCL2; NR1I2; TUBB1 |
| 5 | Sepantronium bromide | `3.3692` | screened_candidate | BIRC5 |
| 6 | Staurosporine | `3.1810` | screened_candidate | CDK2; GSK3B; PIK3CG; PRKCA; SYK and others |
| 7 | SN-38 | `2.5704` | screened_candidate | TOP1 |
| 8 | Vinblastine | `2.2527` | screened_candidate | JUN; TUBA1A; TUBB family |
| 9 | Camptothecin | `2.2271` | indication_expansion | TOP1 |
| 10 | Paclitaxel | `1.8365` | indication_expansion | BCL2; NR1I2; TUBB1 |

## External Validation, ADMET, KG 결과

외부검증은 TCGA-THCA expression/clinical/survival을 사용했다.

| QC | 값 |
|---|---:|
| External target gene match rate mean | `0.9052` |
| Target expressed count in top15 | `10` |
| Clinical data available | `true` |
| ADMET candidate count | `15` |
| ADMET assay count | `22` |
| ADMET Approved / Candidate / Caution | `1 / 2 / 12` |
| KG Tier 1 / Tier 2 / Tier 3 / Excluded | `1 / 2 / 1 / 11` |

최종 comprehensive validation 결과:

| Tier | Drug | Ensemble score | ADMET category | Final category |
|---|---|---:|---|---|
| Tier 1 | Vinorelbine | `1.6937` | Candidate | Thyroid indication expansion candidate |
| Tier 2 | Staurosporine | `3.1810` | Caution | True repurposing or exploratory candidate |
| Tier 2 | Vinblastine | `2.2527` | Approved | True repurposing or exploratory candidate |
| Tier 3 | Romidepsin | `3.9073` | Candidate | Thyroid indication expansion candidate |

## 재현 명령

```bash
python3 scripts/01_acquire_datasets.py --config config/thyroid_pipeline_config.json --manifest config/data_manifest.thyroid_raw_20260420.json
make build-model-ready
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
| Model input QC | `reports/qc_step5_model_inputs.json` |
| Random3 metrics | `results/*_metrics_summary.csv` |
| GroupCV stress test | `results/groupcv_stress_test/numeric_strong_context_smiles_ExtraTrees_groupcv.json` |
| Ensemble top30 | `results/ensemble/thyroid_ensemble_top30_drugs.csv` |
| External validation | `external_validation/top15_validated.csv` |
| ADMET candidates | `admet/final_drug_candidates.csv` |
| Final candidates | `phase5_final_results/final_comprehensive_candidates.csv` |
| Final report | `phase5_final_results/FINAL_REPORT.html` |
