# Thyroid Cancer Pipeline Actual Run Summary 2026-04-20

이 문서는 `say2-4team/thyroid_raw/`에 모은 실제 원천 데이터와 로컬 추가 확보 source를 사용해 갑상선암 약물 재창출 파이프라인을 end-to-end로 실행한 결과를 요약한다. 현재 요약은 `v3_missingness_enhanced` 기준이다.

## 결론

현재 thyroid cancer screened-response 기반 추천 파이프라인은 최종 후보 추천까지 실행 가능한 상태다.

- 원천 screened response는 GDSC2 기반 `4,037`개 drug-cell line response pair, `16`개 THCA cell line, `295`개 screened drug로 시작했다.
- Primary 모델 입력은 SMILES missing drug `52`개를 제거한 뒤 `3,387`개 row, `16`개 cell line, `243`개 SMILES-valid drug로 구성했다.
- 이번 v3에서는 LINCS MCF7 signature를 재매칭해 LINCS coverage를 `101 / 243`에서 `115 / 243`으로 올렸다.
- CRISPR가 없는 `7`개 cell line은 DepMap 24Q2 expression/CNV/mutation fallback feature로 `7 / 7` 모두 보강했다.
- 가장 좋은 random sample 3-fold OOF는 `Numeric + Strong Context + SMILES / LightGBM`이며, Spearman `0.8835`, RMSE `1.1264`, R2 `0.8481`이다.
- GroupCV stress test는 v2 `0.4515`에서 v3 `0.5112`로 개선되어, 결측 보강이 random split보다 unseen-drug stress에서 더 의미 있게 작동했다.
- 외부검증, ADMET, KG/clinical validation까지 통과한 최종 우선 후보는 `Vinorelbine`, `Staurosporine`, `Vinblastine`, `Romidepsin`이다.

주의할 점은 random sample 3-fold 성능이 여전히 cell-line/drug overlap 영향을 받을 수 있다는 것이다. 따라서 모델 성능 수치는 후보 ranking 용도로 보고, 실제 최종 추천은 TCGA-THCA 외부검증, ADMET, KG/clinical evidence를 통과한 Tier 후보를 기준으로 해석한다.

## Source Coverage

원천 데이터는 `s3://say2-4team/thyroid_raw/` 기준으로 정리했고, 큰 보강 source는 로컬에만 내려받아 사용했다.

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

추가 로컬 source:

| Source | 용도 |
|---|---|
| `data/source_staging/lincs/lincs_mcf7.parquet` | LINCS bridge 후보의 실제 MCF7 perturbation signature 재매칭 |
| `data/source_staging/depmap/OmicsExpressionProteinCodingGenesTPMLogp1_24Q2.csv` | CRISPR-missing cell line expression fallback |
| `data/source_staging/depmap/OmicsCNGene_24Q2.csv` | CRISPR-missing cell line CNV fallback |
| `data/source_staging/depmap/OmicsSomaticMutations_24Q2.csv` | CRISPR-missing cell line mutation fallback |

## Model-Ready Input QC

`scripts/02_build_model_ready_from_thyroid_raw.py`가 source staging에서 실제 모델 입력 파일을 생성한다.

| 입력/검사 | 값 |
|---|---:|
| Response table | `3,387` rows |
| Sample feature table | `16 x 5,360` |
| Drug feature table | `243 x 1,563` |
| Drug annotations | `243 x 11` |
| Numeric input dimension | `6,559` |
| SMILES SVD dimension | `64` |
| Strong context dimension | `32` |
| Numeric + SMILES dimension | `6,623` |
| Numeric + Strong Context + SMILES dimension | `6,655` |
| SMILES present / parse OK | `243 / 243` |
| Target gene present after SMILES filter | `212 / 243` |
| LINCS direct signature drugs | `101 / 243` |
| LINCS recovered signature drugs | `14 / 243` |
| LINCS total signature drugs | `115 / 243` |
| CRISPR feature cell lines | `9 / 16` |
| CRISPR-missing cell lines with omics fallback | `7 / 7` |
| Selected non-CRISPR omics fallback features | `1,234` |
| TCGA-THCA target expression genes written | `296` |
| TCGA-THCA clinical samples | `572` |

SMILES filter 결과:

| 기준 | Rows | Cell lines | Drugs |
|---|---:|---:|---:|
| Filter 전 screened response | `4,037` | `16` | `295` |
| Filter 후 primary model input | `3,387` | `16` | `243` |
| 제거된 SMILES missing set | `650` rows | `0` cell-line loss | `52` drugs |

## Missingness Enhancement

이번 v3에서 실제 반영한 보강은 두 축이다.

| 축 | v2 상태 | v3 조치 | v3 결과 |
|---|---:|---|---:|
| LINCS drug signature | `101 / 243` | local LINCS metadata bridge 16개를 MCF7 signature parquet와 재매칭 | `115 / 243` |
| CRISPR-missing cell line | `7`개 결측 | DepMap expression/CNV/mutation fallback feature 추가 | `7 / 7` 보강 |
| SMILES-missing drug | `52`개 제거 | local DrugBank/ChEMBL/LINCS name match 재검토 | 복구 가능 `0`개 |

LINCS bridge 후보 `16`개 중 `14`개는 실제 MCF7 signature row를 찾아 평균 signature로 복구했다. 복구되지 않은 후보는 `Dabrafenib`, `IWP-2`이며, 사유는 현재 로컬 MCF7 signature row 부재다.

## Random Sample 3-Fold OOF 성능

주 성능 비교는 random sample 3-fold OOF 기준이다.

### Numeric

| Model | Spearman | RMSE | R2 | NDCG@20 |
|---|---:|---:|---:|---:|
| LightGBM | `0.8776` | `1.1587` | `0.8392` | `0.9692` |
| RandomForest | `0.8680` | `1.2172` | `0.8226` | `0.9685` |
| XGBoost | `0.8617` | `1.2398` | `0.8160` | `0.9646` |
| ExtraTrees | `0.8587` | `1.2459` | `0.8141` | `0.9602` |
| LightGBM_DART | `0.8483` | `1.5652` | `0.7067` | `0.9636` |
| FlatMLP | `0.8172` | `1.4821` | `0.7370` | `0.9110` |

### Numeric + SMILES

| Model | Spearman | RMSE | R2 | NDCG@20 |
|---|---:|---:|---:|---:|
| LightGBM | `0.8833` | `1.1353` | `0.8457` | `0.9654` |
| XGBoost | `0.8765` | `1.1895` | `0.8306` | `0.9688` |
| RandomForest | `0.8765` | `1.1815` | `0.8328` | `0.9680` |
| LightGBM_DART | `0.8641` | `1.5133` | `0.7258` | `0.9609` |
| ExtraTrees | `0.8629` | `1.2335` | `0.8178` | `0.9599` |
| FlatMLP | `0.8253` | `1.4597` | `0.7448` | `0.9122` |

### Numeric + Strong Context + SMILES

| Model | Spearman | RMSE | R2 | NDCG@20 |
|---|---:|---:|---:|---:|
| LightGBM | `0.8835` | `1.1264` | `0.8481` | `0.9689` |
| XGBoost | `0.8767` | `1.1825` | `0.8326` | `0.9663` |
| RandomForest | `0.8763` | `1.1826` | `0.8325` | `0.9679` |
| LightGBM_DART | `0.8638` | `1.5120` | `0.7263` | `0.9599` |
| ExtraTrees | `0.8635` | `1.2295` | `0.8190` | `0.9651` |
| FlatMLP | `0.8287` | `1.4382` | `0.7523` | `0.9223` |

## v2 대비 v3 변화

v3 보강은 random3 최고 성능을 크게 올리기보다는, 결측 때문에 약했던 unseen-drug stress를 개선하는 효과가 컸다.

| 기준 | Rows | Drugs | Primary input best model | Spearman | RMSE | R2 |
|---|---:|---:|---|---:|---:|---:|
| v2 SMILES-required | `3,387` | `243` | Numeric + Strong Context + SMILES / LightGBM | `0.8831` | `1.1226` | `0.8491` |
| v3 missingness-enhanced | `3,387` | `243` | Numeric + Strong Context + SMILES / LightGBM | `0.8835` | `1.1264` | `0.8481` |

GroupCV stress:

| 기준 | Input set | Model | OOF Spearman | RMSE | R2 | NDCG@20 |
|---|---|---|---:|---:|---:|---:|
| v2 SMILES-required | Numeric + Strong Context + SMILES | ExtraTrees | `0.4515` | `2.7152` | `0.1172` | `0.6264` |
| v3 missingness-enhanced | Numeric + Strong Context + SMILES | ExtraTrees | `0.5112` | `2.6134` | `0.1822` | `0.7222` |

해석:

- random split에서는 이미 SMILES filter만으로 성능이 높았기 때문에 v3 추가 보강의 상승폭은 작다.
- GroupCV에서는 LINCS와 cell-line omics fallback이 unseen-drug stress에 도움을 줬다.
- 따라서 v3는 단순 점수 상승보다 결측 구조 완화와 보수적 검증 안정성 개선으로 해석하는 것이 맞다.

## Ensemble Top Candidates

앙상블은 `Numeric + Strong Context + SMILES` 입력셋 모델 OOF 성능을 기준으로 weighted ensemble을 만든 뒤, 같은 약물명이 여러 ID로 반복되는 경우 normalized drug name 기준으로 중복 제거했다.

| Rank | Drug | Score | Classification | Main targets |
|---:|---|---:|---|---|
| 1 | Bortezomib | `4.0652` | indication_expansion | PROTEASOME; PRSS1; PSMB1; PSMB5 |
| 2 | Dactinomycin | `3.9444` | screened_candidate | TOP1; TOP2A |
| 3 | Docetaxel | `3.5305` | indication_expansion | BCL2; NR1I2; TUBB1 |
| 4 | Sepantronium bromide | `2.7171` | screened_candidate | BIRC5 |
| 5 | Vinblastine | `2.6649` | screened_candidate | JUN; TUBA1A; TUBB family |
| 6 | Staurosporine | `2.6504` | screened_candidate | CDK2; GSK3B; PIK3CG; PRKCA; SYK and others |
| 7 | Romidepsin | `2.5290` | indication_expansion | HDAC1; HDAC2; HDAC3; HDAC4; HDAC6; HDAC8 |
| 8 | Vinorelbine | `2.4301` | indication_expansion | BCL2; MAP4; MAPK1; TUBB |
| 9 | SN-38 | `2.3780` | screened_candidate | TOP1 |
| 10 | Daporinad | `2.1312` | screened_candidate | NAMPT |

모델 raw score 상위는 Bortezomib, Dactinomycin, Docetaxel이지만, ADMET와 KG/clinical validation을 통과한 최종 Tier는 아래 후보로 좁혀진다.

## External Validation, ADMET, KG 결과

외부검증은 TCGA-THCA expression/clinical/survival을 사용했다.

| QC | 값 |
|---|---:|
| External target gene match rate mean | `0.8764` |
| Target expressed count in top15 | `10` |
| Clinical data available | `true` |
| Top-k precision p@5 / p@10 / p@20 | `0.00 / 0.00 / 0.05` |
| ADMET candidate count | `15` |
| ADMET assay count | `22` |
| ADMET Approved / Candidate / Caution | `1 / 2 / 12` |
| KG Tier 1 / Tier 2 / Tier 3 / Excluded | `1 / 2 / 1 / 11` |

최종 comprehensive validation 결과:

| Tier | Drug | Ensemble score | ADMET category | Knowledge score | Final category |
|---|---|---:|---|---:|---|
| Tier 1 | Vinorelbine | `2.4301` | Candidate | `4.5` | Thyroid indication expansion candidate |
| Tier 2 | Staurosporine | `2.6504` | Caution | `3.6` | True repurposing or exploratory candidate |
| Tier 2 | Vinblastine | `2.6649` | Approved | `3.2` | True repurposing or exploratory candidate |
| Tier 3 | Romidepsin | `2.5290` | Candidate | `2.8` | Thyroid indication expansion candidate |

Excluded 상위 raw-score 후보:

| Drug | Ensemble score | Exclusion reason signal |
|---|---:|---|
| Bortezomib | `4.0652` | ADMET Caution, KG score `2.4` |
| Dactinomycin | `3.9444` | ADMET Caution, KG score `2.4` |
| Docetaxel | `3.5305` | ADMET Caution, KG score `1.9` |
| Sepantronium bromide | `2.7171` | ADMET Caution, KG score `1.9` |

## Missingness 정책

현재 정책은 hard filter와 soft indicator를 분리한다.

| 결측 축 | 정책 | 이유 |
|---|---|---|
| SMILES missing drug | primary pool에서 제거 | 구조 표현, ADMET, 최종 검증 모두에서 필수성이 큼 |
| LINCS missing drug | 제거하지 않고 availability flag 유지 | hard filter 시 후보 pool이 과도하게 작아짐 |
| CRISPR missing cell line | 제거하지 않고 fallback feature 보강 | hard filter 시 `16`개 cell line이 `9`개로 줄어 sample bottleneck이 심해짐 |

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
| LINCS recovery review | `reports/missingness/lincs_signature_recovery_review.csv` |
| Model input QC | `reports/qc_step5_model_inputs.json` |
| Random3 metrics | `results/*_metrics_summary.csv` |
| GroupCV stress test | `results/groupcv_stress_test/numeric_strong_context_smiles_ExtraTrees_groupcv.json` |
| Ensemble top30 | `results/ensemble/thyroid_ensemble_top30_drugs.csv` |
| External validation | `external_validation/top15_validated.csv` |
| ADMET candidates | `admet/final_drug_candidates.csv` |
| Final candidates | `phase5_final_results/final_comprehensive_candidates.csv` |
| Final report | `phase5_final_results/FINAL_REPORT.html` |
