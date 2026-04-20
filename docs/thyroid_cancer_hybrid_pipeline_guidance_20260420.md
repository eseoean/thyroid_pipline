# 갑상선암(Thyroid Cancer) 확장 프로젝트용 Hybrid Pipeline 지침

생성일: 2026-04-20  
대상 프로젝트: Thyroid Cancer Drug Repurposing Pipeline  
기준 프로토콜: Choi protocol v1 + say2 BRCA exact slim/strong context/SMILES pipeline

## 1. 목적

이 문서는 기존 BRCA 약물 재창출 파이프라인을 갑상선암(Thyroid cancer)으로 확장할 때 사용할 작업 지침이다.

이번 확장에서는 두 프로토콜을 섞어서 사용한다.

| 영역 | 따를 기준 | 이유 |
|---|---|---|
| Numeric base 입력셋 생성 | Choi protocol v1 | 6,366 row BRCA 기준에서 안정적으로 검증된 numeric baseline 구성 방식 |
| SMILES 표현 | say2 BRCA pipeline | SMILES SVD 64차원 방식과 기존 대시보드/모델 코드 호환성 |
| Context 표현 | say2 strong context pipeline | 단순 물성 context가 아니라 질환/경로/target-resolution/bridge-strength 기반 context 사용 |
| 모델 학습 | say2 BRCA pipeline | random sample 3-fold OOF, 다중 모델, weighted ensemble, diversity 분석 흐름 유지 |
| 외부검증/ADMET/최종 Tier | Choi protocol v1 | 모델 출력 이후 METABRIC-equivalent 검증, ADMET, knowledge validation, Tier 분류 구조가 명확함 |

중요 원칙:

- 최종 추천 파이프라인은 **반응값이 있는 screened drug-response 데이터만** 사용한다.
- 반응값이 없는 unscreened 후보는 메인 추천 결과에 섞지 않는다.
- Unscreened 확장은 별도 실험 브랜치 또는 appendix/future work로 분리한다.

## 2. 전체 구조

갑상선암 확장용 전체 흐름은 다음과 같다.

```text
Step 0. 질환 정의 및 데이터 inventory
Step 1. Thyroid cancer drug-response subset 생성
Step 2. Choi 기준 numeric base 입력셋 생성
Step 3. say2 기준 SMILES SVD 생성
Step 4. say2 기준 strong context 생성
Step 5. 모델 학습 및 OOF 생성
Step 6. 앙상블 및 diversity 분석
Step 7. 갑상선암 외부검증
Step 8. ADMET 안전성 평가
Step 9. Knowledge/clinical validation
Step 10. 최종 후보 Tier 분류 및 리포트
```

## 3. 질환 정의

갑상선암은 데이터 소스마다 이름이 다르게 들어갈 수 있으므로, 초기 단계에서 alias를 고정한다.

권장 disease key:

- `THCA`
- `thyroid`
- `thyroid cancer`
- `thyroid carcinoma`
- `papillary thyroid carcinoma`
- `follicular thyroid carcinoma`
- `medullary thyroid carcinoma`
- `anaplastic thyroid carcinoma`

가능하면 subtype도 별도 컬럼으로 남긴다.

| Subtype | 설명 |
|---|---|
| PTC | Papillary thyroid carcinoma |
| FTC | Follicular thyroid carcinoma |
| MTC | Medullary thyroid carcinoma |
| ATC | Anaplastic thyroid carcinoma |
| PDTC | Poorly differentiated thyroid carcinoma |
| Unknown | subtype 정보 없음 |

QC:

- Disease alias 매칭 결과를 `disease_subset_qc.csv`로 저장한다.
- 원본 disease label, normalized disease label, subtype, cell line ID를 모두 남긴다.
- subtype이 불명확한 경우 임의 보정하지 말고 `Unknown`으로 둔다.

## 4. Step 1: Drug-response subset 생성

목표는 thyroid cancer cell line × screened drug pair 단위의 supervised 학습 테이블을 만드는 것이다.

필수 컬럼:

- `sample_id`
- `cell_line_name`
- `canonical_drug_id`
- `drug_name`
- `IC50` 또는 `LN_IC50`
- `disease_label`
- `thyroid_subtype`

QC:

- 전체 row 수
- unique cell line 수
- unique drug 수
- cell line × drug 중복 row 수
- label missing 수
- label 분포 요약: min, q1, median, q3, max
- drug별 샘플 수 분포
- cell line별 약물 수 분포

권장 최소 기준:

- unique cell line 수가 너무 작으면 random split 성능이 매우 낙관적으로 보일 수 있다.
- unique drug 수가 적으면 Top 후보 다양성이 떨어진다.
- 특정 cell line 하나가 대부분 row를 차지하면 split leakage 가능성이 커진다.

산출물:

- `data/thyroid_response_pairs.parquet`
- `reports/qc_step1_response_subset.json`
- `reports/qc_step1_response_subset.html`

## 5. Step 2: Choi 기준 numeric base 입력셋 생성

이 단계는 Choi protocol v1의 numeric-only 입력셋 철학을 따른다.

Choi 기준:

```text
X_numeric.npy = numeric base only
약물 물성 5개는 numeric이 아니라 context 후보로 분리
```

BRCA Choi protocol에서는 numeric-only가 `5524`개였다. 갑상선암에서도 가능한 한 같은 feature family를 유지한다.

Choi에서 numeric에서 분리했던 5개 물성:

- `drug_desc_hba`
- `drug_desc_hbd`
- `drug_desc_heavy_atoms`
- `drug_desc_ring_count`
- `drug_desc_rot_bonds`

이번 hybrid에서는 이 5개를 numeric base에 넣지 않는다.  
단, 우리 strong context를 따르므로 이 5개를 그대로 context로 쓰지도 않는다. 필요하면 QC/보조 분석용 metadata로만 보관한다.

Numeric base 권장 구성:

- sample CRISPR / gene dependency feature
- drug Morgan fingerprint numeric feature
- LINCS numeric score
- target overlap / target expression / target pathway score
- 기타 수치형 feature 중 leakage가 아닌 것

제외해야 할 컬럼:

- ID 컬럼: `sample_id`, `canonical_drug_id`
- raw text 컬럼: SMILES raw, drug name raw
- label 컬럼: `IC50`, `LN_IC50`, response
- disease label 자체
- 이번 기준에서 분리할 5개 drug physicochemical context 후보

QC:

- feature count가 예상과 맞는지 확인한다.
- all-zero 컬럼 수를 기록한다.
- NaN/inf 값을 0 또는 median으로 처리한 경우 처리 전후 수를 기록한다.
- variance 0 컬럼을 제거했는지 여부를 명시한다.
- train/test split 전에 scaling/encoding을 fit하지 않는다.
- feature name list를 반드시 저장한다.

산출물:

- `data/X_numeric.npy`
- `data/y_train.npy`
- `data/numeric_feature_names.json`
- `reports/qc_step2_numeric_base.json`

주의:

BRCA say2 exact-slim은 numeric이 `5529`개였고, Choi는 `5524`개였다. 이 차이는 위 5개 물성 feature를 numeric에 포함하느냐의 차이다.  
갑상선암 hybrid에서는 **Choi 기준 numeric 5524 방향**을 따른다.

## 6. Step 3: say2 기준 SMILES SVD 생성

SMILES 처리는 say2 BRCA pipeline을 따른다.

권장 방식:

```text
canonical_smiles
→ character n-gram TF-IDF
→ TruncatedSVD 64 components
→ drug별 embedding 생성
→ row별 canonical_drug_id로 매핑
```

기본 설정:

- analyzer: `char`
- ngram range: `(2, 4)`
- SVD dimension: `64`
- random_state: 고정

QC:

- SMILES 보유율
- RDKit parse 성공률
- canonicalization 성공률
- missing SMILES drug list
- duplicate SMILES drug list
- SVD explained variance
- row 매핑 후 all-zero SMILES vector 수

산출물:

- `data/X_smiles_svd64.npy`
- `data/smiles_feature_names.json`
- `data/drug_smiles_qc.csv`
- `reports/qc_step3_smiles.json`

주의:

- SMILES가 없는 약물은 모델 학습에서 제거하지 말고, 우선 zero vector 또는 missing flag 정책을 고정한다.
- 단, ADMET 단계에서는 SMILES가 없으면 안전성 평가가 불가능하므로 별도 `NO_SMILES`로 표시한다.

## 7. Step 4: say2 기준 strong context 생성

Context는 Choi의 물성 context가 아니라 say2의 strong context를 따른다.

권장 strong context family:

- `TCGA_DESC`
- `PATHWAY_NAME_NORMALIZED`
- `classification`
- `drug_bridge_strength`
- `stage3_resolution_status`

BRCA say2에서는 32차원 one-hot context를 사용했다.

갑상선암에서는 다음처럼 변환한다.

| Context | 갑상선암 적용 |
|---|---|
| `TCGA_DESC` | `THCA` 또는 `UNCLASSIFIED` |
| `PATHWAY_NAME_NORMALIZED` | drug target/pathway annotation에서 생성 |
| `classification` | target gene resolution 상태 |
| `drug_bridge_strength` | DrugBank/ChEMBL/LINCS 등 multi-source 연결 강도 |
| `stage3_resolution_status` | target/pathway resolution 보정 상태 |

QC:

- context별 unique class 수
- one-hot dimension
- unknown/unclassified 비율
- pathway missing 비율
- target gene resolution 실패율
- bridge source별 coverage

산출물:

- `data/X_strong_context.npy`
- `data/strong_context_feature_names.json`
- `data/context_mapping_table.csv`
- `reports/qc_step4_strong_context.json`

중요:

Choi numeric 기준을 유지하면 최종 full input dimension은 아래처럼 된다.

```text
X_full = X_numeric(Choi-style 5524)
       + X_strong_context(say2-style 32)
       + X_smiles_svd64(64)
       = 5620 dimensions
```

따라서 기존 BRCA exact-slim `5625`와 차원이 다르다. 새 갑상선암 모델은 새 입력 차원 `5620` 기준으로 학습한다. 기존 BRCA 모델 artifact를 그대로 재사용하면 안 된다.

## 8. Step 5: 최종 입력셋 구성

세 가지 입력셋을 만든다.

| 입력셋 | 구성 | 목적 |
|---|---|---|
| A. Numeric | Choi-style numeric base | baseline |
| B. Numeric + SMILES | numeric + SMILES SVD64 | structure signal 추가 |
| C. Numeric + Strong Context + SMILES | numeric + say2 strong context + SMILES | 최종 후보 선정용 |

파일:

- `data/X_numeric.npy`
- `data/X_numeric_smiles.npy`
- `data/X_numeric_strong_context_smiles.npy`
- `data/y_train.npy`
- `data/row_metadata.parquet`

QC:

- 세 입력셋 row 수가 모두 같은지 확인한다.
- `y_train` row 수와 같은지 확인한다.
- `row_metadata`의 `sample_id`, `canonical_drug_id` 순서가 matrix row 순서와 같은지 확인한다.
- feature name 수와 matrix column 수가 같은지 확인한다.
- dtype은 가능하면 `float32`로 통일한다.

## 9. Step 6: 모델 학습

모델 학습은 say2 BRCA pipeline을 따른다.

기본 평가:

- random sample 3-fold OOF
- 각 fold별 train/validation metric 저장
- OOF prediction 저장
- model별 Spearman, Pearson, RMSE, MAE, R2, NDCG@20 계산

권장 모델:

ML:

- LightGBM
- LightGBM_DART
- XGBoost
- CatBoost
- ExtraTrees
- RandomForest

DL:

- FlatMLP
- ResidualMLP
- CrossAttention
- TabNet
- WideDeep
- FTTransformer 또는 TabTransformer

입력셋별 실행:

- Numeric only
- Numeric + SMILES
- Numeric + Strong Context + SMILES

QC:

- fold별 row 수 균형
- label 분포 fold별 비교
- drug/cell line coverage fold별 비교
- train metric과 validation metric gap
- OOF prediction NaN 여부
- model별 prediction variance
- model별 residual 분포

추가 권장 QC:

- GroupCV by `canonical_drug_id`를 보조 stress-test로 1회 실행한다.
- GroupCV 성능은 최종 대시보드의 primary metric으로 쓰지 않더라도, leakage risk 설명용으로 저장한다.

산출물:

- `results/{input_set}_{model}_random3.json`
- `results/{input_set}_oof/{model}.npy`
- `results/{input_set}_metrics_summary.csv`
- `reports/qc_step6_training.json`

## 10. Step 7: 앙상블 및 diversity

앙상블은 say2 pipeline 방식으로 수행한다.

권장 방식:

- OOF 성능 기반 weighted ensemble
- 후보 모델 간 prediction correlation 계산
- residual correlation 계산
- mean absolute prediction gap 계산
- 성능이 높아도 diversity가 너무 낮으면 보수적으로 해석

앙상블 후보:

- 전체 모델 weighted ensemble
- ML-only ensemble
- DL-only ensemble
- selected ensemble
- view-split ensemble은 후속 실험으로 분리

QC:

- ensemble OOF metric이 best single보다 개선되는지 확인
- diversity가 너무 낮으면 단순 평균의 의미가 약함
- fold별 ensemble 성능 편차 확인
- 특정 모델 하나가 weight 대부분을 차지하는지 확인

산출물:

- `results/thyroid_ensemble_results.json`
- `results/thyroid_ensemble_top30_drugs.csv`
- `results/thyroid_ensemble_diversity.csv`
- `reports/qc_step7_ensemble_diversity.json`

## 11. Step 8: 외부검증

모델 출력 이후 검증은 Choi protocol v1을 따른다.  
BRCA의 METABRIC 역할을 갑상선암에서는 THCA 외부 cohort가 대신한다.

외부검증 구성:

Method A. Target expression validation

- 후보 약물 target gene이 thyroid cancer cohort에서 발현되는지 확인
- 환자 중 target gene expression이 global median 이상인 비율 계산
- `target_expressed = pct_patients_expressing >= 0.30`

Method B. Survival / recurrence stratification

- target gene expression high vs low 그룹 분리
- OS, DFS, PFI, RFS 중 사용 가능한 endpoint 사용
- Mann-Whitney U 또는 log-rank test 사용
- p-value와 direction 기록

Method C. Known thyroid drug precision

Known thyroid cancer drug list를 positive control로 둔다.

예시 known/control 후보:

- Lenvatinib
- Sorafenib
- Cabozantinib
- Vandetanib
- Selpercatinib
- Pralsetinib
- Dabrafenib
- Trametinib
- Larotrectinib
- Entrectinib

QC:

- target gene alias 매핑률
- expression matrix gene match rate
- clinical endpoint missing rate
- survival 분석 가능 patient 수
- known thyroid drug hit rate
- topK precision: P@5, P@10, P@20

산출물:

- `external_validation/thyroid_target_expression.csv`
- `external_validation/thyroid_survival_validation.csv`
- `external_validation/thyroid_known_drug_precision.csv`
- `external_validation/top15_validated.csv`
- `reports/qc_step8_external_validation.json`

## 12. Step 9: ADMET 안전성 평가

ADMET은 Choi protocol v1의 Step 4/Phase 4 구조를 따른다.

기본 방식:

- TDC ADMET 22개 assay 사용
- Morgan fingerprint radius 2, 2048 bits
- Tanimoto similarity 기반 nearest neighbor matching

매칭 기준:

| Match type | 기준 |
|---|---|
| exact | similarity >= 0.99 |
| close analog | similarity >= 0.85 |
| analog | similarity >= 0.70 |
| no match | similarity < 0.70 |

QC:

- 후보별 SMILES 보유 여부
- assay별 valid molecule 수
- 후보별 exact/close/analog/no_match 수
- toxicity flag: Ames, DILI, hERG
- ADMET coverage
- NO_SMILES 후보 분리

산출물:

- `admet/admet_detailed_candidates.csv`
- `admet/admet_summary.json`
- `admet/final_drug_candidates.csv`
- `reports/qc_step9_admet.json`

## 13. Step 10: Knowledge / Clinical Validation

최종 knowledge validation은 Choi protocol v1의 Phase 5 구조를 따른다.

평가 축:

| 축 | 의미 |
|---|---|
| Drug-target evidence | 약물이 target을 실제로 조절하는 근거 |
| Target-thyroid relevance | target이 thyroid cancer와 관련 있는지 |
| Clinical evidence | thyroid cancer 임상시험/승인/연구 여부 |
| Mechanism rationale | 기전이 thyroid cancer biology와 맞는지 |
| Safety profile | ADMET 및 알려진 안전성 |

갑상선암 관련 주요 biology 예시:

- MAPK signaling
- BRAF V600E
- RET fusion / RET mutation
- NTRK fusion
- RAS mutation
- VEGFR/angiogenesis
- PI3K/AKT/mTOR
- DNA damage / cell-cycle

QC:

- evidence source별 hit count
- clinical trial evidence 여부
- target-disease relation evidence 여부
- mechanism rationale 수동 검토 여부
- 자동 점수와 수동 검토 결과 분리

산출물:

- `knowledge_validation/thyroid_knowledge_validation_results.json`
- `knowledge_validation/validation_summary.csv`
- `phase5_final_results/final_comprehensive_candidates.csv`
- `phase5_final_results/tier1_high_confidence.csv`
- `phase5_final_results/FINAL_REPORT.md`

## 14. 최종 후보 선정 규칙

최종 후보는 screened drug-response 기반 후보만 사용한다.

권장 Tier:

| Tier | 의미 |
|---|---|
| Tier 1 | 모델 성능, 외부검증, ADMET, knowledge validation이 모두 강한 후보 |
| Tier 2 | 일부 검증은 강하지만 추가 전임상 검증이 필요한 후보 |
| Tier 3 | 탐색 후보 또는 안전성/근거가 제한적인 후보 |
| Excluded | SMILES 없음, ADMET fail, target 근거 부족, label 불안정 |

최종 리포트에는 다음을 반드시 구분한다.

- Known thyroid cancer positive control
- Thyroid indication expansion candidate
- True repurposing candidate
- Exploratory candidate
- Excluded / caution candidate

## 15. 필수 QC 체크리스트

### Data QC

- [ ] thyroid disease alias mapping 확인
- [ ] subtype 분포 확인
- [ ] row 중복 제거 확인
- [ ] label missing 제거/처리 확인
- [ ] label outlier 기준 기록
- [ ] cell line 수와 drug 수 기록

### Feature QC

- [ ] numeric feature count 확인
- [ ] feature name과 matrix column 수 일치 확인
- [ ] NaN/inf 처리 전후 수 기록
- [ ] all-zero/constant feature 수 기록
- [ ] SMILES coverage 기록
- [ ] strong context unknown 비율 기록

### Split QC

- [ ] random3 fold row 수 균형
- [ ] fold별 label 분포 비교
- [ ] fold별 drug/cell line coverage 비교
- [ ] GroupCV stress-test 결과 별도 저장

### Model QC

- [ ] 모델별 OOF 저장
- [ ] fold별 metric 저장
- [ ] train-validation gap 확인
- [ ] ensemble weight 편중 확인
- [ ] diversity 지표 저장

### External Validation QC

- [ ] target gene match rate
- [ ] expression coverage
- [ ] survival endpoint missing rate
- [ ] known thyroid drug P@K
- [ ] ADMET coverage
- [ ] toxicity flag
- [ ] clinical/knowledge evidence 수동 검토

## 16. 권장 디렉터리 구조

```text
thyroid_project/
  data/
    thyroid_response_pairs.parquet
    X_numeric.npy
    X_smiles_svd64.npy
    X_strong_context.npy
    X_numeric_smiles.npy
    X_numeric_strong_context_smiles.npy
    y_train.npy
    row_metadata.parquet
    numeric_feature_names.json
    smiles_feature_names.json
    strong_context_feature_names.json

  results/
    random3/
    groupcv_stress_test/
    oof/
    ensemble/

  external_validation/
    thyroid_target_expression.csv
    thyroid_survival_validation.csv
    top15_validated.csv

  admet/
    admet_detailed_candidates.csv
    admet_summary.json
    final_drug_candidates.csv

  knowledge_validation/
    thyroid_knowledge_validation_results.json
    validation_summary.csv

  phase5_final_results/
    final_comprehensive_candidates.csv
    tier1_high_confidence.csv
    FINAL_REPORT.md

  reports/
    qc_step1_response_subset.json
    qc_step2_numeric_base.json
    qc_step3_smiles.json
    qc_step4_strong_context.json
    qc_step6_training.json
    qc_step7_ensemble_diversity.json
    qc_step8_external_validation.json
    qc_step9_admet.json
```

## 17. 구현 시 주의할 점

1. BRCA 모델 artifact를 갑상선암에 그대로 사용하지 않는다.

입력 차원, disease distribution, target relevance가 달라지므로 갑상선암용 모델은 새로 학습한다.

2. Choi numeric 기준과 say2 exact-slim 기준을 섞을 때 차원을 명확히 기록한다.

갑상선암 hybrid full input은 기본적으로 `5524 + 32 + 64 = 5620` 방향이다.  
만약 약물 물성 5개를 다시 numeric에 넣으면 기존 say2 BRCA처럼 `5625`가 되지만, 이는 이번 지침의 기본값이 아니다.

3. 최종 추천에는 반응값 없는 후보를 섞지 않는다.

Unscreened 후보는 hypothesis generation이며, supervised 성능 검증 대상이 아니다.

4. random3 성능과 GroupCV 성능을 혼동하지 않는다.

본 hybrid의 모델 학습/대시보드는 say2 random3 OOF 중심으로 가되, GroupCV는 leakage stress-test로 별도 보고한다.

5. 외부검증은 질환별 equivalent cohort를 사용한다.

BRCA의 METABRIC을 그대로 쓰면 안 된다. 갑상선암에서는 TCGA-THCA 또는 thyroid cancer expression/clinical cohort를 사용한다.

## 18. 발표/문서용 요약 문장

> 본 갑상선암 확장 프로젝트는 Choi protocol v1의 numeric base 생성 방식을 유지하면서, say2 BRCA pipeline의 SMILES SVD 및 strong context 표현을 결합한 hybrid 입력셋을 사용한다. 모델 학습은 say2 방식의 random sample 3-fold OOF 및 weighted ensemble로 수행하고, 모델 출력 이후의 외부검증, ADMET 안전성 평가, knowledge validation 및 최종 Tier 분류는 Choi protocol v1의 Phase 3~5 구조를 따른다. 최종 추천은 실제 약물 반응값이 존재하는 screened drug set에 한정하며, 반응값이 없는 unscreened 후보는 별도 참고 실험으로 분리한다.

