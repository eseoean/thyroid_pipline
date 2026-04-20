# Thyroid Cancer Hybrid Drug Repurposing Pipeline

이 레포는 갑상선암(thyroid cancer) 약물 재창출 프로젝트를 위한 end-to-end 파이프라인 구현입니다.

기준은 `docs/thyroid_cancer_hybrid_pipeline_guidance_20260420.md`에 정리된 hybrid 규칙입니다.

- Numeric base 입력셋은 Choi protocol v1 방향을 따릅니다.
- SMILES 표현과 strong context는 say2 BRCA 파이프라인 방식을 따릅니다.
- 모델 학습은 random sample 3-fold OOF 기준으로 수행합니다.
- 모델 출력 이후 외부검증, ADMET, knowledge validation, 최종 Tier 분류는 Choi protocol 흐름을 따릅니다.
- 최종 추천은 screened drug-response, 즉 실제 반응값이 있는 약물-샘플 pair에 한정합니다.

## 빠른 실행

원천 데이터 준비 상태를 먼저 확인합니다.

```bash
make inventory
```

데모 데이터로 전체 파이프라인을 확인하려면 다음을 실행합니다.

```bash
make run-demo
```

실제 데이터를 `config/thyroid_pipeline_config.json`에 지정된 경로에 넣은 뒤 실행하려면 다음을 실행합니다.

```bash
make run
```

S3 또는 다른 로컬 디렉터리에서 데이터를 내려받아야 하면 `docs/DATASET_ACQUISITION.md`를 참고해 `config/data_manifest.json`을 채운 뒤 `scripts/01_acquire_datasets.py`를 실행합니다. 이 스크립트는 로컬 다운로드/복사만 수행하며 S3 업로드는 하지 않습니다.

## 입력 데이터

기본 설정 파일은 `config/thyroid_pipeline_config.json`입니다.

필수 입력은 다음 형태를 권장합니다.

```text
data/raw/thyroid_response_pairs.csv
data/raw/sample_features.csv
data/raw/drug_features.csv
data/raw/drug_annotations.csv
data/external/thyroid_expression.csv
data/external/thyroid_clinical.csv
data/admet/*.csv
```

`thyroid_response_pairs.csv` 필수 컬럼:

- `sample_id`
- `canonical_drug_id`
- `LN_IC50` 또는 `IC50`
- `disease_label`
- `cell_line_name`
- `drug_name`

`drug_features.csv` 권장 컬럼:

- `canonical_drug_id`
- `drug_name`
- `canonical_smiles`
- `drug_morgan_*`
- `drug__lincs_score`
- `drug_desc_hba`, `drug_desc_hbd`, `drug_desc_heavy_atoms`, `drug_desc_ring_count`, `drug_desc_rot_bonds`

위 5개 물성 컬럼은 Choi 기준에 맞춰 numeric base에서는 제외됩니다.

## 산출물

주요 산출물은 다음 위치에 생성됩니다.

```text
data/X_numeric.npy
data/X_smiles_svd64.npy
data/X_strong_context.npy
data/X_numeric_smiles.npy
data/X_numeric_strong_context_smiles.npy
data/y_train.npy
data/row_metadata.csv

results/random3/
results/oof/
results/ensemble/thyroid_ensemble_top30_drugs.csv
results/ensemble/thyroid_ensemble_diversity.csv

external_validation/top15_validated.csv
admet/final_drug_candidates.csv
knowledge_validation/validation_summary.csv
phase5_final_results/final_comprehensive_candidates.csv
phase5_final_results/FINAL_REPORT.md
phase5_final_results/FINAL_REPORT.html
```

## 단계별 흐름

1. 질환 alias와 subtype을 정규화하여 thyroid cancer screened response subset을 만듭니다.
2. Choi 기준 numeric base를 생성하고 5개 물성 컬럼을 numeric에서 제외합니다.
3. SMILES를 char n-gram TF-IDF와 SVD 64차원으로 변환합니다.
4. `TCGA_DESC`, pathway, classification, bridge strength, stage3 status 기반 strong context를 생성합니다.
5. Numeric, Numeric+SMILES, Numeric+StrongContext+SMILES 세 입력셋을 구성합니다.
6. random sample 3-fold OOF로 LightGBM, LightGBM_DART, XGBoost, ExtraTrees, RandomForest, FlatMLP를 학습합니다.
7. OOF Spearman 기반 weighted ensemble과 diversity를 계산합니다.
8. THCA 외부 cohort expression/clinical 파일로 target expression과 survival signal을 확인합니다.
9. ADMET assay 파일을 Morgan fingerprint nearest-neighbor 방식으로 평가합니다.
10. Knowledge/clinical validation으로 Tier 1/2/3/Excluded를 분류하고 최종 리포트를 생성합니다.

## 실제 데이터 연결 체크리스트

- `config/thyroid_pipeline_config.json`의 입력 경로가 실제 파일을 가리키는지 확인합니다.
- response table은 갑상선암 cell line과 screened drug pair만 포함하거나, `disease_label`로 alias 필터링 가능해야 합니다.
- 외부검증 expression은 첫 컬럼이 gene symbol이고 나머지 컬럼이 환자 sample이어야 합니다.
- ADMET assay 파일은 `smiles`, `label` 컬럼을 권장합니다.
- 실행 후 `reports/qc_step*.json`에서 row 수, feature 수, SMILES coverage, context unknown 비율, fold QC, model gap을 확인합니다.
