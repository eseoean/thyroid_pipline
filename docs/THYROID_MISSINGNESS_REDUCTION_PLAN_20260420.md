# Thyroid Missingness Reduction Plan 2026-04-20

이 문서는 thyroid cancer screened-response 파이프라인에서 확인된 결측 구조를 줄이기 위한 실행 계획과 2026-04-20 v3 적용 결과를 정리한다.

## 현재 결측 구조

SMILES hard filter 적용 전후:

| 기준 | Rows | Cell lines | Drugs |
|---|---:|---:|---:|
| 기존 전체 screened response | `4,037` | `16` | `295` |
| SMILES valid drug만 유지 | `3,387` | `16` | `243` |
| 제거된 SMILES missing drug | `650` rows | `0` cell-line loss | `52` drugs |

SMILES filter 이후에도 남는 주요 결측:

| 결측 축 | 현재 상태 | 해석 |
|---|---:|---|
| CRISPR sample coverage | `9 / 16` cell lines, fallback 적용 후 `16 / 16` | CRISPR 자체는 9개지만 non-CRISPR omics로 7개 보강 |
| LINCS drug coverage | direct `101 / 243`, recovered 포함 `115 / 243` | drug perturbation view가 개선됐지만 여전히 절반 미만 |
| Target gene coverage | `212 / 243` SMILES-valid drugs | strong context/KG evidence 약한 약물 존재 |

## 정책 결정

### 1. SMILES missing drug는 primary pool에서 제거

SMILES가 없는 약물은 primary screened-response 추천 후보에서 제외한다.

이유:

- SMILES 0-vector가 `Numeric+SMILES` 입력셋의 의미를 흐린다.
- 제거해도 cell line 수는 `16`개 그대로 유지된다.
- BRCA exact slim과 동일하게 SMILES-valid drug 중심의 추천 pool을 만들 수 있다.
- 제거 대상 52개는 내부 screening code, 숫자 compound, weak bridge compound가 많아 post-validation 근거 확보도 어렵다.

구현:

- `config/thyroid_pipeline_config.json`
  - `data_filters.require_valid_smiles_drugs = true`
- `scripts/02_build_model_ready_from_thyroid_raw.py`
  - source-to-model-ready 단계에서 SMILES parse OK 약물만 response/drug/annotation output에 유지
- `thyroid_pipeline/core.py`
  - downstream safety filter로 SMILES valid row만 한 번 더 확인

### 2. LINCS missing drug는 제거하지 않음

LINCS까지 hard filter로 적용하면 `243`개 drug가 `101`개로 줄어든다. 최종 후보권에도 LINCS missing 약물이 포함되므로, LINCS는 필수 조건이 아니라 보조 view로 처리해야 한다.

권장 처리:

- `has_lincs_signature` availability flag 추가
- LINCS feature 결측은 imputation + missing indicator로 처리
- LINCS가 있는 약물만 쓰는 별도 view model을 만들고 전체 모델과 ensemble
- LINCS missing drug를 구조 유사체 기반으로 보강하되, 실제 measured signature와 imputed signature를 구분

### 3. CRISPR missing cell line도 제거하지 않음

CRISPR 있는 cell line만 남기면 `16`개 cell line이 `9`개로 줄어든다. thyroid는 이미 sample 축이 얇기 때문에 CRISPR hard filter는 위험하다.

권장 처리:

- `sample_has_crispr` availability flag 유지
- CRISPR 없는 cell line은 expression/CNV/mutation/lineage feature로 보강
- CRISPR-rich model과 all-sample model을 분리한 view-split ensemble 구성

## 적용 결과 요약

v3에서 실제 반영한 보강은 다음과 같다.

| 보강 항목 | 입력 source | 결과 |
|---|---|---:|
| LINCS bridge 재매칭 | `lincs_pert_info_basic_20260406.parquet`, `lincs_mcf7.parquet` | `16`개 bridge 후보 중 `14`개 recovered |
| Cell line expression fallback | DepMap 24Q2 `OmicsExpressionProteinCodingGenesTPMLogp1` | `13 / 16` cell line matched |
| Cell line CNV fallback | DepMap 24Q2 `OmicsCNGene` | `16 / 16` cell line matched |
| Cell line mutation fallback | DepMap 24Q2 `OmicsSomaticMutations` | `16 / 16` cell line matched |
| CRISPR-missing cell line 보강 | expression/CNV/mutation fallback union | `7 / 7` covered |

모델 영향:

| 기준 | Primary random3 Spearman | GroupCV Spearman | 해석 |
|---|---:|---:|---|
| v2 SMILES-required | `0.8831` | `0.4515` | SMILES filter로 random split 성능 개선 |
| v3 missingness-enhanced | `0.8835` | `0.5112` | random split 유지, GroupCV stress 개선 |

즉 v3 보강은 최고 random3 점수를 크게 올리는 작업이라기보다, 결측 축을 줄여 unseen-drug stress 안정성을 높이는 작업이었다.

## 추가 확보 우선순위

## Local Source Audit 결과

`scripts/03_audit_missingness_sources.py`로 현재 local staging source에서 즉시 복구 가능한 범위를 확인했다.

| Audit 항목 | 결과 |
|---|---:|
| SMILES 제거 drug | `52` |
| 제거 drug 중 local DrugBank/ChEMBL/LINCS name match로 SMILES 복구 가능 | `0` |
| SMILES-valid drug | `243` |
| LINCS signature drug total | `115 / 243` |
| LINCS direct signature drug | `101 / 243` |
| LINCS recovered signature drug | `14 / 243` |
| LINCS local name/SMILES bridge 후보 remaining | `2` |
| Cell lines | `16` |
| CRISPR 있는 cell line | `9 / 16` |
| CRISPR는 없지만 DepMap model bridge가 있는 cell line | `7 / 7` |
| CRISPR missing이면서 omics fallback으로 보강된 cell line | `7 / 7` |

해석:

- 제거된 52개 약물은 현재 local source만으로는 SMILES 복구가 어렵다. PubChem/ChEMBL web/API 추가 조회가 필요하다.
- LINCS는 direct signature가 101개였고, local LINCS metadata bridge 후보 16개 중 실제 MCF7 signature가 있는 14개를 복구했다.
- CRISPR missing 7개 cell line은 DepMap expression/CNV/mutation fallback으로 모두 보강했다.
- 남은 LINCS bridge 후보 2개는 현재 로컬 MCF7 signature row가 없어 별도 source 또는 analog imputation이 필요하다.

Audit 산출물:

- `reports/missingness/smiles_removed_drug_recovery_review.csv`
- `reports/missingness/lincs_drug_mapping_review.csv`
- `reports/missingness/cellline_feature_coverage_review.csv`
- `reports/missingness/missingness_source_audit_summary.json`

### Priority A. Cell line feature 보강 - 적용 완료

목표:

- CRISPR가 없는 7개 cell line의 sample biology 결측을 줄인다.
- CRISPR 자체가 없더라도 expression, mutation, CNV, lineage/subtype feature를 확보한다.

대상 cell line:

| Missing CRISPR cell line |
|---|
| `CGTH-W-1` |
| `HTC-C3` |
| `K5` |
| `KMH-2` |
| `ML-1` |
| `RO82-W-1` |
| `TT` |

후보 source:

- DepMap 최신 public release
  - Model annotation
  - OmicsExpression
  - OmicsSomaticMutations
  - OmicsCNGene
  - CRISPRGeneDependency, 가능한 경우 최신 버전 재확인
- Cell Model Passports / Sanger GDSC cell model data
  - Sanger model ID, COSMIC ID, tissue, subtype, mutation, CNV bridge
- CCLE legacy expression / mutation / copy-number tables
  - DepMap에서 누락된 legacy cell line alias 보강

매칭 key 우선순위:

1. `SANGER_MODEL_ID`
2. `COSMIC_ID`
3. `DepMap ModelID`
4. normalized `CELL_LINE_NAME`
5. alias/synonym table

산출물:

- `reports/missingness/cellline_feature_coverage.csv`
- `reports/missingness/cellline_alias_mapping_review.csv`
- `data/raw/sample_features_enhanced.csv`
- `reports/qc_cellline_feature_enhancement.json`

현재 구현 산출물:

- `data/raw/sample_features.csv`
- `reports/qc_source_to_model_ready_20260420.json`
- `reports/missingness/cellline_feature_coverage_review.csv`

### Priority B. LINCS drug signature 보강 - 1차 적용 완료

목표:

- 현재 `101 / 243` drug만 LINCS에 연결되는 문제를 개선한다.
- 직접 measured LINCS match와 구조 유사체 기반 imputed match를 구분한다.

후보 source:

- LINCS/CLUE compound metadata
  - `pert_id`, `pert_iname`, `canonical_smiles`, `InChIKey`, aliases
- 현재 staging의 LINCS tables
  - `lincs_pert_info_basic_20260406.parquet`
  - `lincs_drug_signature_normalized.parquet`
- DrugBank / ChEMBL / PubChem bridge
  - drug name, synonym, ChEMBL ID, DrugBank ID, InChIKey, canonical SMILES

매칭 방식:

1. exact ID match
2. normalized name match
3. synonym match
4. InChIKey first block match
5. canonical SMILES exact match
6. Morgan similarity nearest neighbor

Nearest-neighbor 보강 기준:

| Similarity | 처리 |
|---:|---|
| `>= 0.99` | exact structural substitute |
| `>= 0.85` | close analog imputation |
| `0.70 - 0.85` | weak analog, 별도 flag 필요 |
| `< 0.70` | impute 금지 |

산출물:

- `reports/missingness/lincs_drug_mapping_review.csv`
- `reports/missingness/lincs_coverage_before_after.json`
- `data/raw/drug_features_enhanced.csv`
- `data/raw/lincs_signature_imputation_flags.csv`

현재 구현 산출물:

- `data/raw/drug_features.csv`
- `reports/missingness/lincs_signature_recovery_review.csv`
- `reports/missingness/lincs_drug_mapping_review.csv`
- `reports/qc_source_to_model_ready_20260420.json`

### Priority C. Target/KG context 보강

목표:

- target gene missing drug를 줄여 strong context와 KG validation 품질을 높인다.

후보 source:

- ChEMBL mechanism/action table
- DrugBank target table
- OpenTargets target-disease association
- PubChem synonym and compound identifier bridge
- UniProt target synonym mapping

처리 방식:

- GDSC putative target을 gene symbol로 정규화
- non-gene term과 pathway term은 target gene에서 제외
- ambiguous target은 review table로 분리
- final target source priority를 명시

산출물:

- `reports/missingness/target_mapping_review.csv`
- `reports/missingness/target_coverage_before_after.json`
- `data/raw/drug_annotations_enhanced.csv`

## 모델링 적용안

### v2 primary

적용 완료한 기준:

- SMILES valid drug만 사용
- LINCS/CRISPR missing은 제거하지 않음
- availability flags 유지

확인된 결과:

- SMILES 0-vector 제거
- BRCA exact slim과 더 유사한 후보 pool
- sample 16개 유지
- `Numeric + Strong Context + SMILES / LightGBM` random3 Spearman `0.8370`에서 `0.8831`로 개선
- 같은 기준 RMSE `1.2689`에서 `1.1226`으로 개선
- GroupCV Spearman은 `0.4481`에서 `0.4515`로 거의 유지됐지만, RMSE/R2는 악화되어 sample-axis 일반화 병목은 남아 있음

### v3 missingness-enhanced

적용 완료한 기준:

- SMILES valid drug만 primary pool 유지
- DepMap 24Q2 expression/CNV/mutation fallback feature 사용
- LINCS direct signature와 MCF7 recovered signature를 함께 사용
- direct/recovered/unknown availability flag를 모델에 명시

확인된 결과:

- `sample_features`: `16 x 5,360`
- `drug_features`: `243 x 1,563`
- `X_numeric_strong_context_smiles`: `3,387 x 6,655`
- `Numeric + Strong Context + SMILES / LightGBM` random3 Spearman `0.8835`
- `Numeric + Strong Context + SMILES / ExtraTrees` GroupCV Spearman `0.5112`

### v4 view-split ensemble

결측 축별 모델을 분리한다.

| View | Training rows | 목적 |
|---|---|---|
| all-view | SMILES-valid 전체 `3,387` rows | primary ranking 안정성 |
| LINCS-view | LINCS direct/imputed confidence high subset | drug perturbation signal 강화 |
| CRISPR-view | CRISPR-rich sample subset | dependency biology signal 강화 |
| context-view | target/pathway/KG context 중심 | evidence-aware reranking |

최종 점수:

```text
final_score =
  primary_model_score
  + lincs_view_bonus
  + crispr_view_bonus
  + target_kg_bonus
  - missingness_penalty
```

## 다음 실행 순서

1. `v2_smiles_required` 모델 재학습 완료
2. v1 vs v2 성능/후보 변화 비교 완료
3. source staging에서 LINCS/DepMap/DrugBank/ChEMBL 추가 매칭 가능성 스캔 완료
4. `reports/missingness/*_review.csv` 생성 완료
5. LINCS local bridge 후보 16개를 수동/규칙 기반으로 검토
6. DepMap expression/CNV/mutation fallback source를 붙여 CRISPR missing 7개 cell line 보강
7. 보강 가능한 source만 `data/raw/*_enhanced.csv`로 반영
8. v3 재학습
9. view-split ensemble 실험
