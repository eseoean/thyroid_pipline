# 갑상선암 약물 재창출 Hybrid Pipeline 최종 리포트

- 프로젝트: thyroid_pipline
- 질환: thyroid cancer (THCA)
- 기준: Choi numeric base + say2 SMILES/strong context/random3 ensemble + Choi external/ADMET/KG validation
- 최종 추천은 screened drug-response 데이터에 한정한다.

## 1. 실행 요약

- Numeric policy: choi_5524_direction
- SMILES dimension: 64
- Strong context dimension: 32
- CV: random sample 3-fold OOF

## 2. 모델 성능

| input_set | model | spearman | pearson | rmse | mae | r2 | ndcg_at_20 | train_oof_spearman_gap |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| numeric | LightGBM | 0.8776 | 0.9165 | 1.1587 | 0.8854 | 0.8392 | 0.9692 | 0.0758 |
| numeric | RandomForest | 0.868 | 0.9078 | 1.2172 | 0.9165 | 0.8226 | 0.9685 | 0.1046 |
| numeric | XGBoost | 0.8617 | 0.9077 | 1.2398 | 0.9749 | 0.816 | 0.9646 | 0.0489 |
| numeric | ExtraTrees | 0.8587 | 0.9023 | 1.2459 | 0.9191 | 0.8141 | 0.9602 | 0.1335 |
| numeric | LightGBM_DART | 0.8483 | 0.8936 | 1.5652 | 1.2515 | 0.7067 | 0.9636 | 0.0569 |
| numeric | FlatMLP | 0.8172 | 0.8591 | 1.4821 | 1.1417 | 0.737 | 0.911 | 0.1505 |
| numeric_smiles | LightGBM | 0.8833 | 0.9199 | 1.1353 | 0.863 | 0.8457 | 0.9654 | 0.079 |
| numeric_smiles | XGBoost | 0.8765 | 0.9146 | 1.1895 | 0.9226 | 0.8306 | 0.9688 | 0.0487 |
| numeric_smiles | RandomForest | 0.8765 | 0.9134 | 1.1815 | 0.8937 | 0.8328 | 0.968 | 0.0969 |
| numeric_smiles | LightGBM_DART | 0.8641 | 0.904 | 1.5133 | 1.2102 | 0.7258 | 0.9609 | 0.0561 |
| numeric_smiles | ExtraTrees | 0.8629 | 0.9044 | 1.2335 | 0.9108 | 0.8178 | 0.9599 | 0.1297 |
| numeric_smiles | FlatMLP | 0.8253 | 0.8634 | 1.4597 | 1.1131 | 0.7448 | 0.9122 | 0.1388 |
| numeric_strong_context_smiles | LightGBM | 0.8835 | 0.9212 | 1.1264 | 0.8609 | 0.8481 | 0.9689 | 0.0789 |
| numeric_strong_context_smiles | XGBoost | 0.8767 | 0.9158 | 1.1825 | 0.9188 | 0.8326 | 0.9663 | 0.0481 |
| numeric_strong_context_smiles | RandomForest | 0.8763 | 0.9132 | 1.1826 | 0.8947 | 0.8325 | 0.9679 | 0.0972 |
| numeric_strong_context_smiles | LightGBM_DART | 0.8638 | 0.9048 | 1.512 | 1.2094 | 0.7263 | 0.9599 | 0.0563 |
| numeric_strong_context_smiles | ExtraTrees | 0.8635 | 0.905 | 1.2295 | 0.908 | 0.819 | 0.9651 | 0.1292 |
| numeric_strong_context_smiles | FlatMLP | 0.8287 | 0.8677 | 1.4382 | 1.0997 | 0.7523 | 0.9223 | 0.1262 |

## 3. 앙상블 및 Diversity

- 앙상블 OOF Spearman: 0.8821965117740939
- 앙상블 OOF RMSE: 1.1713416576385498
- 모델 가중치: {'LightGBM': 0.17014973972721734, 'XGBoost': 0.16883644801315215, 'RandomForest': 0.16876126633429683, 'LightGBM_DART': 0.1663534145925931, 'ExtraTrees': 0.16629326106133865, 'FlatMLP': 0.15960587027140197}

## 4. 최종 후보

| drug_name | canonical_drug_id | target_genes | pathway | ensemble_score | admet_category | drug_target_evidence | target_thyroid_relevance | clinical_evidence | mechanism_rationale | safety_profile | target_expression_bonus | knowledge_score | tier | final_category |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Vinorelbine | 2048 | BCL2;MAP4;MAPK1;TUBB | Mitosis | 2.430091619491577 | Candidate | 1.0 | 1.0 | 0.4 | 1.0 | 0.6 | 0.5 | 4.5 | Tier 1 | Thyroid indication expansion candidate |
| Staurosporine | 1034 | CDK2;CHRM1;CSK;GSK3B;IKBKB;ITK;LCK;MAPKAPK2;PDPK1;PIK3CG;PIM1;PRKCA;PRKCQ;SYK;ZAP70 | RTK signaling | 2.650351047515869 | Caution | 1.0 | 1.0 | 0.4 | 1.0 | 0.2 | 0.0 | 3.6 | Tier 2 | True repurposing or exploratory candidate |
| Vinblastine | 1004 | JUN;TUBA1A;TUBB;TUBB2A;TUBD1;TUBE1;TUBG1 | Mitosis | 2.6648731231689453 | Approved | 1.0 | 0.0 | 0.4 | 0.3 | 1.0 | 0.5 | 3.2 | Tier 2 | True repurposing or exploratory candidate |
| Romidepsin | 1817 | HDAC1;HDAC2;HDAC3;HDAC4;HDAC6;HDAC8 | Chromatin histone acetylation | 2.5289530754089355 | Candidate | 1.0 | 0.0 | 0.4 | 0.3 | 0.6 | 0.5 | 2.8 | Tier 3 | Thyroid indication expansion candidate |
| Bortezomib | 1191 | PROTEASOME;PRSS1;PSMB1;PSMB5 | Protein stability and degradation | 4.065159797668457 | Caution | 1.0 | 0.0 | 0.4 | 0.3 | 0.2 | 0.5 | 2.4 | Excluded | Thyroid indication expansion candidate |
| Dactinomycin | 1911 | TOP1;TOP2A | Other | 3.944403648376465 | Caution | 1.0 | 0.0 | 0.4 | 0.3 | 0.2 | 0.5 | 2.4 | Excluded | True repurposing or exploratory candidate |
| SN-38 | 1494 | TOP1 | DNA replication | 2.378016471862793 | Caution | 1.0 | 0.0 | 0.4 | 0.3 | 0.2 | 0.5 | 2.4 | Excluded | True repurposing or exploratory candidate |
| Daporinad | 1248 | NAMPT | Metabolism | 2.1311774253845215 | Caution | 1.0 | 0.0 | 0.4 | 0.3 | 0.2 | 0.5 | 2.4 | Excluded | True repurposing or exploratory candidate |
| Camptothecin | 1003 | TOP1 | DNA replication | 1.9975780248641968 | Caution | 1.0 | 0.0 | 0.4 | 0.3 | 0.2 | 0.5 | 2.4 | Excluded | Thyroid indication expansion candidate |
| MG-132 | 1862 | CAPN1;PROTEASOME | Protein stability and degradation | 1.3018019199371338 | Caution | 1.0 | 0.0 | 0.4 | 0.3 | 0.2 | 0.5 | 2.4 | Excluded | True repurposing or exploratory candidate |
| Dinaciclib | 1180 | CDK1;CDK2;CDK5;CDK9 | Cell cycle | 1.2586849927902222 | Caution | 1.0 | 0.0 | 0.4 | 0.3 | 0.2 | 0.5 | 2.4 | Excluded | True repurposing or exploratory candidate |
| Docetaxel | 1007 | BCL2;NR1I2;TUBB1 | Mitosis | 3.530460834503174 | Caution | 1.0 | 0.0 | 0.4 | 0.3 | 0.2 | 0.0 | 1.9 | Excluded | Thyroid indication expansion candidate |
| Sepantronium bromide | 1941 | BIRC5 | Apoptosis regulation | 2.7170848846435547 | Caution | 1.0 | 0.0 | 0.4 | 0.3 | 0.2 | 0.0 | 1.9 | Excluded | True repurposing or exploratory candidate |
| Paclitaxel | 1080 | BCL2;NR1I2;TUBB1 | Mitosis | 1.797881245613098 | Caution | 1.0 | 0.0 | 0.4 | 0.3 | 0.2 | 0.0 | 1.9 | Excluded | Thyroid indication expansion candidate |
| Luminespib | 1559 | HSP90 | Protein stability and degradation | 1.57468581199646 | Caution | 1.0 | 0.0 | 0.4 | 0.3 | 0.2 | 0.0 | 1.9 | Excluded | True repurposing or exploratory candidate |

## 5. QC 요약

- 각 단계별 상세 QC는 `reports/qc_step*.json`에 저장했다.
- 외부검증 결과는 `external_validation/`, ADMET 결과는 `admet/`, knowledge validation 결과는 `phase5_final_results/`에 저장했다.

## 6. 산출물 경로

- `data/X_numeric.npy`
- `data/X_numeric_smiles.npy`
- `data/X_numeric_strong_context_smiles.npy`
- `results/random3/`
- `results/ensemble/thyroid_ensemble_top30_drugs.csv`
- `external_validation/top15_validated.csv`
- `admet/final_drug_candidates.csv`
- `phase5_final_results/final_comprehensive_candidates.csv`
