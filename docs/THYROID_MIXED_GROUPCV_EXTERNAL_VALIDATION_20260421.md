# Thyroid Mixed GroupCV External Validation - 2026-04-21

## 목적

최신 모델 후보인 `pan-cancer LINCS + CrossAttention + ResidualMLP + LightGBM` GroupCV 앙상블 Top30을 기준으로 외부검증, ADMET, KG/Tier 검증을 다시 수행했다.

## 입력 후보

- Variant: `mixed_groupcv_pan_lincs`
- Source: `/Users/skku_aws2_18/team4_project/pre_project/thyroid_pipline/results/mixed_groupcv/crossattention_residualmlp_lightgbm_groupcv/mixed_groupcv_ensemble_top30_drugs.csv`
- Model candidate source: `CrossAttention + ResidualMLP + LightGBM`
- Input features: `numeric + strong context + smiles + pan-cancer LINCS`
- CV basis: `GroupCV by canonical_drug_id`

## 최신 앙상블 Top15

| rank | drug_name | ensemble_score | mean_pred_ln_ic50 | target_genes | PATHWAY_NAME_NORMALIZED | classification |
| --- | --- | --- | --- | --- | --- | --- |
| 1.0000 | Dactinomycin | 4.3256 | -4.3256 | TOP1;TOP2A | Other | screened_candidate |
| 2.0000 | Docetaxel | 3.3417 | -3.3417 | BCL2;NR1I2;TUBB1 | Mitosis | indication_expansion |
| 3.0000 | Topotecan | 2.8960 | -2.8960 | TOP1;TOP1MT | DNA replication | screened_candidate |
| 4.0000 | Vinblastine | 2.3322 | -2.3322 | JUN;TUBA1A;TUBB;TUBB2A;TUBD1;TUBE1;TUBG1 | Mitosis | screened_candidate |
| 5.0000 | Epirubicin | 1.8480 | -1.8480 | ANTHRACYCLINE;TOP2A;TOP2B | DNA replication | indication_expansion |
| 6.0000 | Mitoxantrone | 1.6478 | -1.6478 | TOP2;TOP2A;TOP2B | DNA replication | indication_expansion |
| 7.0000 | Camptothecin | 0.8951 | -0.8951 | TOP1 | DNA replication | indication_expansion |
| 8.0000 | Lestaurtinib | 0.5204 | -0.5204 | FLT3;JAK2;NTRK1;NTRK2;NTRK3 | Other kinases | screened_candidate |
| 9.0000 | Staurosporine | 0.4478 | -0.4478 | CDK2;CHRM1;CSK;GSK3B;IKBKB;ITK;LCK;MAPKAPK2;PDPK1;PIK3CG;PIM1;PRKCA;PRKCQ;SYK;ZAP70 | RTK signaling | screened_candidate |
| 10.0000 | Piperlongumine | 0.2498 | -0.2498 |  | Other | screened_candidate |
| 11.0000 | MG-132 | 0.1759 | -0.1759 | CAPN1;PROTEASOME | Protein stability and degradation | screened_candidate |
| 12.0000 | Paclitaxel | 0.1459 | -0.1459 | BCL2;NR1I2;TUBB1 | Mitosis | indication_expansion |
| 13.0000 | BX795 | 0.0037 | -0.0037 | AURKB;AURKC;IKK;TBK1 | Other kinases | screened_candidate |
| 14.0000 | Teniposide | -0.0667 | 0.0667 | TOP2A;TOP2B | DNA replication | screened_candidate |
| 15.0000 | Vinorelbine | -0.1815 | 0.1815 | BCL2;MAP4;MAPK1;TUBB | Mitosis | indication_expansion |

## Step 8. TCGA-THCA 외부검증 요약

- Expression genes: `296`
- Patient count: `572`
- Mean target gene match rate: `0.9126`
- Target expressed in Top15: `9 / 15`
- Known thyroid positive-control precision: `{"p_at_5": 0.0, "p_at_10": 0.0, "p_at_20": 0.0}`

| rank | drug_name | target_match_genes | target_expression_pct | target_expressed | survival_p_value | survival_direction | known_thyroid_control |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1.0000 | Dactinomycin | TOP1;TOP2A | 0.7745 | 1.0000 | 0.6133 | high_expression_shorter_survival | 0.0000 |
| 2.0000 | Docetaxel | BCL2;NR1I2;TUBB1 | 0.0000 | 0.0000 | 0.9628 | high_expression_longer_survival | 0.0000 |
| 3.0000 | Topotecan | TOP1;TOP1MT | 0.9843 | 1.0000 | 0.2657 | high_expression_longer_survival | 0.0000 |
| 4.0000 | Vinblastine | JUN;TUBA1A;TUBB;TUBB2A;TUBD1;TUBE1;TUBG1 | 0.9860 | 1.0000 | 0.1362 | high_expression_longer_survival | 0.0000 |
| 5.0000 | Epirubicin | TOP2A;TOP2B | 0.7990 | 1.0000 | 0.4246 | high_expression_shorter_survival | 0.0000 |
| 6.0000 | Mitoxantrone | TOP2A;TOP2B | 0.7990 | 1.0000 | 0.4246 | high_expression_shorter_survival | 0.0000 |
| 7.0000 | Camptothecin | TOP1 | 0.9930 | 1.0000 | 0.6496 | high_expression_longer_survival | 0.0000 |
| 8.0000 | Lestaurtinib | FLT3;JAK2;NTRK1;NTRK2;NTRK3 | 0.0017 | 0.0000 | 0.9591 | high_expression_longer_survival | 0.0000 |
| 9.0000 | Staurosporine | CDK2;CHRM1;CSK;GSK3B;IKBKB;ITK;LCK;MAPKAPK2;PDPK1;PIK3CG;PIM1;PRKCA;PRKCQ;SYK;ZAP70 | 0.0647 | 0.0000 | 0.7976 | high_expression_longer_survival | 0.0000 |
| 10.0000 | Piperlongumine |  |  | 0.0000 |  | no_clinical_data | 0.0000 |
| 11.0000 | MG-132 | CAPN1 | 0.9983 | 1.0000 | 0.0839 | high_expression_longer_survival | 0.0000 |
| 12.0000 | Paclitaxel | BCL2;NR1I2;TUBB1 | 0.0000 | 0.0000 | 0.9628 | high_expression_longer_survival | 0.0000 |
| 13.0000 | BX795 | AURKB;AURKC;TBK1 | 0.0000 | 0.0000 | 0.2416 | high_expression_shorter_survival | 0.0000 |
| 14.0000 | Teniposide | TOP2A;TOP2B | 0.7990 | 1.0000 | 0.4246 | high_expression_shorter_survival | 0.0000 |
| 15.0000 | Vinorelbine | BCL2;MAP4;MAPK1;TUBB | 0.9965 | 1.0000 | 0.0680 | high_expression_longer_survival | 0.0000 |

## Step 9. ADMET 요약

- Candidate count: `15`
- Assay count: `22`
- Category counts: `{"Candidate": 8, "Caution": 6, "Approved": 1}`
- Toxicity flag count: `6`
- Low-confidence toxic signal count: `10`

| drug_name | target_expressed | admet_coverage | toxicity_flags | low_confidence_toxic_signals | admet_category |
| --- | --- | --- | --- | --- | --- |
| Dactinomycin | 1.0000 | 0.2273 | dili | herg | Caution |
| Docetaxel | 0.0000 | 0.5455 | dili | herg | Caution |
| Topotecan | 1.0000 | 0.4545 |  | ames | Candidate |
| Vinblastine | 1.0000 | 0.5455 |  |  | Approved |
| Epirubicin | 1.0000 | 0.6818 | ames;dili |  | Caution |
| Mitoxantrone | 1.0000 | 0.5455 | ames;dili |  | Caution |
| Camptothecin | 1.0000 | 0.2273 |  | ames | Candidate |
| Lestaurtinib | 0.0000 | 0.0000 |  | ames;dili;herg | Candidate |
| Staurosporine | 0.0000 | 0.0455 |  | ames;dili;herg | Candidate |
| Piperlongumine | 0.0000 | 0.1364 |  |  | Candidate |
| MG-132 | 1.0000 | 0.1364 |  | dili;herg | Candidate |
| Paclitaxel | 0.0000 | 0.5455 | dili | ames;herg | Caution |
| BX795 | 0.0000 | 0.0455 |  | dili;herg | Candidate |
| Teniposide | 1.0000 | 0.8182 | ames;dili | herg | Caution |
| Vinorelbine | 1.0000 | 0.4545 |  |  | Candidate |

## Step 10. KG/Tier 검증 요약

- Candidate count: `15`
- Tier counts: `{"Excluded": 6, "Tier 3": 5, "Tier 2": 3, "Tier 1": 1}`
- Category counts: `{"True repurposing or exploratory candidate": 9, "Thyroid indication expansion candidate": 6}`

| drug_name | ensemble_score | admet_category | knowledge_score | tier | final_category |
| --- | --- | --- | --- | --- | --- |
| Vinorelbine | -0.1815 | Candidate | 4.5000 | Tier 1 | Thyroid indication expansion candidate |
| Lestaurtinib | 0.5204 | Candidate | 4.0000 | Tier 2 | True repurposing or exploratory candidate |
| Staurosporine | 0.4478 | Candidate | 4.0000 | Tier 2 | True repurposing or exploratory candidate |
| Vinblastine | 2.3322 | Approved | 3.2000 | Tier 2 | True repurposing or exploratory candidate |
| Topotecan | 2.8960 | Candidate | 2.8000 | Tier 3 | True repurposing or exploratory candidate |
| Camptothecin | 0.8951 | Candidate | 2.8000 | Tier 3 | Thyroid indication expansion candidate |
| MG-132 | 0.1759 | Candidate | 2.8000 | Tier 3 | True repurposing or exploratory candidate |
| BX795 | 0.0037 | Candidate | 2.3000 | Tier 3 | True repurposing or exploratory candidate |
| Piperlongumine | 0.2498 | Candidate | 1.3000 | Tier 3 | True repurposing or exploratory candidate |
| Dactinomycin | 4.3256 | Caution | 2.4000 | Excluded | True repurposing or exploratory candidate |
| Epirubicin | 1.8480 | Caution | 2.4000 | Excluded | Thyroid indication expansion candidate |
| Mitoxantrone | 1.6478 | Caution | 2.4000 | Excluded | Thyroid indication expansion candidate |
| Teniposide | -0.0667 | Caution | 2.4000 | Excluded | True repurposing or exploratory candidate |
| Docetaxel | 3.3417 | Caution | 1.9000 | Excluded | Thyroid indication expansion candidate |
| Paclitaxel | 0.1459 | Caution | 1.9000 | Excluded | Thyroid indication expansion candidate |

## 산출물

- `external_validation/mixed_groupcv_pan_lincs/top15_validated.csv`
- `external_validation/mixed_groupcv_pan_lincs/thyroid_target_expression.csv`
- `external_validation/mixed_groupcv_pan_lincs/thyroid_survival_validation.csv`
- `admet/mixed_groupcv_pan_lincs/final_drug_candidates.csv`
- `knowledge_validation/mixed_groupcv_pan_lincs/validation_summary.csv`
- `phase5_final_results/mixed_groupcv_pan_lincs/final_comprehensive_candidates.csv`
- `phase5_final_results/mixed_groupcv_pan_lincs/tier1_high_confidence.csv`
- `reports/mixed_groupcv_pan_lincs/qc_step8_external_validation.json`
- `reports/mixed_groupcv_pan_lincs/qc_step9_admet.json`
- `reports/mixed_groupcv_pan_lincs/qc_step10_knowledge_validation.json`
