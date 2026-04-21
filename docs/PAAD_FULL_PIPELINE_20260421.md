# PAAD Full Drug Repurposing Pipeline - 2026-04-21

## Run scope

- Cancer type: pancreatic cancer / TCGA-PAAD
- Primary input set: `numeric_strong_context_smiles_pan_lincs`
- CV checks: random sample 4-fold and drug GroupCV 4-fold
- Model families: ML, DL, and ML+DL ensembles
- Training label: GDSC2 PAAD `LN_IC50`
- Source prefix: `s3://say2-4team/PAAD_raw/`

## Input data

- Response rows after SMILES filter: `6219`
- PAAD cell lines: `29`
- Drugs: `243`
- Feature matrix shape: `[6219, 2534]`
- LINCS numeric features: `513`

## Random sample 4-fold ensemble

| ensemble | spearman | pearson | rmse | mae | r2 | ndcg_at_20 |
| --- | --- | --- | --- | --- | --- | --- |
| top4_spearman_weighted | 0.9911 | 0.9902 | 0.3926 | 0.2014 | 0.9799 | 0.9699 |
| spearman_weighted | 0.9883 | 0.9884 | 0.4574 | 0.2965 | 0.9727 | 0.9791 |
| equal_weight | 0.9882 | 0.9884 | 0.4590 | 0.2980 | 0.9725 | 0.9791 |

## Drug GroupCV 4-fold ensemble

| ensemble | spearman | pearson | rmse | mae | r2 | ndcg_at_20 |
| --- | --- | --- | --- | --- | --- | --- |
| top4_spearman_weighted | 0.7512 | 0.7432 | 1.8676 | 1.3923 | 0.5453 | 0.8471 |
| spearman_weighted | 0.7246 | 0.7233 | 1.9162 | 1.4212 | 0.5213 | 0.8558 |
| equal_weight | 0.7170 | 0.7176 | 1.9314 | 1.4345 | 0.5137 | 0.8558 |

## Random sample member models

| member | spearman | pearson | rmse | mae | r2 | ndcg_at_20 |
| --- | --- | --- | --- | --- | --- | --- |
| DL_CrossAttention | 0.9920 | 0.9913 | 0.3702 | 0.2129 | 0.9821 | 0.9837 |
| ML_LightGBM | 0.9898 | 0.9893 | 0.4082 | 0.2351 | 0.9783 | 0.9823 |
| ML_ExtraTrees | 0.9862 | 0.9843 | 0.4987 | 0.2482 | 0.9676 | 0.9486 |
| ML_RandomForest | 0.9845 | 0.9834 | 0.5156 | 0.2849 | 0.9653 | 0.9537 |
| DL_WideDeep | 0.9790 | 0.9812 | 0.5348 | 0.3459 | 0.9627 | 0.9813 |
| ML_XGBoost | 0.9762 | 0.9816 | 0.5714 | 0.4169 | 0.9574 | 0.9746 |
| DL_ResidualMLP | 0.9732 | 0.9745 | 0.6333 | 0.4352 | 0.9477 | 0.9595 |
| ML_LightGBM_DART | 0.9700 | 0.9743 | 1.1630 | 0.9427 | 0.8237 | 0.9416 |
| ML_FlatMLP | 0.9455 | 0.9536 | 0.8358 | 0.6024 | 0.9089 | 0.9407 |

## Drug GroupCV member models

| member | spearman | pearson | rmse | mae | r2 | ndcg_at_20 |
| --- | --- | --- | --- | --- | --- | --- |
| ML_LightGBM_DART | 0.7398 | 0.7262 | 2.0616 | 1.6028 | 0.4459 | 0.8491 |
| ML_LightGBM | 0.7307 | 0.7229 | 1.9143 | 1.3875 | 0.5222 | 0.8511 |
| ML_XGBoost | 0.7285 | 0.7236 | 1.9115 | 1.4088 | 0.5236 | 0.8477 |
| DL_CrossAttention | 0.7132 | 0.7192 | 1.9368 | 1.4441 | 0.5110 | 0.8816 |
| ML_RandomForest | 0.7105 | 0.6827 | 2.0456 | 1.4768 | 0.4545 | 0.8499 |
| ML_ExtraTrees | 0.6970 | 0.6950 | 1.9922 | 1.4544 | 0.4826 | 0.8524 |
| DL_ResidualMLP | 0.5655 | 0.6113 | 2.2705 | 1.6877 | 0.3279 | 0.8382 |
| ML_FlatMLP | 0.5496 | 0.5837 | 2.3734 | 1.7790 | 0.2656 | 0.7702 |
| DL_WideDeep | 0.5289 | 0.5707 | 2.4232 | 1.8420 | 0.2345 | 0.7879 |

## GroupCV ensemble Top 15

| rank | drug_name | ensemble_score | mean_pred_ln_ic50 | screened_rows | target_genes | PATHWAY_NAME_NORMALIZED | classification |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1.0000 | Dactinomycin | 3.4108 | -3.4108 | 24.0000 | TOP1;TOP2A | Other | screened_candidate |
| 2.0000 | Docetaxel | 3.0918 | -3.0918 | 24.0000 | BCL2;NR1I2;TUBB1 | Mitosis | paad_standard_or_known |
| 3.0000 | MG-132 | 1.7314 | -1.7314 | 29.0000 | CAPN1;PROTEASOME | Protein stability and degradation | screened_candidate |
| 4.0000 | Topotecan | 1.0373 | -1.0373 | 24.0000 | TOP1;TOP1MT | DNA replication | paad_standard_or_known |
| 5.0000 | Vinblastine | 0.3808 | -0.3808 | 28.0000 | JUN;TUBA1A;TUBB;TUBB2A;TUBD1;TUBE1;TUBG1 | Mitosis | paad_trial_or_indication |
| 6.0000 | Teniposide | -0.2337 | 0.2337 | 24.0000 | TOP2A;TOP2B | DNA replication | screened_candidate |
| 7.0000 | Epirubicin | -0.3472 | 0.3472 | 29.0000 | ANTHRACYCLINE;TOP2A;TOP2B | DNA replication | paad_trial_or_indication |
| 8.0000 | Mitoxantrone | -0.4549 | 0.4549 | 24.0000 | TOP2;TOP2A;TOP2B | DNA replication | paad_trial_or_indication |
| 9.0000 | Lestaurtinib | -0.5889 | 0.5889 | 27.0000 | FLT3;JAK2;NTRK1;NTRK2;NTRK3 | Other kinases | paad_trial_or_indication |
| 10.0000 | OSI-027 | -0.5906 | 0.5906 | 29.0000 | MTOR | PI3K/MTOR signaling | screened_candidate |
| 11.0000 | Bortezomib | -0.7052 | 0.7052 | 29.0000 | PROTEASOME;PRSS1;PSMB1;PSMB5 | Protein stability and degradation | paad_trial_or_indication |
| 12.0000 | BI-2536 | -0.8184 | 0.8184 | 27.0000 | PLK1;PLK2;PLK3 | Cell cycle | screened_candidate |
| 13.0000 | Vinorelbine | -0.8925 | 0.8925 | 29.0000 | BCL2;MAP4;MAPK1;TUBB | Mitosis | paad_trial_or_indication |
| 14.0000 | BMS-345541 | -0.9423 | 0.9423 | 25.0000 | IKK-1;IKK-2 | Other | screened_candidate |
| 15.0000 | Staurosporine | -1.0570 | 1.0570 | 29.0000 | CDK2;CHRM1;CSK;GSK3B;IKBKB;ITK;LCK;MAPKAPK2;PDPK1;PIK3CG;PIM1;PRKCA;PRKCQ;SYK;ZAP70 | RTK signaling | paad_trial_or_indication |

## External validation snapshot

| rank | drug_name | target_match_genes | target_expression_pct | target_expressed | survival_p_value | clinical_trial_mention_count | prism_mean_log2fc | opentargets_max_score |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1.0000 | Dactinomycin | TOP1;TOP2A | 0.9344 | 1.0000 | 0.5066 | 0.0000 | -5.1549 | 0.4620 |
| 2.0000 | Docetaxel | BCL2;NR1I2;TUBB1 | 0.0000 | 0.0000 | 0.5166 | 849.0000 | -1.8692 | 0.4049 |
| 3.0000 | MG-132 | CAPN1 | 1.0000 | 1.0000 | 0.1763 | 0.0000 | -3.2282 | 0.0148 |
| 4.0000 | Topotecan | TOP1;TOP1MT | 0.9454 | 1.0000 | 0.4739 | 58.0000 | -1.5127 | 0.4620 |
| 5.0000 | Vinblastine | JUN;TUBA1A;TUBB;TUBB2A;TUBD1;TUBE1;TUBG1 | 0.9836 | 1.0000 | 0.4226 | 11.0000 | -0.3885 | 0.4049 |
| 6.0000 | Teniposide | TOP2A;TOP2B | 0.9235 | 1.0000 | 0.3792 | 0.0000 | -2.4304 | 0.2513 |
| 7.0000 | Epirubicin | TOP2A;TOP2B | 0.9235 | 1.0000 | 0.3792 | 56.0000 | -3.4320 | 0.2513 |
| 8.0000 | Mitoxantrone | TOP2A;TOP2B | 0.9235 | 1.0000 | 0.3792 | 12.0000 | -2.9996 | 0.2513 |
| 9.0000 | Lestaurtinib | FLT3;JAK2;NTRK1;NTRK2;NTRK3 | 0.0000 | 0.0000 | 0.8770 | 0.0000 | -2.8677 | 0.4643 |
| 10.0000 | OSI-027 | MTOR | 0.3169 | 1.0000 | 0.7752 | 0.0000 | -0.3381 | 0.3188 |
| 11.0000 | Bortezomib | PRSS1;PSMB1;PSMB5 | 1.0000 | 1.0000 | 0.6895 | 182.0000 | -4.0230 | 0.8075 |
| 12.0000 | BI-2536 | PLK1;PLK2;PLK3 | 0.8361 | 1.0000 | 0.0149 | 0.0000 | -2.0015 | 0.1118 |
| 13.0000 | Vinorelbine | BCL2;MAP4;MAPK1;TUBB | 0.9836 | 1.0000 | 0.1227 | 26.0000 | -1.7224 | 0.4049 |
| 14.0000 | BMS-345541 |  |  | 0.0000 |  | 0.0000 | -0.1259 |  |
| 15.0000 | Staurosporine | CDK2;CHRM1;CSK;GSK3B;IKBKB;ITK;LCK;MAPKAPK2;PDPK1;PIK3CG;PIM1;PRKCA;PRKCQ;SYK;ZAP70 | 0.1475 | 0.0000 | 0.6304 | 27.0000 | -3.3232 | 0.3047 |

## Final tiered candidates

| drug_name | tier | knowledge_score | final_category | admet_category | ensemble_score | target_expressed | clinical_trial_mention_count | prism_mean_log2fc |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Topotecan | Tier 1 | 4.1684 | Known pancreatic cancer control/support | Candidate | 1.0373 | 1.0000 | 58.0000 | -1.5127 |
| Vinblastine | Tier 1 | 3.6825 | Pancreatic clinical-trial supported candidate | Approved | 0.3808 | 1.0000 | 11.0000 | -0.3885 |
| Vinorelbine | Tier 1 | 3.6162 | Pancreatic clinical-trial supported candidate | Candidate | -0.8925 | 1.0000 | 26.0000 | -1.7224 |
| MG-132 | Tier 1 | 3.0314 | Exploratory repurposing candidate | Candidate | 1.7314 | 1.0000 | 0.0000 | -3.2282 |
| BI-2536 | Tier 1 | 3.0285 | Exploratory repurposing candidate | Candidate | -0.8184 | 1.0000 | 0.0000 | -2.0015 |
| Bortezomib | Tier 2 | 3.4908 | Pancreatic clinical-trial supported candidate | Caution | -0.7052 | 1.0000 | 182.0000 | -4.0230 |
| Epirubicin | Tier 2 | 3.2013 | Pancreatic clinical-trial supported candidate | Caution | -0.3472 | 1.0000 | 56.0000 | -3.4320 |
| Mitoxantrone | Tier 2 | 3.1347 | Pancreatic clinical-trial supported candidate | Caution | -0.4549 | 1.0000 | 12.0000 | -2.9996 |
| Docetaxel | Tier 2 | 2.8729 | Known pancreatic cancer control/support | Caution | 3.0918 | 0.0000 | 849.0000 | -1.8692 |
| Dactinomycin | Tier 2 | 2.8120 | Exploratory repurposing candidate | Caution | 3.4108 | 1.0000 | 0.0000 | -5.1549 |
| Staurosporine | Tier 2 | 2.7713 | Pancreatic clinical-trial supported candidate | Candidate | -1.0570 | 0.0000 | 27.0000 | -3.3232 |
| Lestaurtinib | Tier 2 | 2.3310 | Exploratory repurposing candidate | Candidate | -0.5889 | 0.0000 | 0.0000 | -2.8677 |
| Teniposide | Tier 2 | 2.2680 | Exploratory repurposing candidate | Caution | -0.2337 | 1.0000 | 0.0000 | -2.4304 |
| OSI-027 | Tier 3 | 2.0378 | Exploratory repurposing candidate | Candidate | -0.5906 | 1.0000 | 0.0000 | -0.3381 |
| BMS-345541 | Tier 3 | 0.5963 | Exploratory repurposing candidate | Candidate | -0.9423 | 0.0000 | 0.0000 | -0.1259 |

## Key outputs

- `results/paad/ml/random4/numeric_strong_context_smiles_pan_lincs/`
- `results/paad/ml/groupcv4_drug/numeric_strong_context_smiles_pan_lincs/`
- `results/paad/dl/random4/numeric_strong_context_smiles_pan_lincs/`
- `results/paad/dl/groupcv4_drug/numeric_strong_context_smiles_pan_lincs/`
- `results/paad/ensemble/random4/numeric_strong_context_smiles_pan_lincs/`
- `results/paad/ensemble/groupcv4_drug/numeric_strong_context_smiles_pan_lincs/`
- `external_validation/paad/groupcv4_drug/top50_external_validation.csv`
- `phase5_final_results/paad/groupcv4_drug/final_comprehensive_candidates.csv`
