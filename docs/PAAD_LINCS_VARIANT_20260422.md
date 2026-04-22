# PAAD-Focused LINCS Variant - 2026-04-22

## Scope

- Cancer type: pancreatic cancer / TCGA-PAAD
- Input set: `numeric_strong_context_smiles_paad_lincs`
- Baseline comparator: `numeric_strong_context_smiles_pan_lincs`
- CV checks: random sample 4-fold and drug GroupCV 4-fold
- Model families: ML, DL, and ML+DL ensembles

## Exact PAAD LINCS Availability

- Pancreas-like LINCS cells in metadata: `['YAPC', 'YAPC.311']`
- `sig_info` rows for those cells: `0`
- `sig_info` trt_cp rows for those cells: `0`
- `inst_info` rows for those cells: `0`
- Exact PAAD-cell LINCS available: `False`

Because exact pancreatic-cell LINCS signatures were absent in the available source metadata, this run uses a PAAD-focused LINCS feature selection:
pan-cancer drug perturbation values are retained only for LINCS genes supported by PAAD priority biology, candidate drug targets, TCGA-PAAD expression coverage, or PAAD DepMap feature availability.

## Input QC

- Source shape: `[6219, 2534]`
- Output shape: `[6219, 2110]`
- Source LINCS gene features: `512`
- Selected PAAD-focused LINCS gene features: `88`
- Kept LINCS flag features: `1`
- Selected feature table: `/Users/skku_aws2_18/team4_project/pre_project/thyroid_pipline/reports/paad/paad_lincs_selected_features_20260422.csv`

## Variant Comparison

| variant | cv | ensemble | spearman | pearson | rmse | mae | r2 | ndcg_at_20 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| pan_cancer_lincs | random4 | top4_spearman_weighted | 0.9911 | 0.9902 | 0.3926 | 0.2014 | 0.9799 | 0.9699 |
| paad_focused_lincs | random4 | top4_spearman_weighted | 0.9910 | 0.9899 | 0.3991 | 0.2047 | 0.9792 | 0.9695 |
| pan_cancer_lincs | groupcv4_drug | top4_spearman_weighted | 0.7512 | 0.7432 | 1.8676 | 1.3923 | 0.5453 | 0.8471 |
| paad_focused_lincs | groupcv4_drug | top4_spearman_weighted | 0.7374 | 0.7199 | 1.9292 | 1.4338 | 0.5148 | 0.8455 |

## Random sample 4-fold ensemble

| ensemble | spearman | pearson | rmse | mae | r2 | ndcg_at_20 |
| --- | --- | --- | --- | --- | --- | --- |
| top4_spearman_weighted | 0.9910 | 0.9899 | 0.3991 | 0.2047 | 0.9792 | 0.9695 |
| spearman_weighted | 0.9887 | 0.9886 | 0.4521 | 0.2917 | 0.9734 | 0.9747 |
| equal_weight | 0.9886 | 0.9886 | 0.4535 | 0.2933 | 0.9732 | 0.9747 |

## Drug GroupCV 4-fold ensemble

| ensemble | spearman | pearson | rmse | mae | r2 | ndcg_at_20 |
| --- | --- | --- | --- | --- | --- | --- |
| top4_spearman_weighted | 0.7374 | 0.7199 | 1.9292 | 1.4338 | 0.5148 | 0.8455 |
| spearman_weighted | 0.7198 | 0.7193 | 1.9274 | 1.4358 | 0.5157 | 0.8539 |
| equal_weight | 0.7082 | 0.7119 | 1.9468 | 1.4556 | 0.5059 | 0.8516 |

## Random sample member models

| member | spearman | pearson | rmse | mae | r2 | ndcg_at_20 |
| --- | --- | --- | --- | --- | --- | --- |
| DL_CrossAttention | 0.9919 | 0.9912 | 0.3749 | 0.2156 | 0.9817 | 0.9801 |
| ML_LightGBM | 0.9893 | 0.9888 | 0.4172 | 0.2433 | 0.9773 | 0.9752 |
| ML_ExtraTrees | 0.9862 | 0.9841 | 0.5017 | 0.2503 | 0.9672 | 0.9526 |
| ML_RandomForest | 0.9840 | 0.9820 | 0.5378 | 0.2972 | 0.9623 | 0.9530 |
| DL_WideDeep | 0.9820 | 0.9835 | 0.5018 | 0.3210 | 0.9672 | 0.9746 |
| DL_ResidualMLP | 0.9767 | 0.9772 | 0.6013 | 0.4080 | 0.9529 | 0.9671 |
| ML_XGBoost | 0.9752 | 0.9807 | 0.5862 | 0.4334 | 0.9552 | 0.9753 |
| ML_LightGBM_DART | 0.9679 | 0.9731 | 1.1697 | 0.9487 | 0.8216 | 0.9454 |
| ML_FlatMLP | 0.9475 | 0.9572 | 0.8033 | 0.5813 | 0.9159 | 0.9432 |

## Drug GroupCV member models

| member | spearman | pearson | rmse | mae | r2 | ndcg_at_20 |
| --- | --- | --- | --- | --- | --- | --- |
| ML_LightGBM_DART | 0.7387 | 0.7197 | 2.0669 | 1.6340 | 0.4430 | 0.8490 |
| ML_LightGBM | 0.7314 | 0.7190 | 1.9254 | 1.4129 | 0.5167 | 0.8540 |
| ML_XGBoost | 0.7284 | 0.7146 | 1.9385 | 1.4294 | 0.5101 | 0.8488 |
| ML_RandomForest | 0.7163 | 0.6860 | 2.0259 | 1.4593 | 0.4649 | 0.8596 |
| DL_CrossAttention | 0.7115 | 0.7191 | 1.9370 | 1.3961 | 0.5109 | 0.9239 |
| ML_ExtraTrees | 0.6828 | 0.6783 | 2.0370 | 1.4911 | 0.4590 | 0.8527 |
| DL_ResidualMLP | 0.5365 | 0.5957 | 2.3090 | 1.7583 | 0.3049 | 0.8091 |
| ML_FlatMLP | 0.5211 | 0.5650 | 2.4156 | 1.8695 | 0.2392 | 0.8595 |
| DL_WideDeep | 0.4739 | 0.5431 | 2.4257 | 1.8641 | 0.2329 | 0.7742 |

## GroupCV ensemble Top 15

| rank | drug_name | ensemble_score | mean_pred_ln_ic50 | screened_rows | target_genes | PATHWAY_NAME_NORMALIZED | classification |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1.0000 | Dactinomycin | 3.3077 | -3.3077 | 24.0000 | TOP1;TOP2A | Other | screened_candidate |
| 2.0000 | Docetaxel | 3.0993 | -3.0993 | 24.0000 | BCL2;NR1I2;TUBB1 | Mitosis | paad_standard_or_known |
| 3.0000 | Topotecan | 1.4081 | -1.4081 | 24.0000 | TOP1;TOP1MT | DNA replication | paad_standard_or_known |
| 4.0000 | Temsirolimus | 0.5911 | -0.5911 | 26.0000 | MTOR | PI3K/MTOR signaling | paad_trial_or_indication |
| 5.0000 | Staurosporine | 0.0608 | -0.0608 | 29.0000 | CDK2;CHRM1;CSK;GSK3B;IKBKB;ITK;LCK;MAPKAPK2;PDPK1;PIK3CG;PIM1;PRKCA;PRKCQ;SYK;ZAP70 | RTK signaling | paad_trial_or_indication |
| 6.0000 | Lestaurtinib | -0.1837 | 0.1837 | 27.0000 | FLT3;JAK2;NTRK1;NTRK2;NTRK3 | Other kinases | paad_trial_or_indication |
| 7.0000 | MG-132 | -0.5094 | 0.5094 | 29.0000 | CAPN1;PROTEASOME | Protein stability and degradation | screened_candidate |
| 8.0000 | Vinorelbine | -0.5181 | 0.5181 | 29.0000 | BCL2;MAP4;MAPK1;TUBB | Mitosis | paad_trial_or_indication |
| 9.0000 | Vinblastine | -0.7537 | 0.7537 | 28.0000 | JUN;TUBA1A;TUBB;TUBB2A;TUBD1;TUBE1;TUBG1 | Mitosis | paad_trial_or_indication |
| 10.0000 | Teniposide | -0.7940 | 0.7940 | 24.0000 | TOP2A;TOP2B | DNA replication | screened_candidate |
| 11.0000 | Epirubicin | -1.0386 | 1.0386 | 29.0000 | ANTHRACYCLINE;TOP2A;TOP2B | DNA replication | paad_trial_or_indication |
| 12.0000 | Selumetinib | -1.2170 | 1.2170 | 24.0000 | MAP2K1;MAP2K2 | ERK MAPK signaling | paad_trial_or_indication |
| 13.0000 | Dasatinib | -1.2613 | 1.2613 | 29.0000 | ABL;ABL1;ABL2;BCR;BTK;CSK;EPHA2;EPHA5;EPHB4;EPHRINS;FGR;FRK;FYN;HSPA8;KIT;LCK;LYN;MAP3K20;MAPK14;NR4A3;PDGFR;PDGFRB;PPAT;SRC;STAT5B;YES1 | Other kinases | paad_trial_or_indication |
| 14.0000 | Dactolisib | -1.3083 | 1.3083 | 29.0000 | MTOR;PIK3CG | PI3K/MTOR signaling | paad_trial_or_indication |
| 15.0000 | Camptothecin | -1.3409 | 1.3409 | 29.0000 | TOP1 | DNA replication | paad_trial_or_indication |

## Key outputs

- `data/paad/X_numeric_strong_context_smiles_paad_lincs.npy`
- `data/paad/numeric_strong_context_smiles_paad_lincs_feature_names.json`
- `reports/paad/qc_paad_lincs_variant_20260422.json`
- `reports/paad/paad_lincs_selected_features_20260422.csv`
- `results/paad/ml/random4/numeric_strong_context_smiles_paad_lincs/`
- `results/paad/ml/groupcv4_drug/numeric_strong_context_smiles_paad_lincs/`
- `results/paad/dl/random4/numeric_strong_context_smiles_paad_lincs/`
- `results/paad/dl/groupcv4_drug/numeric_strong_context_smiles_paad_lincs/`
- `results/paad/ensemble/random4/numeric_strong_context_smiles_paad_lincs/`
- `results/paad/ensemble/groupcv4_drug/numeric_strong_context_smiles_paad_lincs/`
