# PAAD YAPC LINCS Variant - 2026-04-22

## Scope

- Cancer type: pancreatic cancer / TCGA-PAAD
- LINCS source: GSE70138 Phase II Level5
- Cell filter: `YAPC`, `YAPC.311`
- Perturbation filter: `pert_type == trt_cp`
- Input set: `numeric_strong_context_smiles_yapc_lincs`
- Baseline comparator: `numeric_strong_context_smiles_pan_lincs`
- CV checks: random sample 4-fold and drug GroupCV 4-fold
- Model families: ML, DL, and ML+DL ensembles

## Signature QC

- Target cells: `['YAPC', 'YAPC.311']`
- YAPC/YAPC.311 signature rows: `11001`
- YAPC/YAPC.311 trt_cp signature rows: `10061`
- Unique YAPC perturbagens: `1568`
- Mapped PAAD GDSC drugs: `84`
- Mapped Level5 signatures: `1141`
- Selected YAPC LINCS features: `512`
- Signature table: `/Users/skku_aws2_18/team4_project/pre_project/thyroid_pipline/data/source_staging/lincs_gse70138/lincs_drug_signature_yapc_phase2_20260422.parquet`

## Input QC

- Base no-LINCS shape: `[6219, 2021]`
- Output shape: `[6219, 2534]`
- Row-level YAPC LINCS coverage: `0.3595`
- YAPC LINCS drugs: `84`

## Variant Comparison

| variant | cv | ensemble | spearman | pearson | rmse | mae | r2 | ndcg_at_20 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| pan_cancer_lincs | random4 | top4_spearman_weighted | 0.9911 | 0.9902 | 0.3926 | 0.2014 | 0.9799 | 0.9699 |
| paad_focused_lincs | random4 | top4_spearman_weighted | 0.9910 | 0.9899 | 0.3991 | 0.2047 | 0.9792 | 0.9695 |
| yapc_cell_lincs | random4 | top4_spearman_weighted | 0.9910 | 0.9904 | 0.3878 | 0.2010 | 0.9804 | 0.9839 |
| pan_cancer_lincs | groupcv4_drug | top4_spearman_weighted | 0.7512 | 0.7432 | 1.8676 | 1.3923 | 0.5453 | 0.8471 |
| paad_focused_lincs | groupcv4_drug | top4_spearman_weighted | 0.7374 | 0.7199 | 1.9292 | 1.4338 | 0.5148 | 0.8455 |
| yapc_cell_lincs | groupcv4_drug | top4_spearman_weighted | 0.7732 | 0.7736 | 1.7740 | 1.3437 | 0.5897 | 0.8523 |

## Random sample 4-fold ensemble

| ensemble | spearman | pearson | rmse | mae | r2 | ndcg_at_20 |
| --- | --- | --- | --- | --- | --- | --- |
| top4_spearman_weighted | 0.9910 | 0.9904 | 0.3878 | 0.2010 | 0.9804 | 0.9839 |
| spearman_weighted | 0.9888 | 0.9885 | 0.4530 | 0.2942 | 0.9732 | 0.9797 |
| equal_weight | 0.9887 | 0.9885 | 0.4544 | 0.2957 | 0.9731 | 0.9798 |

## Drug GroupCV 4-fold ensemble

| ensemble | spearman | pearson | rmse | mae | r2 | ndcg_at_20 |
| --- | --- | --- | --- | --- | --- | --- |
| top4_spearman_weighted | 0.7732 | 0.7736 | 1.7740 | 1.3437 | 0.5897 | 0.8523 |
| spearman_weighted | 0.7479 | 0.7639 | 1.7924 | 1.3490 | 0.5811 | 0.8556 |
| equal_weight | 0.7408 | 0.7605 | 1.8019 | 1.3563 | 0.5767 | 0.8548 |

## Random sample member models

| member | spearman | pearson | rmse | mae | r2 | ndcg_at_20 |
| --- | --- | --- | --- | --- | --- | --- |
| DL_CrossAttention | 0.9922 | 0.9916 | 0.3692 | 0.2119 | 0.9822 | 0.9839 |
| ML_LightGBM | 0.9897 | 0.9892 | 0.4091 | 0.2333 | 0.9782 | 0.9868 |
| ML_ExtraTrees | 0.9867 | 0.9852 | 0.4820 | 0.2431 | 0.9697 | 0.9529 |
| ML_RandomForest | 0.9840 | 0.9844 | 0.5001 | 0.2918 | 0.9674 | 0.9693 |
| DL_WideDeep | 0.9811 | 0.9821 | 0.5221 | 0.3316 | 0.9645 | 0.9796 |
| ML_XGBoost | 0.9761 | 0.9817 | 0.5681 | 0.4155 | 0.9579 | 0.9895 |
| DL_ResidualMLP | 0.9737 | 0.9743 | 0.6373 | 0.4321 | 0.9470 | 0.9303 |
| ML_LightGBM_DART | 0.9691 | 0.9751 | 1.1567 | 0.9475 | 0.8256 | 0.9695 |
| ML_FlatMLP | 0.9505 | 0.9571 | 0.8031 | 0.5781 | 0.9159 | 0.9208 |

## Drug GroupCV member models

| member | spearman | pearson | rmse | mae | r2 | ndcg_at_20 |
| --- | --- | --- | --- | --- | --- | --- |
| ML_XGBoost | 0.7518 | 0.7488 | 1.8368 | 1.3717 | 0.5601 | 0.8463 |
| ML_LightGBM | 0.7482 | 0.7464 | 1.8440 | 1.3692 | 0.5567 | 0.8479 |
| ML_LightGBM_DART | 0.7453 | 0.7427 | 2.0244 | 1.6375 | 0.4657 | 0.8377 |
| DL_CrossAttention | 0.7437 | 0.7630 | 1.7969 | 1.3268 | 0.5791 | 0.9279 |
| ML_ExtraTrees | 0.7242 | 0.7434 | 1.8529 | 1.3911 | 0.5524 | 0.8485 |
| ML_RandomForest | 0.7143 | 0.7177 | 1.9382 | 1.4595 | 0.5102 | 0.8608 |
| DL_ResidualMLP | 0.6026 | 0.6810 | 2.1196 | 1.6033 | 0.4143 | 0.8469 |
| ML_FlatMLP | 0.5644 | 0.6377 | 2.2174 | 1.7361 | 0.3590 | 0.8478 |
| DL_WideDeep | 0.5628 | 0.6454 | 2.2005 | 1.6876 | 0.3687 | 0.8945 |

## GroupCV ensemble Top 15

| rank | drug_name | ensemble_score | mean_pred_ln_ic50 | screened_rows | target_genes | PATHWAY_NAME_NORMALIZED | classification |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1.0000 | Docetaxel | 3.0223 | -3.0223 | 24.0000 | BCL2;NR1I2;TUBB1 | Mitosis | paad_standard_or_known |
| 2.0000 | MG-132 | 2.2226 | -2.2226 | 29.0000 | CAPN1;PROTEASOME | Protein stability and degradation | screened_candidate |
| 3.0000 | Dactinomycin | 2.1889 | -2.1889 | 24.0000 | TOP1;TOP2A | Other | screened_candidate |
| 4.0000 | Bortezomib | 0.8777 | -0.8777 | 29.0000 | PROTEASOME;PRSS1;PSMB1;PSMB5 | Protein stability and degradation | paad_trial_or_indication |
| 5.0000 | Camptothecin | 0.4422 | -0.4422 | 29.0000 | TOP1 | DNA replication | paad_trial_or_indication |
| 6.0000 | Temsirolimus | 0.2923 | -0.2923 | 26.0000 | MTOR | PI3K/MTOR signaling | paad_trial_or_indication |
| 7.0000 | Lestaurtinib | 0.2663 | -0.2663 | 27.0000 | FLT3;JAK2;NTRK1;NTRK2;NTRK3 | Other kinases | paad_trial_or_indication |
| 8.0000 | Romidepsin | 0.2120 | -0.2120 | 24.0000 | HDAC1;HDAC2;HDAC3;HDAC4;HDAC6;HDAC8 | Chromatin histone acetylation | paad_trial_or_indication |
| 9.0000 | BI-2536 | 0.1523 | -0.1523 | 27.0000 | PLK1;PLK2;PLK3 | Cell cycle | screened_candidate |
| 10.0000 | Mitoxantrone | 0.1354 | -0.1354 | 24.0000 | TOP2;TOP2A;TOP2B | DNA replication | paad_trial_or_indication |
| 11.0000 | Buparlisib | -0.3531 | 0.3531 | 29.0000 | PI3KALPHA;PI3KBETA;PI3KDELTA;PI3KGAMMA;PIK3CA;PIK3CB;PIK3CD;PIK3CG | PI3K/MTOR signaling | paad_trial_or_indication |
| 12.0000 | Pevonedistat | -0.4463 | 0.4463 | 29.0000 | NAE;NAE1;NEDD8;UBA3 | Other | screened_candidate |
| 13.0000 | Paclitaxel | -0.6071 | 0.6071 | 29.0000 | BCL2;NR1I2;TUBB1 | Mitosis | paad_standard_or_known |
| 14.0000 | Vinblastine | -0.6115 | 0.6115 | 28.0000 | JUN;TUBA1A;TUBB;TUBB2A;TUBD1;TUBE1;TUBG1 | Mitosis | paad_trial_or_indication |
| 15.0000 | SN-38 | -0.6173 | 0.6173 | 29.0000 | TOP1 | DNA replication | paad_trial_or_indication |

## Key outputs

- `data/source_staging/lincs_gse70138/lincs_drug_signature_yapc_phase2_20260422.parquet`
- `data/paad/X_numeric_strong_context_smiles_yapc_lincs.npy`
- `data/paad/numeric_strong_context_smiles_yapc_lincs_feature_names.json`
- `reports/paad/qc_paad_yapc_lincs_signature_20260422.json`
- `reports/paad/qc_paad_yapc_lincs_input_20260422.json`
- `results/paad/ensemble/paad_yapc_lincs_variant_comparison.csv`
- `results/paad/ensemble/groupcv4_drug/numeric_strong_context_smiles_yapc_lincs/ensemble_metrics.csv`
