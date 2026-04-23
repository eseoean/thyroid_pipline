# PAAD External Validation Source Acquisition

Date: 2026-04-23

This manifest records the PAAD external validation and post-processing data sources staged locally and prepared for upload to `s3://say2-4team/PAAD_raw/`.

## Newly acquired public sources

| Axis | Source | Local path | S3 target | Verification |
|---|---|---|---|---|
| Normal pancreas baseline | GTEx v8 Pancreas sample-level TPM | `thyroid_pipline/data/paad_source_staging/gtex/gene_tpm_2017-06-05_v8_pancreas.gct.gz` | `s3://say2-4team/PAAD_raw/gtex/` | gzip OK, 56,203 lines |
| Normal tissue baseline | GTEx v8 median tissue TPM | `thyroid_pipline/data/paad_source_staging/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz` | `s3://say2-4team/PAAD_raw/gtex/` | gzip OK, 56,203 lines |
| GEO external validation | GSE62452 PAAD tumor/adjacent matrix | `thyroid_pipline/data/paad_source_staging/geo/GSE62452/GSE62452_series_matrix.txt.gz` | `s3://say2-4team/PAAD_raw/geo/GSE62452/` | gzip OK, 33,364 lines |
| GEO annotation | GPL6244 Affymetrix HuGene 1.0 ST annotation for GSE62452 probe-to-gene mapping | `thyroid_pipline/data/paad_source_staging/geo/GPL6244/GPL6244.annot.gz` | `s3://say2-4team/PAAD_raw/geo/GPL6244/` | gzip OK, uploaded 2026-04-23 |
| GEO external validation | GSE71729 PAAD tumor/stroma subtype matrix | `thyroid_pipline/data/paad_source_staging/geo/GSE71729/GSE71729_series_matrix.txt.gz` | `s3://say2-4team/PAAD_raw/geo/GSE71729/` | gzip OK, 19,820 lines |
| CPTAC validation | cBioPortal CPTAC GDC PAAD archive and extracted files | `thyroid_pipline/data/paad_source_staging/cptac_paad/` | `s3://say2-4team/PAAD_raw/cptac_paad/` | tar OK, extracted; TPM 40,797 lines, mutation 14,777 lines |
| Proteomics validation | LinkedOmics CPTAC-PDAC tumor proteomics | `thyroid_pipline/data/paad_source_staging/cptac_paad/proteomics_gene_level_MD_abundance_tumor.cct` | `s3://say2-4team/PAAD_raw/cptac_paad/` | 11,663 lines |
| Proteomics validation | LinkedOmics CPTAC-PDAC normal proteomics | `thyroid_pipline/data/paad_source_staging/cptac_paad/proteomics_gene_level_MD_abundance_normal.cct` | `s3://say2-4team/PAAD_raw/cptac_paad/` | 11,663 lines |
| Proteomics validation | LinkedOmics CPTAC-PDAC clinical table | `thyroid_pipline/data/paad_source_staging/cptac_paad/clinical_table_140.tsv` | `s3://say2-4team/PAAD_raw/cptac_paad/` | 141 lines |
| cBioPortal validation | TCGA PAAD PanCancer Atlas 2018 archive and extracted files | `thyroid_pipline/data/paad_source_staging/cbioportal/paad_tcga_pan_can_atlas_2018*` | `s3://say2-4team/PAAD_raw/cbioportal/paad_tcga_pan_can_atlas_2018/` | tar OK, extracted; expression 20,532 lines, mutation 31,395 lines, RPPA 199 lines |
| Safety filtering | SIDER drug names | `thyroid_pipline/data/paad_source_staging/sider/drug_names.tsv` | `s3://say2-4team/PAAD_raw/sider/` | 34,759 bytes |
| Safety filtering | SIDER MedDRA all side effects | `thyroid_pipline/data/paad_source_staging/sider/meddra_all_se.tsv.gz` | `s3://say2-4team/PAAD_raw/sider/` | gzip OK, 309,849 lines |
| Safety filtering | SIDER MedDRA frequency | `thyroid_pipline/data/paad_source_staging/sider/meddra_freq.tsv.gz` | `s3://say2-4team/PAAD_raw/sider/` | gzip OK, 291,632 lines |
| Network validation | STRING human protein links v12.0 | `thyroid_pipline/data/paad_source_staging/string/9606.protein.links.v12.0.txt.gz` | `s3://say2-4team/PAAD_raw/string/` | gzip OK, 13,715,405 lines |
| Network validation | STRING human protein info v12.0 | `thyroid_pipline/data/paad_source_staging/string/9606.protein.info.v12.0.txt.gz` | `s3://say2-4team/PAAD_raw/string/` | gzip OK, uploaded 2026-04-23 |
| Pathway validation | MSigDB Hallmark GMT v7.5 local reference | `thyroid_pipline/data/paad_source_staging/msigdb/h.all.v7.5.symbols.gmt` | `s3://say2-4team/PAAD_raw/msigdb/` | 48,244 bytes |

## Already present in PAAD raw source area

| Axis | Source | Existing location |
|---|---|---|
| Disease cohort | TCGA-PAAD RNA TPM, clinical, survival | `s3://say2-4team/PAAD_raw/additional_sources/tcga_paad/` |
| Drug-response / external evidence | DepMap repurposing / PRISM-like files | `s3://say2-4team/PAAD_raw/depmap/` |
| Clinical relevance | ClinicalTrials PAAD cancer drug query | `s3://say2-4team/PAAD_raw/additional_sources/clinical_trials/` |
| Toxicity / ADMET | ADMET files | `s3://say2-4team/PAAD_raw/admet/` |
| Disease-target relevance | OpenTargets files | `s3://say2-4team/PAAD_raw/opentargets/` |

## Excluded or blocked

| Source | Status | Reason |
|---|---|---|
| cBioPortal `paad_msk_2025` archive | Excluded from upload | DataHub asset URL returned a 111 byte XML AccessDenied document instead of a gzip archive. |
| GEO GPL20769 annotation | Excluded from upload | Download endpoint returned a small XML/error document instead of a valid gzip annotation; GSE71729 was parsed directly because the matrix rows are already gene symbols. |
| COSMIC | Not acquired | COSMIC requires licensed/login-gated access. It should be added manually if the team has credentials and redistribution permission. |

## Intended validation use

1. Build PAAD disease signatures using TCGA-PAAD plus GTEx Pancreas.
2. Check GEO reproducibility with GSE62452 and GSE71729.
3. Validate pathway/protein-level consistency with CPTAC-PDAC proteomics and cBioPortal CPTAC PAAD.
4. Use STRING and MSigDB for network/pathway support of target and pair features.
5. Use ADMET plus SIDER plus ClinicalTrials for post-model safety and clinical relevance filtering.
