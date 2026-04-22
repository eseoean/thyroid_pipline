# Thyroid Raw Source Inventory 2026-04-20

Bucket prefix: `s3://say2-4team/thyroid_raw/`

This inventory records the source data collected for the thyroid cancer hybrid drug-repurposing pipeline. The rule used here is:

- Reuse existing `say2-4team` sources when they already cover the pipeline requirement.
- Download only thyroid-specific missing sources directly from public endpoints.
- Upload source data to `thyroid_raw/`; do not upload derived model outputs in this step.

## S3 Upload Verification

- Verified prefix: `s3://say2-4team/thyroid_raw/`
- Verified object count: `134`
- Verified total size: `1,973,365,038` bytes, about `1.97 GB`
- Verification command: `aws s3 ls s3://say2-4team/thyroid_raw/ --recursive --summarize`

## Thyroid Screened Response Coverage

The screened-response layer uses GDSC2 because it contains thyroid cancer cell lines with measured drug-response labels.

- THCA cell lines in GDSC2 annotation: `16`
- THCA screened response rows: `4,037`
- Unique screened drugs: `295`
- Primary model-ready rows after SMILES filter: `3,387`
- Primary model-ready SMILES-valid drugs: `243`
- SMILES-missing drugs excluded from primary pool: `52`
- THCA cell lines with labels: `16`
- Cell lines: `K5`, `FTC-133`, `RO82-W-1`, `TT2609-C02`, `ML-1`, `TT`, `ASH-3`, `HTC-C3`, `IHH-4`, `KMH-2`, `CAL-62`, `BHT-101`, `B-CPAP`, `8505C`, `8305C`, `CGTH-W-1`

This keeps the recommendation pipeline anchored on drugs with actual screened response values. The primary pool now additionally requires valid SMILES because structure features and ADMET validation are core downstream requirements, while KG/API and clinical evidence remain post-model validation layers.

## Source From say2-4team

| Source | Original S3 source | thyroid_raw destination | Main role |
|---|---|---|---|
| GDSC2 screened drug response and annotations | `s3://say2-4team/pipeline_bundle/BRCA/stage1/basic_preprocessing_20260406/gdsc/` | `s3://say2-4team/thyroid_raw/GDSC/` | thyroid screened response labels, cell-line annotations, drug/pathway annotations |
| GDSC2 original response CSV | `s3://say2-4team/raw_data/GDSC2/GDSC2-dataset.csv` | `s3://say2-4team/thyroid_raw/GDSC/GDSC2-dataset.csv` | source-level provenance for screened response values |
| BRCA GDSC drug feature catalog reference | local BRCA pipeline cache `drug_features_catalog.parquet` | `s3://say2-4team/thyroid_raw/GDSC/drug_features_catalog_brca_reference_20260420.parquet` | SMILES bridge for GDSC drug IDs used by the thyroid model-ready builder |
| DepMap/CCLE model and CRISPR/repurposing features | `s3://say2-4team/pipeline_bundle/BRCA/stage1/basic_preprocessing_20260406/depmap/` | `s3://say2-4team/thyroid_raw/depmap/` | sample features, model bridge, optional repurposing matrix reference |
| DrugBank processed source | `s3://say2-4team/pipeline_bundle/BRCA/stage1/basic_preprocessing_20260406/drugbank/` | `s3://say2-4team/thyroid_raw/drugbank/` | drug metadata, synonyms, groups, targets |
| TDC ADMET processed assays | `s3://say2-4team/20260408_new_pre_project_biso/20260408_pre_project_biso_myprotocol/data/admet/` | `s3://say2-4team/thyroid_raw/admet/` | ADMET safety nearest-neighbor evaluation |
| LINCS metadata/signature summaries | `s3://say2-4team/20260408_new_pre_project_biso/20260408_pre_project_biso_myprotocol/data/lincs/` | `s3://say2-4team/thyroid_raw/L1000/` | drug perturbation metadata and signature features |
| LINCS normalized drug signature | `s3://say2-4team/20260408_new_pre_project_biso/20260408_pre_project_biso_myprotocol/data/lincs_drug_signature_normalized.parquet` | `s3://say2-4team/thyroid_raw/L1000/lincs_drug_signature_normalized.parquet` | drug-level LINCS numeric features |
| OpenTargets disease/target associations | `s3://say2-4team/20260408_new_pre_project_biso/20260408_pre_project_biso_myprotocol/data/opentargets/` | `s3://say2-4team/thyroid_raw/opentargets/` | target-thyroid relevance and KG evidence |
| ChEMBL compound and mechanism subset | `s3://say2-4team/20260408_new_pre_project_biso/20260408_pre_project_biso_myprotocol/data/chembl/` selected compound/mechanism/target tables | `s3://say2-4team/thyroid_raw/chembl/` | SMILES, compound metadata, drug-target mechanism evidence |

## Local Enhancement Sources

아래 파일은 v3 결측 보강을 위해 로컬 staging에 내려받아 사용했다. 생성 산출물은 S3에 업로드하지 않았다.

| Source | Local path | Main role |
|---|---|---|
| LINCS MCF7 raw signature parquet | `data/source_staging/lincs/lincs_mcf7.parquet` | LINCS bridge 후보의 실제 perturbation signature 복구 |
| DepMap 24Q2 expression | `data/source_staging/depmap/OmicsExpressionProteinCodingGenesTPMLogp1_24Q2.csv` | CRISPR missing cell line expression fallback |
| DepMap 24Q2 CNV | `data/source_staging/depmap/OmicsCNGene_24Q2.csv` | CRISPR missing cell line CNV fallback |
| DepMap 24Q2 mutation | `data/source_staging/depmap/OmicsSomaticMutations_24Q2.csv` | CRISPR missing cell line mutation fallback |

## External Downloads

| Source | URL | thyroid_raw destination | Main role |
|---|---|---|---|
| TCGA-THCA STAR TPM expression matrix | `https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-THCA.star_tpm.tsv.gz` | `s3://say2-4team/thyroid_raw/additional_sources/tcga_thca/TCGA-THCA.star_tpm.tsv.gz` | external target expression validation |
| TCGA-THCA clinical table | `https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-THCA.clinical.tsv.gz` | `s3://say2-4team/thyroid_raw/additional_sources/tcga_thca/TCGA-THCA.clinical.tsv.gz` | external clinical covariates |
| TCGA-THCA survival table | `https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-THCA.survival.tsv.gz` | `s3://say2-4team/thyroid_raw/additional_sources/tcga_thca/TCGA-THCA.survival.tsv.gz` | survival/recurrence validation endpoint |
| GENCODE v36 annotation GTF | `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz` | `s3://say2-4team/thyroid_raw/additional_sources/gencode/gencode.v36.annotation.gtf.gz` | Ensembl gene ID to gene symbol mapping for TCGA matrix |
| ClinicalTrials.gov thyroid cancer drug studies API dump | `https://clinicaltrials.gov/api/v2/studies?query.cond=Thyroid%20Cancer&query.intr=drug&pageSize=1000&format=json` | `s3://say2-4team/thyroid_raw/additional_sources/clinical_trials/clinicaltrials_thyroid_cancer_drug_20260420.json` | clinical/KG evidence for thyroid drug candidates |

ClinicalTrials.gov dump check:

- Studies collected: `570`
- `nextPageToken`: absent, meaning the requested page covered all returned studies for this query.

## Pipeline Roles Covered

| Pipeline role | Covered by |
|---|---|
| `screened_response_labels` | GDSC2 label, cell-line, and drug annotation tables under `GDSC/` |
| `sample_features` | DepMap model tables and CRISPR matrices under `depmap/` |
| `drug_features` | GDSC drug annotation, DrugBank, ChEMBL compound master, LINCS normalized drug signature |
| `strong_context` | GDSC pathway/target annotation, DrugBank/ChEMBL targets, LINCS, OpenTargets |
| `external_validation` | TCGA-THCA expression, clinical, survival, and GENCODE mapping |
| `admet` | TDC ADMET assays under `admet/` |
| `kg_api_validation` | OpenTargets, DrugBank, ChEMBL, and ClinicalTrials.gov thyroid cancer drug studies |

## Model-Ready Build Check

`scripts/02_build_model_ready_from_thyroid_raw.py` converts the staged `thyroid_raw/` sources into the actual pipeline input files.

```bash
python3 scripts/01_acquire_datasets.py --config config/thyroid_pipeline_config.json --manifest config/data_manifest.thyroid_raw_20260420.json
python3 scripts/02_build_model_ready_from_thyroid_raw.py --config config/thyroid_pipeline_config.json
python3 scripts/run_all.py --config config/thyroid_pipeline_config.json
```

Latest model-ready QC from the actual source build:

| Check | Value |
|---|---:|
| Response rows before SMILES filter | `4,037` |
| Response rows after SMILES filter | `3,387` |
| THCA cell lines | `16` |
| Screened drugs before SMILES filter | `295` |
| SMILES-valid screened drugs | `243` |
| Removed SMILES-missing drugs | `52` |
| Sample feature table | `16 x 4,114` |
| Sample feature table after v3 fallback | `16 x 5,360` |
| Drug feature table | `243 x 1,563` |
| SMILES present / parse OK | `243 / 243` |
| Target gene present after SMILES filter | `212 / 243` |
| LINCS direct matched drugs | `101 / 243` |
| LINCS recovered matched drugs | `14 / 243` |
| LINCS total matched drugs | `115 / 243` |
| CRISPR feature cell lines | `9 / 16` |
| CRISPR-missing cell lines with omics fallback | `7 / 7` |
| External TCGA-THCA expression genes written | `296` |
| External TCGA-THCA samples | `572` |
| ADMET assays written | `22` |

## Notes

- No raw TCGA-THCA individual GDC file manifest was found in `say2-4team`; Xena/GDC hub precompiled THCA matrices were downloaded directly instead.
- The source-to-model-ready builder converts these source files into `thyroid_response_pairs.csv`, `sample_features.csv`, `drug_features.csv`, `drug_annotations.csv`, `thyroid_expression.csv`, `thyroid_clinical.csv`, and `data/admet/tdc/*.csv`.
- The final candidate recommendation should continue to use screened-response drugs with valid SMILES as the primary ranking set. Unscreened repurposing candidates can be kept as a separate exploratory layer if needed, but should not be mixed into the main screened-response recommendation score without clear labeling.
- LINCS and CRISPR missingness are intentionally kept as soft availability features rather than hard filters because hard filtering would shrink the thyroid sample/drug axes too aggressively.
- v3 applies LINCS recovered signature flags and non-CRISPR omics fallback flags so model inputs retain whether a feature was direct, recovered, or unavailable.
