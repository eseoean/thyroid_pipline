# Thyroid Raw Source Inventory 2026-04-20

Bucket prefix: `s3://say2-4team/thyroid_raw/`

This inventory records the source data collected for the thyroid cancer hybrid drug-repurposing pipeline. The rule used here is:

- Reuse existing `say2-4team` sources when they already cover the pipeline requirement.
- Download only thyroid-specific missing sources directly from public endpoints.
- Upload source data to `thyroid_raw/`; do not upload derived model outputs in this step.

## S3 Upload Verification

- Verified prefix: `s3://say2-4team/thyroid_raw/`
- Verified object count: `130`
- Verified total size: `1,937,035,410` bytes, about `1.94 GB`
- Verification command: `aws s3 ls s3://say2-4team/thyroid_raw/ --recursive --summarize`

## Thyroid Screened Response Coverage

The screened-response layer uses GDSC2 because it contains thyroid cancer cell lines with measured drug-response labels.

- THCA cell lines in GDSC2 annotation: `16`
- THCA screened response rows: `4,037`
- Unique screened drugs: `295`
- THCA cell lines with labels: `16`
- Cell lines: `K5`, `FTC-133`, `RO82-W-1`, `TT2609-C02`, `ML-1`, `TT`, `ASH-3`, `HTC-C3`, `IHH-4`, `KMH-2`, `CAL-62`, `BHT-101`, `B-CPAP`, `8505C`, `8305C`, `CGTH-W-1`

This keeps the recommendation pipeline anchored on drugs with actual screened response values, while KG/API and clinical evidence remain post-model validation layers.

## Source From say2-4team

| Source | Original S3 source | thyroid_raw destination | Main role |
|---|---|---|---|
| GDSC2 screened drug response and annotations | `s3://say2-4team/pipeline_bundle/BRCA/stage1/basic_preprocessing_20260406/gdsc/` | `s3://say2-4team/thyroid_raw/source_from_say2/gdsc/` | thyroid screened response labels, cell-line annotations, drug/pathway annotations |
| DepMap/CCLE model and CRISPR/repurposing features | `s3://say2-4team/pipeline_bundle/BRCA/stage1/basic_preprocessing_20260406/depmap/` | `s3://say2-4team/thyroid_raw/source_from_say2/depmap/` | sample features, model bridge, optional repurposing matrix reference |
| DrugBank processed source | `s3://say2-4team/pipeline_bundle/BRCA/stage1/basic_preprocessing_20260406/drugbank/` | `s3://say2-4team/thyroid_raw/source_from_say2/drugbank/` | drug metadata, synonyms, groups, targets |
| TDC ADMET processed assays | `s3://say2-4team/20260408_new_pre_project_biso/20260408_pre_project_biso_myprotocol/data/admet/` | `s3://say2-4team/thyroid_raw/source_from_say2/admet/` | ADMET safety nearest-neighbor evaluation |
| LINCS metadata/signature summaries | `s3://say2-4team/20260408_new_pre_project_biso/20260408_pre_project_biso_myprotocol/data/lincs/` | `s3://say2-4team/thyroid_raw/source_from_say2/lincs/` | drug perturbation metadata and signature features |
| LINCS normalized drug signature | `s3://say2-4team/20260408_new_pre_project_biso/20260408_pre_project_biso_myprotocol/data/lincs_drug_signature_normalized.parquet` | `s3://say2-4team/thyroid_raw/source_from_say2/lincs/lincs_drug_signature_normalized.parquet` | drug-level LINCS numeric features |
| OpenTargets disease/target associations | `s3://say2-4team/20260408_new_pre_project_biso/20260408_pre_project_biso_myprotocol/data/opentargets/` | `s3://say2-4team/thyroid_raw/source_from_say2/opentargets/` | target-thyroid relevance and KG evidence |
| ChEMBL compound and mechanism subset | `s3://say2-4team/20260408_new_pre_project_biso/20260408_pre_project_biso_myprotocol/data/chembl/` selected compound/mechanism/target tables | `s3://say2-4team/thyroid_raw/source_from_say2/chembl/` | SMILES, compound metadata, drug-target mechanism evidence |

## External Downloads

| Source | URL | thyroid_raw destination | Main role |
|---|---|---|---|
| TCGA-THCA STAR TPM expression matrix | `https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-THCA.star_tpm.tsv.gz` | `s3://say2-4team/thyroid_raw/external_downloads/tcga_thca/TCGA-THCA.star_tpm.tsv.gz` | external target expression validation |
| TCGA-THCA clinical table | `https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-THCA.clinical.tsv.gz` | `s3://say2-4team/thyroid_raw/external_downloads/tcga_thca/TCGA-THCA.clinical.tsv.gz` | external clinical covariates |
| TCGA-THCA survival table | `https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-THCA.survival.tsv.gz` | `s3://say2-4team/thyroid_raw/external_downloads/tcga_thca/TCGA-THCA.survival.tsv.gz` | survival/recurrence validation endpoint |
| GENCODE v36 annotation GTF | `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz` | `s3://say2-4team/thyroid_raw/external_downloads/gencode/gencode.v36.annotation.gtf.gz` | Ensembl gene ID to gene symbol mapping for TCGA matrix |
| ClinicalTrials.gov thyroid cancer drug studies API dump | `https://clinicaltrials.gov/api/v2/studies?query.cond=Thyroid%20Cancer&query.intr=drug&pageSize=1000&format=json` | `s3://say2-4team/thyroid_raw/external_downloads/clinical_trials/clinicaltrials_thyroid_cancer_drug_20260420.json` | clinical/KG evidence for thyroid drug candidates |

ClinicalTrials.gov dump check:

- Studies collected: `570`
- `nextPageToken`: absent, meaning the requested page covered all returned studies for this query.

## Pipeline Roles Covered

| Pipeline role | Covered by |
|---|---|
| `screened_response_labels` | GDSC2 label, cell-line, and drug annotation tables under `source_from_say2/gdsc/` |
| `sample_features` | DepMap model tables and CRISPR matrices under `source_from_say2/depmap/` |
| `drug_features` | GDSC drug annotation, DrugBank, ChEMBL compound master, LINCS normalized drug signature |
| `strong_context` | GDSC pathway/target annotation, DrugBank/ChEMBL targets, LINCS, OpenTargets |
| `external_validation` | TCGA-THCA expression, clinical, survival, and GENCODE mapping |
| `admet` | TDC ADMET assays under `source_from_say2/admet/` |
| `kg_api_validation` | OpenTargets, DrugBank, ChEMBL, and ClinicalTrials.gov thyroid cancer drug studies |

## Notes

- No raw TCGA-THCA individual GDC file manifest was found in `say2-4team`; Xena/GDC hub precompiled THCA matrices were downloaded directly instead.
- The next build step is to convert these source files into model-ready tables: `thyroid_response_pairs.csv`, `sample_features.csv`, `drug_features.csv`, `drug_annotations.csv`, `thyroid_expression.csv`, `thyroid_clinical.csv`, and ADMET/KG evidence tables.
- The final candidate recommendation should continue to use screened-response drugs as the primary ranking set. Unscreened repurposing candidates can be kept as a separate exploratory layer if needed, but should not be mixed into the main screened-response recommendation score without clear labeling.
