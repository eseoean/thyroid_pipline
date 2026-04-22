# Thyroid S3 Source Layout - 2026-04-22

The thyroid raw source bundle is stored under:

`s3://say2-4team/thyroid_raw/`

For S3 console readability, the source bundle is now exposed directly by
source family at the top level:

| S3 prefix | Source family |
| --- | --- |
| `thyroid_raw/GDSC/` | GDSC2 response and annotation sources |
| `thyroid_raw/depmap/` | DepMap model, CRISPR, and repurposing-screen sources |
| `thyroid_raw/L1000/` | LINCS/L1000 perturbation signature sources |
| `thyroid_raw/admet/` | ADMET assay sources |
| `thyroid_raw/chembl/` | ChEMBL compound, mechanism, and target sources |
| `thyroid_raw/drugbank/` | DrugBank drug, synonym, target, and identifier sources |
| `thyroid_raw/opentargets/` | OpenTargets disease, target, and association sources |
| `thyroid_raw/additional_sources/clinical_trials/` | Thyroid cancer clinical-trial download |
| `thyroid_raw/additional_sources/gencode/` | GENCODE annotation |
| `thyroid_raw/additional_sources/tcga_thca/` | TCGA-THCA expression, clinical, and survival files |
| `thyroid_raw/manifests/` | Raw source inventory and source manifest |

The older nested prefixes are retained for backward compatibility:

- `thyroid_raw/source_from_say2/`
- `thyroid_raw/external_downloads/`

The final thyroid result bundle remains separate:

`s3://say2-4team/20260409_eseo/20260421_thyroid/`

That result bundle contains model inputs, reports, validation outputs, final
candidates, and repo snapshots. It does not duplicate the full raw/staging
source cache.
