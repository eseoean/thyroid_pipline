# Thyroid v2 vs v3 Missingness Enhancement Comparison 2026-04-20

이 문서는 SMILES-required v2와 LINCS/cell-line 결측 보강을 적용한 v3를 비교한다.

## 한 줄 결론

v3는 random sample 3-fold 최고 점수를 거의 유지하면서, 더 보수적인 GroupCV stress 성능을 `0.4515`에서 `0.5112`로 개선했다. 즉 최고점 경쟁보다 결측 구조 완화와 unseen-drug 안정성 개선에 의미가 있다.

## 입력셋 변화

| 항목 | v2 SMILES-required | v3 missingness-enhanced |
|---|---:|---:|
| Response rows | `3,387` | `3,387` |
| Cell lines | `16` | `16` |
| SMILES-valid drugs | `243` | `243` |
| Sample feature table | `16 x 4,114` | `16 x 5,360` |
| Drug feature table | `243 x 1,560` | `243 x 1,563` |
| Numeric input dimension | `5,332` | `6,559` |
| Numeric + Strong Context + SMILES dimension | `5,428` | `6,655` |
| LINCS direct signature drugs | `101` | `101` |
| LINCS recovered signature drugs | `0` | `14` |
| LINCS total signature drugs | `101 / 243` | `115 / 243` |
| CRISPR feature cell lines | `9 / 16` | `9 / 16` |
| CRISPR-missing cell lines with fallback | `0 / 7` | `7 / 7` |

## 성능 변화

Primary input인 `Numeric + Strong Context + SMILES` 기준 비교다.

| 기준 | Model | Spearman | RMSE | R2 | NDCG@20 |
|---|---|---:|---:|---:|---:|
| v2 random3 | LightGBM | `0.8831` | `1.1226` | `0.8491` | `0.9595` |
| v3 random3 | LightGBM | `0.8835` | `1.1264` | `0.8481` | `0.9689` |
| v2 GroupCV | ExtraTrees | `0.4515` | `2.7152` | `0.1172` | `0.6264` |
| v3 GroupCV | ExtraTrees | `0.5112` | `2.6134` | `0.1822` | `0.7222` |

## 후보 변화

| 구분 | v2 | v3 |
|---|---|---|
| Raw ensemble top1 | Bortezomib | Bortezomib |
| Raw ensemble top5 | Bortezomib, Dactinomycin, Docetaxel, Romidepsin, Sepantronium bromide | Bortezomib, Dactinomycin, Docetaxel, Sepantronium bromide, Vinblastine |
| Final Tier 1 | Vinorelbine | Vinorelbine |
| Final Tier 2 | Staurosporine, Vinblastine | Staurosporine, Vinblastine |
| Final Tier 3 | Romidepsin | Romidepsin |

최종 Tier 후보는 크게 흔들리지 않았다. 이는 후단의 ADMET/KG/clinical validation이 raw model score를 그대로 따르지 않고, 안전성과 근거 축을 다시 적용하기 때문이다.

## 구현 변경

변경 파일:

- `scripts/02_build_model_ready_from_thyroid_raw.py`
- `scripts/03_audit_missingness_sources.py`
- `config/thyroid_pipeline_config.json`

핵심 구현:

- `lincs_pert_info_basic_20260406.parquet`에서 normalized name 또는 canonical SMILES로 bridge 후보를 찾는다.
- `lincs_mcf7.parquet`에서 bridge 후보 BRD/perturbagen signature row를 찾아 drug 단위 평균 signature로 복구한다.
- 복구 signature는 `drug__lincs__*` feature로 들어가되, `drug__lincs_signature_source_recovered` flag를 함께 둔다.
- DepMap 24Q2 expression/CNV/mutation wide table을 ModelID 기준으로 thyroid cell line에 붙인다.
- CRISPR 없는 cell line에는 `sample_has_non_crispr_omics_fallback` flag를 둔다.

## 해석

- SMILES 없는 약물은 여전히 primary pool에서 제거하는 것이 맞다. 현재 local source만으로 제거된 52개 drug의 신뢰 가능한 SMILES 복구 후보는 없었다.
- LINCS는 hard filter가 아니라 soft indicator와 recovered feature가 적절하다. hard filter를 걸면 후보 pool이 너무 작아진다.
- CRISPR는 갑상선암에서 cell line 수가 너무 적으므로 hard filter를 걸면 안 된다. fallback omics feature로 sample biology signal을 보강하는 쪽이 안전하다.
- v3의 가장 중요한 성과는 random3 최고점이 아니라 GroupCV stress 개선이다.

## 재현

```bash
make build-model-ready
python3 scripts/03_audit_missingness_sources.py --config config/thyroid_pipeline_config.json
python3 scripts/run_all.py --config config/thyroid_pipeline_config.json
```

주요 확인 파일:

- `reports/qc_source_to_model_ready_20260420.json`
- `reports/missingness/missingness_source_audit_summary.json`
- `results/numeric_strong_context_smiles_metrics_summary.csv`
- `results/groupcv_stress_test/numeric_strong_context_smiles_ExtraTrees_groupcv.json`
- `phase5_final_results/final_comprehensive_candidates.csv`
