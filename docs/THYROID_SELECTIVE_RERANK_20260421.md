# Thyroid-Selective Re-ranking - 2026-04-21

## 목적

pan-cancer LINCS 모델의 상위 후보가 BRCA 추천 후보와 과도하게 겹치는 문제를 줄이기 위해, 모델 예측 점수만이 아니라 thyroid cell line에서의 상대 선택성, thyroid 관련 target/pathway context, broad cytotoxic penalty를 함께 반영한 재랭킹 레이어를 추가했다.

## 왜 필요한가

- 기존 pan-cancer LINCS Top 후보는 BRCA random3 strong-context-smiles Top 후보와 많이 겹쳤다.
- 기존 겹침: Top10 9/10, Top15 14/15, Top30 25/30.
- 이는 오류라기보다 GDSC screened pool에서 전반적으로 강한 항암제와 broad cytotoxic 후보가 공통으로 상위에 올라오는 구조 때문이다.

## Score 정의

`thyroid_selective_score = 0.30 * model_rank_score + 0.30 * selectivity_rank_score + 0.20 * thyroid_response_rank_score + 0.20 * thyroid_context_score - broad_cytotoxic_penalty`

- `model_rank_score`: pan-cancer LINCS ensemble 예측 민감도 rank.
- `selectivity_rank_score`: `pancancer_mean_ln_ic50 - thyroid_mean_ln_ic50` rank. 값이 클수록 thyroid가 pan-cancer 평균보다 더 민감하다.
- `thyroid_response_rank_score`: thyroid screened response의 observed sensitivity rank.
- `thyroid_context_score`: BRAF/RET/NTRK/RAS/MAPK/VEGFR/PI3K/MTOR/CDK 등 thyroid 관련 target/pathway 점수.
- `broad_cytotoxic_penalty`: mitosis, microtubule, DNA replication/damage, proteasome 등 broad cytotoxic mechanism 및 pan-cancer broad activity에 대한 감점.

## 전체 재랭킹 Top 20

| thyroid_selective_rank | drug_name | thyroid_selective_score | ranking_layer | ensemble_score | thyroid_mean_ln_ic50 | pancancer_mean_ln_ic50 | selectivity_delta_ln_ic50 | thyroid_context_score | broad_cytotoxic_flag |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | Trametinib | 0.9603 | thyroid_context_supported | 0.4800 | -0.7924 | 0.2335 | 1.0259 | 1.0000 | 0 |
| 2 | Staurosporine | 0.9550 | thyroid_context_supported | 2.7392 | -3.2425 | -2.7856 | 0.4569 | 0.8500 | 0 |
| 3 | PD0325901 | 0.9050 | thyroid_context_supported | -0.7097 | 0.3465 | 1.3590 | 1.0125 | 0.8500 | 0 |
| 4 | Refametinib | 0.8696 | thyroid_context_supported | -1.8725 | 1.3340 | 2.2903 | 0.9562 | 0.8500 | 0 |
| 5 | Ulixertinib | 0.8409 | thyroid_context_supported | -2.1544 | 2.6246 | 2.8043 | 0.1797 | 0.9500 | 0 |
| 6 | SCH772984 | 0.8298 | thyroid_context_supported | -2.5917 | 2.1448 | 2.4595 | 0.3148 | 0.8500 | 0 |
| 7 | Selumetinib | 0.8276 | thyroid_context_supported | -2.7179 | 2.9058 | 3.9322 | 1.0264 | 0.9500 | 0 |
| 8 | Dasatinib | 0.8052 | screened_selectivity_candidate | -1.0116 | 0.5224 | 1.2985 | 0.7762 | 0.4000 | 0 |
| 9 | Dabrafenib | 0.7987 | thyroid_context_supported | -2.9742 | 3.5183 | 4.2657 | 0.7474 | 1.0000 | 0 |
| 10 | BX795 | 0.7502 | screened_selectivity_candidate | -1.4708 | 2.2385 | 2.6917 | 0.4532 | 0.3000 | 0 |
| 11 | AZD2014 | 0.7375 | screened_selectivity_candidate | -2.5965 | 2.4239 | 2.1336 | -0.2903 | 0.8500 | 0 |
| 12 | VX-11e | 0.7221 | screened_selectivity_candidate | -2.7657 | 2.9459 | 2.7060 | -0.2399 | 0.8500 | 0 |
| 13 | Cediranib | 0.7216 | screened_selectivity_candidate | -2.6095 | 2.6434 | 2.2658 | -0.3776 | 0.9500 | 0 |
| 14 | Podophyllotoxin bromide | 0.7184 | screened_selectivity_candidate | -0.4023 | -0.5751 | -0.4445 | 0.1306 | 0.0000 | 0 |
| 15 | OSI-027 | 0.6931 | screened_selectivity_candidate | -2.5971 | 2.9310 | 2.5583 | -0.3727 | 0.8500 | 0 |
| 16 | PLX-4720 | 0.6867 | thyroid_context_supported | -3.5410 | 4.0699 | 4.2122 | 0.1423 | 0.8500 | 0 |
| 17 | Lestaurtinib | 0.6632 | screened_selectivity_candidate | -0.2120 | 0.1173 | -0.0151 | -0.1324 | 0.8500 | 0 |
| 18 | AZD4547 | 0.6473 | screened_selectivity_candidate | -3.2155 | 3.1526 | 2.7184 | -0.4341 | 0.9500 | 0 |
| 19 | AZD8055 | 0.6448 | screened_selectivity_candidate | -0.2659 | 0.1678 | -0.0366 | -0.2044 | 0.8500 | 0 |
| 20 | BMS-536924 | 0.6397 | screened_selectivity_candidate | -2.0575 | 2.1550 | 2.1095 | -0.0456 | 0.0000 | 0 |

## Broad Cytotoxic 제외 Targeted/Selective Top 20

| thyroid_selective_rank | thyroid_targeted_rank | drug_name | thyroid_selective_score | ranking_layer | ensemble_score | thyroid_mean_ln_ic50 | pancancer_mean_ln_ic50 | selectivity_delta_ln_ic50 | thyroid_context_score | broad_cytotoxic_flag |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 1 | Trametinib | 0.9603 | thyroid_context_supported | 0.4800 | -0.7924 | 0.2335 | 1.0259 | 1.0000 | 0 |
| 2 | 2 | Staurosporine | 0.9550 | thyroid_context_supported | 2.7392 | -3.2425 | -2.7856 | 0.4569 | 0.8500 | 0 |
| 3 | 3 | PD0325901 | 0.9050 | thyroid_context_supported | -0.7097 | 0.3465 | 1.3590 | 1.0125 | 0.8500 | 0 |
| 4 | 4 | Refametinib | 0.8696 | thyroid_context_supported | -1.8725 | 1.3340 | 2.2903 | 0.9562 | 0.8500 | 0 |
| 5 | 5 | Ulixertinib | 0.8409 | thyroid_context_supported | -2.1544 | 2.6246 | 2.8043 | 0.1797 | 0.9500 | 0 |
| 6 | 6 | SCH772984 | 0.8298 | thyroid_context_supported | -2.5917 | 2.1448 | 2.4595 | 0.3148 | 0.8500 | 0 |
| 7 | 7 | Selumetinib | 0.8276 | thyroid_context_supported | -2.7179 | 2.9058 | 3.9322 | 1.0264 | 0.9500 | 0 |
| 8 | 8 | Dasatinib | 0.8052 | screened_selectivity_candidate | -1.0116 | 0.5224 | 1.2985 | 0.7762 | 0.4000 | 0 |
| 9 | 9 | Dabrafenib | 0.7987 | thyroid_context_supported | -2.9742 | 3.5183 | 4.2657 | 0.7474 | 1.0000 | 0 |
| 10 | 10 | BX795 | 0.7502 | screened_selectivity_candidate | -1.4708 | 2.2385 | 2.6917 | 0.4532 | 0.3000 | 0 |
| 11 | 11 | AZD2014 | 0.7375 | screened_selectivity_candidate | -2.5965 | 2.4239 | 2.1336 | -0.2903 | 0.8500 | 0 |
| 12 | 12 | VX-11e | 0.7221 | screened_selectivity_candidate | -2.7657 | 2.9459 | 2.7060 | -0.2399 | 0.8500 | 0 |
| 13 | 13 | Cediranib | 0.7216 | screened_selectivity_candidate | -2.6095 | 2.6434 | 2.2658 | -0.3776 | 0.9500 | 0 |
| 14 | 14 | Podophyllotoxin bromide | 0.7184 | screened_selectivity_candidate | -0.4023 | -0.5751 | -0.4445 | 0.1306 | 0.0000 | 0 |
| 15 | 15 | OSI-027 | 0.6931 | screened_selectivity_candidate | -2.5971 | 2.9310 | 2.5583 | -0.3727 | 0.8500 | 0 |
| 16 | 16 | PLX-4720 | 0.6867 | thyroid_context_supported | -3.5410 | 4.0699 | 4.2122 | 0.1423 | 0.8500 | 0 |
| 17 | 17 | Lestaurtinib | 0.6632 | screened_selectivity_candidate | -0.2120 | 0.1173 | -0.0151 | -0.1324 | 0.8500 | 0 |
| 18 | 18 | AZD4547 | 0.6473 | screened_selectivity_candidate | -3.2155 | 3.1526 | 2.7184 | -0.4341 | 0.9500 | 0 |
| 19 | 19 | AZD8055 | 0.6448 | screened_selectivity_candidate | -0.2659 | 0.1678 | -0.0366 | -0.2044 | 0.8500 | 0 |
| 20 | 20 | BMS-536924 | 0.6397 | screened_selectivity_candidate | -2.0575 | 2.1550 | 2.1095 | -0.0456 | 0.0000 | 0 |

## BRCA 후보와 겹침 변화

### 전체 재랭킹

- Top10 overlap: 1/10
- Top15 overlap: 1/15
- Top30 overlap: 7/30

### Broad Cytotoxic 제외 Targeted/Selective 리스트

- Top10 overlap: 1/10
- Top15 overlap: 1/15
- Top30 overlap: 6/30

## Downstream validation

- External validation rows: 15
- ADMET rows: 15
- Knowledge validation rows: 15
- Tier counts: `{"Tier 2": 10, "Excluded": 3, "Tier 1": 2}`

## 산출물

- `results/thyroid_selective_rerank/thyroid_selective_reranked_top30.csv`
- `results/thyroid_selective_rerank/thyroid_targeted_non_broad_reranked_top30.csv`
- `external_validation/thyroid_selective/top15_validated.csv`
- `admet/thyroid_selective/final_drug_candidates.csv`
- `phase5_final_results/thyroid_selective/final_comprehensive_candidates.csv`
- `reports/qc_thyroid_selective_rerank_20260421.json`
