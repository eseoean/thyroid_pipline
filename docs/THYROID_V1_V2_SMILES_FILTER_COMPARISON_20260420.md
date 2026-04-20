# Thyroid v1 vs v2 SMILES Filter Comparison 2026-04-20

이 문서는 갑상선암 파이프라인에서 SMILES missing drug를 유지했던 v1과, primary pool에서 제거한 v2를 비교한다.

## 변경 요약

| 항목 | v1 all screened drugs | v2 SMILES-required |
|---|---:|---:|
| Response rows | `4,037` | `3,387` |
| Cell lines | `16` | `16` |
| Drugs | `295` | `243` |
| SMILES missing drugs | `52` retained | `52` removed |
| Removed response rows | `0` | `650` |
| Sample-axis loss | `0` | `0` |

## Random Sample 3-Fold OOF

Primary 최종 입력셋인 `Numeric + Strong Context + SMILES` 기준 비교다.

| Version | Best model | Spearman | RMSE | R2 | NDCG@20 |
|---|---|---:|---:|---:|---:|
| v1 all screened drugs | LightGBM | `0.8370` | `1.2689` | `0.7947` | `0.9697` |
| v2 SMILES-required | LightGBM | `0.8831` | `1.1226` | `0.8491` | `0.9595` |

해석:

- SMILES missing drug를 제거하자 Spearman은 `+0.0461`, RMSE는 `-0.1463`, R2는 `+0.0544` 개선됐다.
- NDCG@20은 약간 낮아졌지만, 전체 rank correlation과 회귀 오차가 더 안정적이다.
- Cell line 손실 없이 구조 표현/ADMET 가능한 drug pool로 정리됐다는 점이 가장 큰 장점이다.

## GroupCV Stress Test

GroupCV는 random sample보다 훨씬 보수적인 stress-test다.

| Version | Input set | Model | Spearman | RMSE | R2 | NDCG@20 |
|---|---|---|---:|---:|---:|---:|
| v1 all screened drugs | Numeric + Strong Context + SMILES | ExtraTrees | `0.4481` | `2.5048` | `0.2000` | `0.8018` |
| v2 SMILES-required | Numeric + Strong Context + SMILES | ExtraTrees | `0.4515` | `2.7152` | `0.1172` | `0.6264` |

해석:

- GroupCV Spearman은 거의 유지됐지만 RMSE/R2는 악화됐다.
- 즉, SMILES filter는 noisy structure input 문제를 줄였지만, unseen cell-line/drug generalization 문제 자체를 해결하지는 못했다.
- 갑상선암은 `16`개 cell line만 사용 가능하므로, GroupCV 병목은 sample 수와 feature coverage 부족에서 온다.

## 후보 변화

v2 앙상블 top10은 다음과 같다.

| Rank | Drug | Score | 최종 검증상태 |
|---:|---|---:|---|
| 1 | Bortezomib | `4.1600` | Excluded, ADMET Caution |
| 2 | Dactinomycin | `4.1477` | Excluded, ADMET Caution |
| 3 | Docetaxel | `3.3544` | Excluded, ADMET Caution |
| 4 | Romidepsin | `3.2340` | Tier 3 |
| 5 | Sepantronium bromide | `3.0155` | Excluded, ADMET Caution |
| 6 | Staurosporine | `2.7263` | Tier 2 |
| 7 | SN-38 | `2.6638` | Excluded, ADMET Caution |
| 8 | Vinblastine | `2.5744` | Tier 2 |
| 9 | Camptothecin | `2.1505` | Excluded, ADMET Caution |
| 10 | Vinorelbine | `2.0955` | Tier 1 |

최종 추천은 모델 score 순위가 아니라 external validation, ADMET, KG/clinical validation을 통과한 Tier 기준으로 판단한다.

## 현재 결론

v2를 primary 기준으로 채택하는 것이 타당하다.

- SMILES missing drug를 제거하면 random3 성능과 downstream ADMET 일관성이 개선된다.
- 제거해도 cell line 수가 줄지 않는다.
- LINCS와 CRISPR는 hard filter하면 데이터가 너무 작아지므로, 다음 단계에서는 source 보강과 availability-aware modeling으로 결측률을 줄인다.
