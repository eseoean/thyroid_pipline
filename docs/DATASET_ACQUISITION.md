# Dataset Acquisition Guide

이 파이프라인은 산출물을 S3에 업로드하지 않습니다. 모든 중간 산출물과 최종 결과는 우선 로컬 레포 하위에 저장합니다.

## 1. 데이터 인벤토리

현재 로컬에 필요한 입력 파일이 있는지 확인합니다.

```bash
make inventory
```

S3 prefix도 함께 확인하려면 다음처럼 실행합니다.

```bash
python3 scripts/00_dataset_inventory.py --config config/thyroid_pipeline_config.json --s3-prefix s3://say2-4team/<thyroid-prefix>/
```

결과는 다음 파일에 저장됩니다.

```text
reports/dataset_inventory.csv
reports/dataset_inventory.json
```

## 2. 원천 데이터 내려받기 또는 복사

`config/data_manifest.example.json`을 `config/data_manifest.json`으로 복사한 뒤 `local_path` 또는 `s3_uri`를 채웁니다.

```bash
cp config/data_manifest.example.json config/data_manifest.json
python3 scripts/01_acquire_datasets.py --config config/thyroid_pipeline_config.json --manifest config/data_manifest.json --dry-run
```

dry-run 결과가 맞으면 `--dry-run`을 제거합니다.

```bash
python3 scripts/01_acquire_datasets.py --config config/thyroid_pipeline_config.json --manifest config/data_manifest.json
```

## 3. 필수 원천 데이터

- `data/raw/thyroid_response_pairs.csv`
- `data/raw/sample_features.csv`
- `data/raw/drug_features.csv`
- `data/raw/drug_annotations.csv`

## 4. 선택 검증 데이터

- `data/external/thyroid_expression.csv`
- `data/external/thyroid_clinical.csv`
- `data/admet/*.csv`

외부검증과 ADMET 데이터가 없으면 파이프라인은 해당 단계의 QC에 `no_external_expression` 또는 낮은 ADMET coverage를 기록하고 계속 진행합니다.

