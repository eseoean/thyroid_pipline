.PHONY: inventory acquire build-model-ready audit-missingness seed-demo run-demo run clean

inventory:
	python3 scripts/00_dataset_inventory.py --config config/thyroid_pipeline_config.json

acquire:
	python3 scripts/01_acquire_datasets.py --config config/thyroid_pipeline_config.json --manifest config/data_manifest.example.json --dry-run

build-model-ready:
	python3 scripts/02_build_model_ready_from_thyroid_raw.py --config config/thyroid_pipeline_config.json

audit-missingness:
	python3 scripts/03_audit_missingness_sources.py --config config/thyroid_pipeline_config.json

seed-demo:
	python3 scripts/seed_demo_data.py --config config/thyroid_pipeline_config.json

run-demo: seed-demo
	python3 scripts/run_all.py --config config/thyroid_pipeline_config.json --demo

run:
	python3 scripts/run_all.py --config config/thyroid_pipeline_config.json

clean:
	find data results reports external_validation admet knowledge_validation phase5_final_results -type f ! -name ".gitkeep" -delete
