.PHONY: inventory acquire seed-demo run-demo run clean

inventory:
	python3 scripts/00_dataset_inventory.py --config config/thyroid_pipeline_config.json

acquire:
	python3 scripts/01_acquire_datasets.py --config config/thyroid_pipeline_config.json --manifest config/data_manifest.example.json --dry-run

seed-demo:
	python3 scripts/seed_demo_data.py --config config/thyroid_pipeline_config.json

run-demo: seed-demo
	python3 scripts/run_all.py --config config/thyroid_pipeline_config.json --demo

run:
	python3 scripts/run_all.py --config config/thyroid_pipeline_config.json

clean:
	find data results reports external_validation admet knowledge_validation phase5_final_results -type f ! -name ".gitkeep" -delete
