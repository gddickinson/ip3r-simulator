# IP3R structural simulator. Every target runs in the `ip3r_sim` conda env.
ENV_NAME ?= ip3r_sim
PY := conda run --no-capture-output -n $(ENV_NAME) python

.DEFAULT_GOAL := help
.PHONY: help env fetch sync sync-check params test test-quick checks states gui screenshots sizes lint

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-12s\033[0m %s\n", $$1, $$2}'

env:  ## Create the conda environment
	bash scripts/create_env.sh $(ENV_NAME)

fetch:  ## Download every registry structure into ref/structures
	$(PY) -m ip3r fetch

sync:  ## Re-import curated resources from ../ip3r_genes
	$(PY) scripts/sync_genes.py

ryr:  ## Re-curate the RyR1 resource (sequence, domains, state panel) from UniProt/InterPro/RCSB
	$(PY) scripts/curate_ryr.py

sync-check:  ## Exit 1 if any ip3r_genes source changed since the last sync
	$(PY) scripts/sync_genes.py --check

params:  ## Rebuild the parameter registry (validates every citation)
	$(PY) scripts/build_parameters.py

test:  ## Full test suite (real-data tests skip when data is absent)
	$(PY) -m pytest

test-quick:  ## Tests that need no structures and no ip3r_genes
	$(PY) -m pytest -k "not real and not calibration and not resources_in_step"

checks:  ## Re-derive the ip3r_genes findings (writes data/derived/checks.json)
	$(PY) -m ip3r checks --json data/derived/checks.json --figures data/derived/exhibits

states:  ## Measure the pore of every ITPR3 gating state
	$(PY) -m ip3r states

gui:  ## Launch the application
	$(PY) -m ip3r

screenshots:  ## GUI smoke test; writes docs/img screenshots
	$(PY) scripts/screenshot_app.py --checks

sizes:  ## Fail if any Python file exceeds 500 lines
	@find ip3r scripts tests -name '*.py' -exec wc -l {} + | awk '$$1 > 500 && $$2 != "total" {print "TOO LONG: " $$2; bad=1} END {exit bad}'

lint:  ## Static checks
	conda run -n $(ENV_NAME) ruff check ip3r scripts tests
