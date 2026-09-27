# IP3R structural simulator. Every target runs in the `ip3r_sim` conda env.
ENV_NAME ?= ip3r_sim
PY := conda run --no-capture-output -n $(ENV_NAME) python

.DEFAULT_GOAL := help
.PHONY: help env fetch sync extended sync-check params test test-quick checks states gui screenshots screenshots-full screenshot-groups sizes lint

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-12s\033[0m %s\n", $$1, $$2}'

env:  ## Create the conda environment
	bash scripts/create_env.sh $(ENV_NAME)

fetch:  ## Download every registry structure into ref/structures
	$(PY) -m ip3r fetch

sync:  ## Re-import curated resources from ../ip3r_genes
	$(PY) scripts/sync_genes.py

extended:  ## Re-curate the extended IP3R deposits (every other full-length EM tetramer) from RCSB
	$(PY) scripts/curate_ip3r.py

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

STEPS ?= all

screenshots:  ## GUI smoke test, no findings checks; STEPS=lumen,extras runs only those groups (list: make screenshot-groups)
	$(PY) scripts/screenshot_app.py --steps $(STEPS)

screenshots-full:  ## GUI smoke test with every findings check first (~13 min more); refreshes the findings screenshot
	$(PY) scripts/screenshot_app.py --checks

screenshot-groups:  ## List the smoke test's step groups
	@$(PY) scripts/screenshot_app.py --list

sizes:  ## Fail if any Python file exceeds 500 lines
	@find ip3r scripts tests -name '*.py' -exec wc -l {} + | awk '$$1 > 500 && $$2 != "total" {print "TOO LONG: " $$2; bad=1} END {exit bad}'

lint:  ## Static checks
	conda run -n $(ENV_NAME) ruff check ip3r scripts tests
