# Reproduce every figure and result table reported in the thesis from a clean
# checkout. Designed to run in the Dockerfile environment but works on any
# host with Python 3.10+ and the requirements.txt deps installed.
#
# Usage:
#   make reproduce          # everything
#   make data               # download the phenopacket-store corpus only
#   make attacks            # MIA + k-anonymity ablation
#   make benchmark          # Phenopacket Store retrieval benchmark
#   make figures            # alias for `attacks benchmark`
#   make clean              # remove generated results and figures

PY ?= python
SEED ?= 0
TOP_N_DISEASES ?= 100
MAX_PER_DISEASE ?= 15

.PHONY: reproduce data attacks benchmark figures clean

reproduce: data attacks benchmark
	@echo "[reproduce] All figures and result tables generated."
	@echo "[reproduce]   figures/  : fig10–fig12"
	@echo "[reproduce]   results/  : privacy_attacks/, phenopacket_store/"

data:
	bash scripts/download_phenopacket_store.sh

attacks:
	$(PY) experiments/evaluate_privacy_attacks.py --seed $(SEED)

benchmark: data
	$(PY) experiments/evaluate_phenopacket_store.py \
	    --top-n-diseases $(TOP_N_DISEASES) \
	    --max-per-disease $(MAX_PER_DISEASE) \
	    --seed $(SEED)

figures: attacks benchmark

clean:
	rm -rf figures/fig1[0-2]_*.png figures/fig1[0-2]_*.pdf
	rm -rf results/privacy_attacks results/phenopacket_store
