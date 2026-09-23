#!/usr/bin/env bash
# Create the `ip3r_sim` conda environment.
# Usage:  bash scripts/create_env.sh   (then: conda activate ip3r_sim)
# On a machine that already has the PIEZO1 simulator's env, the fastest route
# is:  conda create -n ip3r_sim --clone piezo1   (which is how it was made).
set -euo pipefail
ENV_NAME="${1:-ip3r_sim}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda create -y -n "$ENV_NAME" -c conda-forge python=3.11 numpy scipy matplotlib pytest ruff
conda run -n "$ENV_NAME" python -m pip install "PyQt6>=6.6" moderngl
echo "==> done. activate with:  conda activate $ENV_NAME"
