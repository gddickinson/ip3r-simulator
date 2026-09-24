#!/bin/bash
# Launch the IP3R Structural Simulator GUI.
# Double-click in Finder, or run ./run_app.command from a terminal.
# Extra arguments are passed through (e.g. ./run_app.command checks).
cd "$(dirname "$0")" || exit 1

CONDA_BASE="$(conda info --base 2>/dev/null)"
for base in "$CONDA_BASE" /opt/anaconda3 /opt/miniconda3 "$HOME/anaconda3" "$HOME/miniconda3" "$HOME/miniforge3"; do
    if [ -n "$base" ] && [ -f "$base/etc/profile.d/conda.sh" ]; then
        source "$base/etc/profile.d/conda.sh"
        break
    fi
done

if ! conda activate ip3r_sim 2>/dev/null; then
    echo "Could not activate the 'ip3r_sim' conda environment."
    echo "Create it with: make env"
    exit 1
fi

exec python -m ip3r "$@"
