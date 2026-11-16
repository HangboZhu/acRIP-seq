#!/usr/bin/env bash
# Run the acRIP-seq pipeline (Snakemake + Singularity).
# This script lives in slurm/; the pipeline root is its parent directory.
# Usage: slurm/run_acRIP.sh [snakemake extra args...]
#   CORES=20 slurm/run_acRIP.sh
#   EXTRA_BIND=/path/to/data:/path/to/data slurm/run_acRIP.sh
set -euo pipefail

PIPE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CORES="${CORES:-8}"
EXTRA_BIND="${EXTRA_BIND:-}"

BIND_ARGS="--bind ${PIPE_DIR}:${PIPE_DIR}"
if [[ -n "${EXTRA_BIND}" ]]; then
    BIND_ARGS="${BIND_ARGS},${EXTRA_BIND}"
fi

cd "${PIPE_DIR}"

snakemake \
    -s Scripts/pipeline.smk \
    --use-singularity \
    --cores "${CORES}" \
    --printshellcmds \
    --singularity-args "--cleanenv ${BIND_ARGS}" \
    "$@"
