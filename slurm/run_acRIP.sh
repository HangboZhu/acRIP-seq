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

if [[ "${EXECUTOR:-local}" == "slurm" ]]; then
    # Cluster execution: every rule job becomes its own sbatch, spread over
    # nodes. core56 (computg*) is excluded from the partition list because
    # computg4 lacks a passwd entry for uid 1175 (singularity aborts there).
    # --jobs caps concurrent slurm jobs; --cores is ignored by the slurm
    # executor but kept harmless. --rerun-triggers mtime can be passed via
    # "$@" to resume without re-running jobs whose rule code changed.
    snakemake \
        -s Scripts/pipeline.smk \
        --use-singularity \
        --executor slurm \
        --jobs "${SLURM_JOBS:-40}" \
        --default-resources \
            "slurm_partition=c128m1T,vmcore128,corexd192" \
            "mem_mb=16000" \
            "tmpdir=/hpcdisk1/apod/shijiachen/work/tmp" \
        --printshellcmds \
        --singularity-args "--cleanenv ${BIND_ARGS}" \
        "$@"
else
snakemake \
    -s Scripts/pipeline.smk \
    --use-singularity \
    --cores "${CORES}" \
    --printshellcmds \
    --singularity-args "--cleanenv ${BIND_ARGS}" \
    "$@"
fi
