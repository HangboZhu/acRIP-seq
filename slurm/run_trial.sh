#!/usr/bin/env bash
# Trial run using only complete IP/Input pairs listed in config.trial.yaml.
# Data live outside the pipeline directory, so bind the data root too.
# This script lives in slurm/; the pipeline root is its parent directory.
set -euo pipefail

SLURM_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPE_DIR="$(cd "${SLURM_DIR}/.." && pwd)"

export ACRIP_CONFIG="${PIPE_DIR}/Scripts/config.trial.yaml"
export CORES="${CORES:-8}"
export PATH="/hpcdisk1/apod/shijiachen/work/conda/acrip_smk8_venv/bin:/hpcdisk1/biosoft/app1/miniforge3/envs/singularity/bin:${PATH}"
export EXTRA_BIND="${EXTRA_BIND:-/hpcdisk1:/hpcdisk1}"

exec "${SLURM_DIR}/run_acRIP.sh" "$@"
