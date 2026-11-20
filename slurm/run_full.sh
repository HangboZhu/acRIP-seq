#!/usr/bin/env bash
# Full human run: every complete IP/Input pair from
# downloaded_complete_pairs.tsv (14 groups / 19 pairs / 38 samples), see
# Scripts/config.full.yaml. Data live outside the pipeline directory, so the
# data root is bound as well. Absolute paths only (slurm spool lesson).
set -euo pipefail

PIPE=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline

export ACRIP_CONFIG="${PIPE}/Scripts/config.full.yaml"
export CORES="${CORES:-16}"
export PATH="/hpcdisk1/apod/shijiachen/work/conda/acrip_smk8_venv/bin:/hpcdisk1/biosoft/app1/miniforge3/envs/singularity/bin:${PATH}"
export EXTRA_BIND="${EXTRA_BIND:-/hpcdisk1:/hpcdisk1}"

exec "${PIPE}/slurm/run_acRIP.sh" "$@"
