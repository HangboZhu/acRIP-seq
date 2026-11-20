#!/usr/bin/env bash
# Mouse full run: every complete IP/Input pair (40 groups / 72 pairs), see
# Scripts/config.mouse.yaml. GC correction is DISABLED (genome: "") because
# the container only bundles the hg38 BSgenome.
set -euo pipefail

PIPE=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline

export ACRIP_CONFIG="${PIPE}/Scripts/config.mouse.yaml"
export CORES="${CORES:-16}"
export PATH="/hpcdisk1/apod/shijiachen/work/conda/acrip_smk8_venv/bin:/hpcdisk1/biosoft/app1/miniforge3/envs/singularity/bin:${PATH}"
export EXTRA_BIND="${EXTRA_BIND:-/hpcdisk1:/hpcdisk1}"

exec "${PIPE}/slurm/run_acRIP.sh" "$@"
