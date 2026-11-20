#!/usr/bin/env bash
#SBATCH --job-name=acRIPsp
#SBATCH --partition=c128m1T,vmcore128,corexd192,core56
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/logs/acRIP_%x.%j.out
#SBATCH --error=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/logs/acRIP_%x.%j.err
#SBATCH --exclude=computg4
#
# Generic per-species run. Driver env vars: SPECIES (lowercase label).
# GC correction is disabled for all non-human species (container only
# bundles the hg38 BSgenome).
set -euo pipefail

PIPE=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline

export ACRIP_CONFIG="${PIPE}/Scripts/config.${SPECIES}.yaml"
export CORES="${CORES:-8}"
export PATH="/hpcdisk1/apod/shijiachen/work/conda/acrip_smk8_venv/bin:/hpcdisk1/biosoft/app1/miniforge3/envs/singularity/bin:${PATH}"
export EXTRA_BIND="${EXTRA_BIND:-/hpcdisk1:/hpcdisk1}"

exec "${PIPE}/slurm/run_acRIP.sh" "$@"
