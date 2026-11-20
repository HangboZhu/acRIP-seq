#!/usr/bin/env bash
#SBATCH --job-name=acRIPmouse
#SBATCH --partition=c128m1T,vmcore128,corexd192,core56
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/logs/acRIPmouse.%j.out
#SBATCH --error=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/logs/acRIPmouse.%j.err
#SBATCH --exclude=computg4
#
# Mouse full run: 40 peak-calling groups / 72 complete IP/Input pairs
# (144 samples). Same resource shape as the human full run (two 8-thread
# STARs in parallel under 128G). GC correction disabled -- container has no
# mouse BSgenome. Submit: sbatch slurm/submit_mouse.sbatch
set -euo pipefail

PIPE=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline
export CORES="${SLURM_CPUS_PER_TASK:-16}"
export EXECUTOR=slurm
export SLURM_JOBS="${SLURM_JOBS:-60}"
"${PIPE}/slurm/run_mouse.sh" --rerun-triggers mtime "$@"
