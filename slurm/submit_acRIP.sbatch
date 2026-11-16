#!/usr/bin/env bash
#SBATCH --job-name=acRIPseq
#SBATCH --partition=c128m1T,vmcore128,corexd192,core56
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --output=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/logs/acRIPseq.%j.out
#SBATCH --error=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/logs/acRIPseq.%j.err
#
# Example Slurm launcher. Adjust --cpus-per-task / --mem / --time to fit
# the current free resources (see AGENTS.md), then from the pipeline root:
#   sbatch slurm/submit_acRIP.sbatch
# The Snakemake --cores value follows SLURM_CPUS_PER_TASK.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export CORES="${SLURM_CPUS_PER_TASK:-20}"
"${SCRIPT_DIR}/run_acRIP.sh"
