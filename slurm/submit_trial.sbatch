#!/usr/bin/env bash
#SBATCH --job-name=acRIPtrial
#SBATCH --partition=c128m1T,vmcore128,corexd192,core56
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/logs/acRIPtrial.%j.out
#SBATCH --error=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/logs/acRIPtrial.%j.err
#
# Trial run: 5 complete human IP/Input pairs (2 groups), see
# Scripts/config.trial.yaml. Adjust resources to the current free hole.
# Submit from the pipeline root with: sbatch slurm/submit_trial.sbatch
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export CORES="${SLURM_CPUS_PER_TASK:-8}"
"${SCRIPT_DIR}/run_trial.sh"
