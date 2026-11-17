#!/usr/bin/env bash
#SBATCH --job-name=acRIPtrial
#SBATCH --partition=c128m1T,vmcore128,corexd192,core56
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
# Restrict to node families where uid 1175 resolves: on computg4 (core56)
# singularity aborts with "user: unknown userid 1175" because the node's
# local user db has no shijiachen entry (job 1866963 died on this before
# running any rule). nodef* (c128m1T) and node25/27 (vmcore128) are proven
# good. Remove this line once the admins fix the computg* user lookups.
#SBATCH --nodelist=nodef1,nodef2,nodef3,nodef4,node25,node27
#SBATCH --output=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/logs/acRIPtrial.%j.out
#SBATCH --error=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/logs/acRIPtrial.%j.err
#
# Trial run: 5 complete human IP/Input pairs (2 groups), see
# Scripts/config.trial.yaml. Adjust resources to the current free hole.
# Submit from the pipeline root with: sbatch slurm/submit_trial.sbatch
set -euo pipefail

# Hardcoded absolute path: slurm runs the batch script from a spool copy
# under /var/spool/slurmd, so BASH_SOURCE-based path derivation breaks
# (job 1865823 died on this). Keep the style of test_sandbox.sbatch.
PIPE=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline
export CORES="${SLURM_CPUS_PER_TASK:-8}"
"${PIPE}/slurm/run_trial.sh"
