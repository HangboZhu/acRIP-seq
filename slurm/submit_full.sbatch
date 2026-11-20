#!/usr/bin/env bash
#SBATCH --job-name=acRIPfull
#SBATCH --partition=c128m1T,vmcore128,corexd192,core56
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/logs/acRIPfull.%j.out
#SBATCH --error=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/logs/acRIPfull.%j.err
# NOTE: this directive must sit in the header block ABOVE the descriptive
# comments below -- slurm stops parsing #SBATCH lines after the first plain
# comment, so a directive placed after them is silently ignored
# (ReqNodeList came out (null) on job 1867864). computg4 lacks a passwd
# entry for uid 1175 (singularity aborts there); computg9 works, so only
# computg4 is excluded for now.
#SBATCH --exclude=computg4
#
# Full human run: 14 peak-calling groups / 19 complete IP/Input pairs /
# 38 samples (~269 GB raw). 16 cores -> two 8-thread STARs in parallel,
# 128G covers 2 x ~35G genome indexes plus overhead. 48h wall time for the
# big MOLM13/HeLa libraries. Node restriction: see --exclude in the header.
# Submit from anywhere with: sbatch slurm/submit_full.sbatch
set -euo pipefail

PIPE=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline
export CORES="${SLURM_CPUS_PER_TASK:-16}"
# Cluster execution: each rule job = its own sbatch spread over nodes.
export EXECUTOR=slurm
export SLURM_JOBS="${SLURM_JOBS:-40}"
"${PIPE}/slurm/run_full.sh" --rerun-triggers mtime "$@"
