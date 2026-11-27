#!/bin/bash
#SBATCH --job-name=zhacRIP-fly
#SBATCH --partition=c128m1T,vmcore128,corexd192,core56
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=7-00:00:00
#SBATCH --output=/hpcdisk1/apod/zhuhb/project/ac4c_yuliao/acRIPseq_pipeline/logs/zhacRIP-fly.%j.out
#SBATCH --error=/hpcdisk1/apod/zhuhb/project/ac4c_yuliao/acRIPseq_pipeline/logs/zhacRIP-fly.%j.err
#SBATCH --exclude=computg4
set -euo pipefail

# sbatch inherits the restricted login-shell PATH (home over quota): restore coreutils.
export PATH="/usr/bin:/bin:/usr/sbin:/sbin:${PATH:-}"

PIPE=/hpcdisk1/apod/zhuhb/project/ac4c_yuliao/acRIPseq_pipeline
SP=fly
WORK="$PIPE/.zhuhb_run/$SP"

mkdir -p "$WORK/home" "$WORK/tmp" "$PIPE/logs"

# /home/zhuhb is over quota -> point HOME/TMPDIR at shared storage
export HOME="$WORK/home"
export TMPDIR="$WORK/tmp"
export ACRIP_CONFIG="$PIPE/Scripts/run_fixed/config.$SP.yaml"
export PATH="/hpcdisk1/apod/shijiachen/work/conda/acrip_smk8_venv/bin:/hpcdisk1/biosoft/app1/miniforge3/envs/singularity/bin:$PATH"

cd "$WORK"
exec snakemake \
    -s "$PIPE/Scripts/pipeline.smk" \
    --directory "$WORK" \
    --use-singularity \
    --executor slurm \
    --jobs 8 \
    --default-resources \
        "slurm_partition=c128m1T,vmcore128,corexd192" \
        "mem_mb=16000" \
        "runtime=180" \
        "tmpdir='$WORK/tmp'" \
    --printshellcmds \
    --singularity-args "--cleanenv --home $WORK/home --env TMPDIR=$WORK/tmp --bind $PIPE:$PIPE,/hpcdisk1:/hpcdisk1" \
    --nolock \
    --rerun-triggers mtime \
    "$@"
