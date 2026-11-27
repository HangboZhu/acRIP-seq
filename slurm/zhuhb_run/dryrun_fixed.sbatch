#!/bin/bash
#SBATCH --job-name=zhdry
#SBATCH --partition=c128m1T,vmcore128,corexd192,core56
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0-02:00:00
#SBATCH --output=/hpcdisk1/apod/zhuhb/project/ac4c_yuliao/acRIPseq_pipeline/logs/zhdry.%j.out
#SBATCH --error=/hpcdisk1/apod/zhuhb/project/ac4c_yuliao/acRIPseq_pipeline/logs/zhdry.%j.err
#SBATCH --exclude=computg4
set -uo pipefail

export PATH="/usr/bin:/bin:/usr/sbin:/sbin:${PATH:-}"
PIPE=/hpcdisk1/apod/zhuhb/project/ac4c_yuliao/acRIPseq_pipeline

for SP in fly rat full mouse; do
    WORK="$PIPE/.zhuhb_run/$SP"
    mkdir -p "$WORK/home" "$WORK/tmp"
    export HOME="$WORK/home"
    export TMPDIR="$WORK/tmp"
    export ACRIP_CONFIG="$PIPE/Scripts/run_fixed/config.$SP.yaml"
    export PATH="/hpcdisk1/apod/shijiachen/work/conda/acrip_smk8_venv/bin:/hpcdisk1/biosoft/app1/miniforge3/envs/singularity/bin:$PATH"
    cd "$WORK"
    echo "########## DRYRUN $SP ##########"
    snakemake -s "$PIPE/Scripts/pipeline.smk" \
        --directory "$WORK" \
        --use-singularity \
        --executor slurm \
        --jobs 8 \
        --default-resources \
            "slurm_partition=c128m1T,vmcore128,corexd192" \
            "mem_mb=16000" \
            "runtime=180" \
        --nolock \
        --rerun-triggers mtime \
        -n --summary
    echo "########## DRYRUN $SP exit=$? ##########"
done
