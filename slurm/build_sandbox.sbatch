#!/usr/bin/env bash
#SBATCH --job-name=acRIPsandbox
#SBATCH --partition=core56,corexd192,c128m1T,vmcore128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/logs/sandbox_build.%j.out
#SBATCH --error=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/logs/sandbox_build.%j.err
set -euo pipefail
PIPE=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline
SING=/hpcdisk1/biosoft/app1/miniforge3/envs/singularity/bin/singularity
export PATH="/usr/sbin:${PATH}"
export SINGULARITY_TMPDIR=/hpcdisk1/apod/shijiachen/work/tmp/singularity
mkdir -p "$SINGULARITY_TMPDIR"
cd "$PIPE"
rm -rf Containers/MeRIPseq_hg38_final.sandbox
"$SING" build --sandbox Containers/MeRIPseq_hg38_final.sandbox Containers/MeRIPseq_hg38_final.sif
echo BUILD_DONE
du -sh Containers/MeRIPseq_hg38_final.sandbox
