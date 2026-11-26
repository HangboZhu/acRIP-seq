#!/bin/bash
#SBATCH --job-name=zhtest
#SBATCH --partition=c128m1T,vmcore128,corexd192,core56
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --exclude=computg4
#SBATCH --output=/hpcdisk1/apod/zhuhb/project/ac4c_yuliao/acRIPseq_pipeline/logs/zhtest.%j.out
#SBATCH --error=/hpcdisk1/apod/zhuhb/project/ac4c_yuliao/acRIPseq_pipeline/logs/zhtest.%j.err
set -euo pipefail
PIPE=/hpcdisk1/apod/zhuhb/project/ac4c_yuliao/acRIPseq_pipeline
export HOME="$PIPE/.zhuhb_run/home/test"
export TMPDIR="$PIPE/.zhuhb_run/tmp"
mkdir -p "$HOME" "$TMPDIR"
export PATH="/hpcdisk1/biosoft/app1/miniforge3/envs/singularity/bin:$PATH"
load_sif="$PIPE/Containers/MeRIPseq_hg38_final.sandbox"
echo "host=$(hostname) uid=$(id -u) HOME=$HOME"
singularity exec --cleanenv --bind "$PIPE:$PIPE,/hpcdisk1:/hpcdisk1" "$load_sif" \
    bash -c 'echo "inside HOME=$HOME TMPDIR=$TMPDIR"; command -v STAR; command -v trim_galore; samtools --version | head -1; Rscript -e "library(exomePeak2); cat(\"exomePeak2 OK\n\")"'
echo TEST_OK
