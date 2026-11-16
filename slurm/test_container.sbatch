#!/usr/bin/env bash
#SBATCH --job-name=acRIPcontainerTest
#SBATCH --partition=core56,corexd192,c128m1T,vmcore128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/logs/container_test.%j.out
#SBATCH --error=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/logs/container_test.%j.err
set -euo pipefail
SIF=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/Containers/MeRIPseq_hg38_final.sif
SING=/hpcdisk1/biosoft/app1/miniforge3/envs/singularity/bin/singularity
"$SING" exec --cleanenv --bind /hpcdisk1:/hpcdisk1 "$SIF" bash -lc '
  echo HOST=$(hostname)
  echo ---
  for t in STAR samtools fastqc trim_galore Rscript python multiqc; do
    printf "%s: " "$t"; command -v "$t" || true
  done
  echo ---
  STAR --version 2>&1 | head -2
  samtools --version 2>&1 | head -2
  trim_galore --version 2>&1 | head -2
  Rscript -e "suppressPackageStartupMessages(library(exomePeak2)); cat(\"exomePeak2\", as.character(packageVersion(\"exomePeak2\")), \"\\n\")" 2>&1 | tail -3
'
