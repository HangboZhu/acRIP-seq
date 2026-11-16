#!/usr/bin/env bash
#SBATCH --job-name=acRIPsandboxTest
#SBATCH --partition=core56,corexd192,c128m1T,vmcore128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/logs/sandbox_test.%j.out
#SBATCH --error=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/logs/sandbox_test.%j.err
set -euo pipefail
PIPE=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline
SING=/hpcdisk1/biosoft/app1/miniforge3/envs/singularity/bin/singularity
"$SING" exec --cleanenv --bind /hpcdisk1:/hpcdisk1 "$PIPE/Containers/MeRIPseq_hg38_final.sandbox" bash -lc '
  echo HOST=$(hostname)
  STAR --version 2>&1 | head -1
  samtools --version 2>&1 | head -1
  Rscript -e "suppressPackageStartupMessages(library(exomePeak2)); cat(\"exomePeak2\", as.character(packageVersion(\"exomePeak2\")), \"\\n\")" 2>&1 | tail -2
  echo SANDBOX_OK
'
