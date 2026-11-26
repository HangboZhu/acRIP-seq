#!/bin/bash
#SBATCH --job-name=zhacRIP-probe
#SBATCH --partition=c128m1T,vmcore128,corexd192
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-00:30:00
#SBATCH --output=/hpcdisk1/apod/zhuhb/project/ac4c_yuliao/acRIPseq_pipeline/logs/probe_ep2.%j.out
#SBATCH --error=/hpcdisk1/apod/zhuhb/project/ac4c_yuliao/acRIPseq_pipeline/logs/probe_ep2.%j.err
#SBATCH --exclude=computg4
set -uo pipefail

export PATH="/usr/bin:/bin:/usr/sbin:/sbin:${PATH:-}"
PIPE=/hpcdisk1/apod/zhuhb/project/ac4c_yuliao/acRIPseq_pipeline
WORK="$PIPE/.zhuhb_run/dryrun"
mkdir -p "$WORK/home" "$WORK/tmp"
export HOME="$WORK/home"
export TMPDIR="$WORK/tmp"
export PATH="/hpcdisk1/biosoft/app1/miniforge3/envs/singularity/bin:$PATH"

OUT="$PIPE/logs/exomepeak2_defaults.txt"
: > "$OUT"

singularity exec \
    --cleanenv --home "$WORK/home" --env TMPDIR="$WORK/tmp" \
    --bind "$PIPE:$PIPE,/hpcdisk1:/hpcdisk1" \
    "$PIPE/Containers/MeRIPseq_hg38_final.sandbox" \
    Rscript -e 'suppressPackageStartupMessages(library(exomePeak2)); cat(as.character(packageVersion("exomePeak2")), "\n"); f <- formals(exomePeak2); for (n in names(f)) cat(n, "=", paste(deparse(f[[n]]), collapse=" "), "\n")' \
    >> "$OUT" 2>&1
echo "exit=$?" >> "$OUT"
