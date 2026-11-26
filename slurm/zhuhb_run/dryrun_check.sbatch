#!/bin/bash
#SBATCH --job-name=zhacRIP-dryrun
#SBATCH --partition=c128m1T,vmcore128,corexd192
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-01:00:00
#SBATCH --output=/hpcdisk1/apod/zhuhb/project/ac4c_yuliao/acRIPseq_pipeline/logs/zhacRIP-dryrun.%j.out
#SBATCH --error=/hpcdisk1/apod/zhuhb/project/ac4c_yuliao/acRIPseq_pipeline/logs/zhacRIP-dryrun.%j.err
#SBATCH --exclude=computg4
set -uo pipefail

# sbatch inherits the restricted login-shell PATH (home over quota): restore coreutils.
export PATH="/usr/bin:/bin:/usr/sbin:/sbin:${PATH:-}"

PIPE=/hpcdisk1/apod/zhuhb/project/ac4c_yuliao/acRIPseq_pipeline
WORK="$PIPE/.zhuhb_run/dryrun"

mkdir -p "$WORK/home" "$WORK/tmp" "$PIPE/logs"
export HOME="$WORK/home"
export TMPDIR="$WORK/tmp"
export PATH="/hpcdisk1/apod/shijiachen/work/conda/acrip_smk8_venv/bin:$PATH"

OUT="$PIPE/logs/dryrun_summary.txt"
: > "$OUT"

for SP in rat fly full mouse; do
    echo "================ $SP ================" >> "$OUT"
    ACRIP_CONFIG="$PIPE/Scripts/zhuhb_run/config.$SP.yaml" \
    snakemake -s "$PIPE/Scripts/pipeline.smk" \
        --directory "$WORK" \
        --nolock \
        -n -p --cores 1 \
        --rerun-triggers mtime >> "$OUT" 2>&1
    echo "exit=$?" >> "$OUT"
done

# DAG for the smallest species (rat) and for human
ACRIP_CONFIG="$PIPE/Scripts/zhuhb_run/config.rat.yaml" \
snakemake -s "$PIPE/Scripts/pipeline.smk" --directory "$WORK" --nolock \
    --dag --rerun-triggers mtime > "$PIPE/logs/dag_rat.dot" 2> "$PIPE/logs/dag_rat.err"
echo "dag_rat exit=$?" >> "$OUT"
