#!/bin/bash
#SBATCH --job-name=zhacRIP-dryrun
#SBATCH --partition=c128m1T,vmcore128,corexd192,core56
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

# Each species gets its own scratch directory: separate .snakemake metadata,
# no cross-talk between configs.
for SP in rat fly full mouse; do
    SPWORK="$WORK/$SP"
    rm -rf "$SPWORK"
    mkdir -p "$SPWORK"

    echo "================ $SP ================" >> "$OUT"
    ACRIP_CONFIG="$PIPE/Scripts/zhuhb_run/config.$SP.yaml" \
    snakemake -s "$PIPE/Scripts/pipeline.smk" \
        --directory "$SPWORK" \
        --nolock \
        -n -p --cores 1 \
        --rerun-triggers mtime >> "$OUT" 2>&1
    echo "exit=$?" >> "$OUT"
done

# DAG for the smallest species (rat) and the largest (human)
for SP in rat full; do
    ACRIP_CONFIG="$PIPE/Scripts/zhuhb_run/config.$SP.yaml" \
    snakemake -s "$PIPE/Scripts/pipeline.smk" --directory "$WORK/$SP" --nolock \
        --dag --rerun-triggers mtime > "$PIPE/logs/dag_$SP.dot" 2> "$PIPE/logs/dag_$SP.err"
    echo "dag_$SP exit=$?" >> "$OUT"
done

# Hard gate: the consensus rules must be present in the planned jobs.
echo "================ checks ================" >> "$OUT"
echo "peak_consensus_job_lines=$(grep -c '^peak_consensus  ' "$OUT")" >> "$OUT"
echo "consensus_summary_job_lines=$(grep -c '^consensus_summary  ' "$OUT")" >> "$OUT"
echo "n_exit_0=$(grep -c '^exit=0$' "$OUT")" >> "$OUT"
echo "n_exit_nonzero=$(grep '^exit=' "$OUT" | grep -vc '^exit=0$')" >> "$OUT"
