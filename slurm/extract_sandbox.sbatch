#!/usr/bin/env bash
#SBATCH --job-name=acRIPextract
#SBATCH --partition=core56,corexd192,c128m1T,vmcore128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/logs/extract_sandbox.%j.out
#SBATCH --error=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/logs/extract_sandbox.%j.err
set -euo pipefail
PIPE=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline
SING=/hpcdisk1/biosoft/app1/miniforge3/envs/singularity/bin/singularity
SQUASH=/hpcdisk1/apod/shijiachen/work/tmp/merip_rootfs.squashfs
echo HOST=$(hostname)
"$SING" sif dump 3 "$PIPE/Containers/MeRIPseq_hg38_final.sif" > "$SQUASH"
ls -lh "$SQUASH"
rm -rf "$PIPE/Containers/MeRIPseq_hg38_final.sandbox"
/usr/sbin/unsquashfs -d "$PIPE/Containers/MeRIPseq_hg38_final.sandbox" "$SQUASH"
rm -f "$SQUASH"
echo EXTRACT_DONE
du -sh "$PIPE/Containers/MeRIPseq_hg38_final.sandbox"
