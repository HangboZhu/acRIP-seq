#!/usr/bin/env bash
#SBATCH --job-name=staridxMouse
#SBATCH --partition=c128m1T,vmcore128,corexd192,core56
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/logs/star_index_mouse.%j.out
#SBATCH --error=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/logs/star_index_mouse.%j.err
#SBATCH --exclude=computg4
#
# Build the GRCm39 (mouse) STAR index with the container's own STAR 2.7.11a
# so the index matches the aligner at run time (same reason as the hg38 one,
# see build_star_index.sbatch). Sources: shared cellranger-arc GRCm39-2024-A
# bundle (GENCODE M33 exon-level GTF, decompressed to
# References/mouse_GRCm39/genes.gtf) + its genome.fa. sjdbOverhang 149 fits
# our 150 bp mouse reads. ~40 min, ~28G index.
set -euo pipefail

PIPE=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline
SING=/hpcdisk1/biosoft/app1/miniforge3/envs/singularity/bin/singularity
SANDBOX="${PIPE}/Containers/MeRIPseq_hg38_final.sandbox"
INDEX_DIR="${PIPE}/References/mouse_GRCm39/STAR_GRCm39_star2.7.11a"
FASTA=/hpcdisk1/XDB_vc/liushy/reference/refdata-cellranger-arc-GRCm39-2024-A/fasta/genome.fa
GTF="${PIPE}/References/mouse_GRCm39/genes.gtf"

for f in "${SANDBOX}" "${FASTA}" "${GTF}"; do
    [[ -e "${f}" ]] || { echo "ERROR: missing ${f}" >&2; exit 1; }
done

rm -rf "${INDEX_DIR}"
mkdir -p "${INDEX_DIR}" "${PIPE}/logs"

echo "=== building mouse STAR index on $(hostname) ==="
"${SING}" exec --cleanenv --bind /hpcdisk1:/hpcdisk1 \
    "${SANDBOX}" \
    STAR --runMode genomeGenerate \
         --runThreadN "${SLURM_CPUS_PER_TASK:-8}" \
         --genomeDir "${INDEX_DIR}" \
         --genomeFastaFiles "${FASTA}" \
         --sjdbGTFfile "${GTF}" \
         --sjdbOverhang 149 \
         --genomeSAsparseD 1 \
         --genomeSAindexNbases 14 \
         --genomeChrBinNbits 18 \
         2>&1 | tee "${PIPE}/logs/star_index_mouse.${SLURM_JOB_ID:-local}.log"

grep -m1 '^versionGenome' "${INDEX_DIR}/genomeParameters.txt"
echo "MOUSE_INDEX_OK"
