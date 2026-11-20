#!/usr/bin/env bash
#SBATCH --job-name=staridxSp
#SBATCH --partition=c128m1T,vmcore128,corexd192,core56
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/logs/star_index_%x.%j.out
#SBATCH --error=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/logs/star_index_%x.%j.err
#SBATCH --exclude=computg4
#
# Generic per-species STAR index builder. Driven by environment variables
# (set by the driver via sbatch --export):
#   SPECIES   rat|zebrafish|fly|worm|rice|arabidopsis
#   FASTA     uncompressed genome FASTA path
#   GTF       uncompressed GTF/GFF path
#   OVERHANG  sjdbOverhang (read length - 1)
set -euo pipefail

PIPE=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline
SING=/hpcdisk1/biosoft/app1/miniforge3/envs/singularity/bin/singularity
SANDBOX="${PIPE}/Containers/MeRIPseq_hg38_final.sandbox"
INDEX_DIR="${PIPE}/References/${SPECIES}/STAR_index_2.7.11a"

case "${SPECIES}" in
    rat)        NB=14 ;;
    zebrafish)  NB=13 ;;
    rice)       NB=12 ;;
    fly)        NB=12 ;;
    worm)       NB=11 ;;
    arabidopsis) NB=12 ;;
    *)          NB=14 ;;
esac

for f in "${SANDBOX}" "${FASTA}" "${GTF}"; do
    [[ -e "${f}" ]] || { echo "ERROR: missing ${f}" >&2; exit 1; }
done

rm -rf "${INDEX_DIR}"
mkdir -p "${INDEX_DIR}" "${PIPE}/logs"

echo "=== building ${SPECIES} STAR index on $(hostname), SAindexNbases=${NB} ==="
"${SING}" exec --cleanenv --bind /hpcdisk1:/hpcdisk1 "${SANDBOX}" \
    STAR --runMode genomeGenerate \
         --runThreadN "${SLURM_CPUS_PER_TASK:-8}" \
         --genomeDir "${INDEX_DIR}" \
         --genomeFastaFiles "${FASTA}" \
         --sjdbGTFfile "${GTF}" \
         --sjdbOverhang "${OVERHANG:-99}" \
         --genomeSAsparseD 1 \
         --genomeSAindexNbases "${NB}" \
         --genomeChrBinNbits 18 \
         2>&1 | tee "${PIPE}/logs/star_index_${SPECIES}.${SLURM_JOB_ID:-local}.log"

grep -m1 '^versionGenome' "${INDEX_DIR}/genomeParameters.txt"
echo "${SPECIES}_INDEX_OK"
