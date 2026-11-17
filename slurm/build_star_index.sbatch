#!/usr/bin/env bash
#SBATCH --job-name=staridx
#SBATCH --partition=c128m1T,vmcore128,corexd192,core56
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/logs/star_index.%j.out
#SBATCH --error=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/logs/star_index.%j.err
#
# Rebuild the hg38 STAR genome index with the container's own STAR (2.7.11a).
# Why: the shared index /hpcdisk1/XDB_vc/huxh/ncRNA_edit/rna_editing_protocol/
# STAR/STAR_genome_index_hg38_v33 was generated with STAR 2.7.0d and is
# rejected by STAR >= 2.7.11a at alignment time:
#   "Genome version: 2.7.0d is INCOMPATIBLE with running STAR version: 2.7.11a"
#   (exit 105, killed the 2026-11-15 trial run at the star_alignment step).
# References are the same FASTA/GTF the old index was built from;
# sjdbOverhang = 125 matches our 126 bp reads; genomeSAsparseD 1 /
# genomeSAindexNbases 14 / genomeChrBinNbits 18 mirror the old build so the
# index stays ~28G and fits in 64G RAM. Runtime ~2-4 h.
#
# NOTE: all paths are hardcoded absolute paths (same style as
# test_sandbox.sbatch / extract_sandbox.sbatch). Do NOT derive paths from
# ${BASH_SOURCE[0]} here: slurm executes the batch script from a spool copy
# under /var/spool/slurmd, so BASH_SOURCE points at the spool dir and all
# derived paths break (that is exactly how job 1865823 failed).
# Submit from the pipeline root with: sbatch slurm/build_star_index.sbatch
set -euo pipefail

PIPE=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline
SING=/hpcdisk1/biosoft/app1/miniforge3/envs/singularity/bin/singularity
SANDBOX="${PIPE}/Containers/MeRIPseq_hg38_final.sandbox"
INDEX_DIR="${PIPE}/References/STAR_hg38_v33_star2.7.11a"
FASTA=/hpcdisk1/XDB_vc/huxh/ncRNA_edit/rna_editing_protocol/ref_hg38/genome/GRCh38.primary_assembly.genome.fa
GTF=/hpcdisk1/XDB_vc/huxh/ncRNA_edit/rna_editing_protocol/ref_hg38/annotation/gencode.v33.primary_assembly.annotation.gtf

for f in "${SANDBOX}" "${FASTA}" "${GTF}"; do
    [[ -e "${f}" ]] || { echo "ERROR: missing ${f}" >&2; exit 1; }
done

# A half-written index from an aborted attempt would poison a rerun
rm -rf "${INDEX_DIR}"
mkdir -p "${INDEX_DIR}" "${PIPE}/logs"

echo "=== building STAR index on $(hostname), this takes ~2-4 h ==="
"${SING}" exec --cleanenv --bind /hpcdisk1:/hpcdisk1 \
    "${SANDBOX}" \
    STAR --runMode genomeGenerate \
         --runThreadN "${SLURM_CPUS_PER_TASK:-8}" \
         --genomeDir "${INDEX_DIR}" \
         --genomeFastaFiles "${FASTA}" \
         --sjdbGTFfile "${GTF}" \
         --sjdbOverhang 125 \
         --genomeSAsparseD 1 \
         --genomeSAindexNbases 14 \
         --genomeChrBinNbits 18 \
         2>&1 | tee "${PIPE}/logs/star_index.${SLURM_JOB_ID:-local}.log"

echo "=== verifying index version (must be 2.7.11a) ==="
grep -m1 '^versionGenome' "${INDEX_DIR}/genomeParameters.txt"
echo "INDEX_BUILD_OK"
