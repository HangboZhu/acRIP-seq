#!/usr/bin/env bash
# Driver: for each non-human species, wait for its reference download to
# finish, prepare (gunzip, GFF3->GTF for rat, read-length probe), generate
# the pairs/config, then submit index build chained with the species run.
# Run on submit01 with: nohup bash species_pipeline_driver.sh &
set -uo pipefail

PIPE=/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline
DL=$PIPE/References/downloads
MAN=/hpcdisk1/apod/shijiachen/yuliao/4/data/metadata/download_manifest.full.tsv.bak
LOG=$PIPE/logs/species_driver.log
mkdir -p "$PIPE/logs" "$PIPE/Scripts/species"

# label|manifest_species|data_dir|fasta.gz|anno.gz|anno_fmt
SPECIES_TABLE="
rat|rat|rat|GRCr8_GCA_036323735.1_genomic.fna.gz|GRCr8_GCF_036323735.1_genomic.gtf.gz|gtf
zebrafish|zebrafish|zebrafish|Danio_rerio.GRCz11.dna.primary_assembly.fa.gz|Danio_rerio.GRCz11.113.gtf.gz|gtf
fly|Drosophila|Drosophila|Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.gz|Drosophila_melanogaster.BDGP6.46.113.gtf.gz|gtf
worm|C. elegans|C_elegans|Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz|Caenorhabditis_elegans.WBcel235.113.gtf.gz|gtf
rice|rice|rice|Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz|Oryza_sativa.IRGSP-1.0.59.gtf.gz|gtf
arabidopsis|arabidopsis|arabidopsis|Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz|Arabidopsis_thaliana.TAIR10.59.gtf.gz|gtf
"

log() { echo "[$(date '+%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

readlen() {  # readlen <data_dir> -> prints most common read length
    local f
    f=$(find "/hpcdisk1/apod/shijiachen/yuliao/4/data/$1/acRIP-seq/bulk" -name "*_1.fastq.gz" 2>/dev/null | head -1)
    [ -z "$f" ] && { echo 100; return; }
    zcat "$f" | head -4004 | awk 'NR%4==2{print length($0)}' | sort | uniq -c | sort -rn | head -1 | awk '{print $2}'
}

prepare_annotation() {  # prepare_annotation <label> <anno.gz> <fmt> -> prints plain path
    local label=$1 anno=$2 fmt=$3 plain=$PIPE/References/$label/annotation.gtf
    mkdir -p "$PIPE/References/$label"
    if [ "$fmt" = gff3 ]; then
        [ -s "$plain" ] && { echo "$plain"; return; }
        /usr/bin/python3 - "$anno" "$plain" <<'PYEOF'
import sys, gzip
src, dst = sys.argv[1], sys.argv[2]
tx2gene = {}
rows = []
with gzip.open(src, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9:
            continue
        attr = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
        if c[2] == "transcript" and "ID" in attr and "Parent" in attr:
            tx2gene[attr["ID"]] = attr["Parent"].replace("gene:", "")
        elif c[2] == "exon" and "Parent" in attr:
            rows.append((c[0], c[3], c[4], c[6], attr["Parent"].replace("transcript:", "")))
with open(dst, "w") as out:
    for chrom, s, e, strand, tx in rows:
        g = tx2gene.get(tx, tx)
        out.write(f'{chrom}\t.\texon\t{s}\t{e}\t.\t{strand}\t.\tgene_id "{g}"; transcript_id "{tx}";\n')
print(f"converted {len(rows)} exons", file=sys.stderr)
PYEOF
        echo "$plain"
    else
        plain=${plain%.gtf}.gtf
        [ -s "$plain" ] || gunzip -c "$anno" > "$plain"
        echo "$plain"
    fi
}

log "driver started"

for round in $(seq 1 240); do
    all_done=1
    while IFS='|' read -r label msp datadir fasta_gz anno_gz fmt; do
        [ -z "$label" ] && continue
        [ -f "$PIPE/Scripts/species/${label}.done" ] && continue
        all_done=0
        if [ ! -s "$DL/$label/$fasta_gz" ] || [ ! -s "$DL/$label/$anno_gz" ]; then
            continue  # 还没下完
        fi
        log "=== ${label}: 参考齐了，开始准备 ==="
        mkdir -p "$PIPE/References/$label"
        FASTA=$PIPE/References/$label/genome.fa
        [ -s "$FASTA" ] || gunzip -c "$DL/$label/$fasta_gz" > "$FASTA"
        GTF=$(prepare_annotation "$label" "$DL/$label/$anno_gz" "$fmt")
        RL=$(readlen "$datadir")
        OH=$((RL - 1))
        log "${label}: fasta=$(du -h "$FASTA" | cut -f1) anno=$(du -h "$GTF" | cut -f1) readlen=${RL}"
        /usr/bin/python3 "$PIPE/Scripts/scripts/make_pairs_tsv.py" --species "$msp" --manifest "$MAN" \
            --out "$PIPE/Scripts/species/${label}_pairs.tsv" >> "$LOG" 2>&1
        if ! /usr/bin/python3 "$PIPE/Scripts/scripts/make_config_from_pairs.py" \
            --pairs-tsv "$PIPE/Scripts/species/${label}_pairs.tsv" \
            --output "$PIPE/Scripts/config.${label}.yaml" --outputdir "./result_${label}" \
            --referencedir "$PIPE/References/${label}/STAR_index_2.7.11a" \
            --gtf "$GTF" --genome-fasta "$FASTA" \
            --sif Containers/MeRIPseq_hg38_final.sandbox --threads 8 --genome "" >> "$LOG" 2>&1; then
            log "${label}: config 生成失败，跳过本轮（不标记 done），下轮重试"
            continue
        fi
        # 复用仍在排队/运行的索引作业，避免 rm -rf 掉正在写的索引
        IDX=""
        if [ -f "$PIPE/Scripts/species/${label}.idxjob" ]; then
            OLD=$(cat "$PIPE/Scripts/species/${label}.idxjob")
            if squeue -h -j "$OLD" >/dev/null 2>&1 && [ -n "$(squeue -h -j "$OLD" 2>/dev/null)" ]; then
                IDX=$OLD
                log "${label}: 复用在途索引作业 ${IDX}"
            fi
        fi
        if [ -z "$IDX" ]; then
            IDX=$(sbatch --parsable --export=ALL,SPECIES=${label},FASTA=${FASTA},GTF=${GTF},OVERHANG=${OH} \
                  "$PIPE/slurm/build_star_index_species.sbatch")
            echo "$IDX" > "$PIPE/Scripts/species/${label}.idxjob"
        fi
        RUN=$(sbatch --parsable --dependency=afterok:${IDX} --export=ALL,SPECIES=${label} \
              "$PIPE/slurm/submit_species.sbatch")
        log "${label}: index job ${IDX} -> run job ${RUN} 已提交"
        touch "$PIPE/Scripts/species/${label}.done"
    done <<< "$SPECIES_TABLE"
    [ "$all_done" = 1 ] && { log "全部物种已提交，driver 退出"; exit 0; }
    sleep 300
done
log "driver 超时退出（12小时）"
