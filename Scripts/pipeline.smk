import os
from pathlib import Path
from datetime import datetime

# Load configuration file (relative to this Snakefile).
# ACRIP_CONFIG can point to an alternative config, e.g. config.trial.yaml.
configfile: os.environ.get("ACRIP_CONFIG", os.path.join(workflow.basedir, "config.yaml"))

# ---------------------------------------------------------------------------
# Path handling: every path in config.yaml is resolved relative to the
# pipeline root (the parent directory of Scripts/) unless already absolute.
# ---------------------------------------------------------------------------
PROJECT_DIR = Path(workflow.basedir).resolve().parent


def resolve_path(value):
    if value is None or value == "":
        return None
    p = Path(str(value))
    return p if p.is_absolute() else (PROJECT_DIR / p).resolve()


OUTPUT_DIR = resolve_path(config.get("outputdir", "./result"))
REFERENCEDIR = resolve_path(config["referencedir"])
GTF = resolve_path(config["gtf"])
GENOME_FASTA = resolve_path(config.get("genome_fasta"))
SIF = str(resolve_path(config["sif"]))
SCAN_MOTIF = resolve_path(config.get("scan_motif"))

# Output sub-directories
QC_DIR = OUTPUT_DIR / "qc"
TRIM_DIR = OUTPUT_DIR / "trim"
BAM_DIR = OUTPUT_DIR / "bam"
MOTIF_DIR = OUTPUT_DIR / "motif_analysis"
SCANMOTIF_DIR = OUTPUT_DIR / "scan_motif"
MULTIQC_DIR = OUTPUT_DIR / "multiqc"
REPORT_DIR = OUTPUT_DIR / "report"
REPORT_DATA_DIR = REPORT_DIR / "data"
REPORT_FIG_DIR = REPORT_DIR / "figures"
LOG_DIR = OUTPUT_DIR / "logs"

for dir_path in [QC_DIR, TRIM_DIR, BAM_DIR, MOTIF_DIR, SCANMOTIF_DIR,
                 MULTIQC_DIR, REPORT_DIR, REPORT_DATA_DIR, REPORT_FIG_DIR, LOG_DIR]:
    dir_path.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Parse sample groups.  A group may contain multiple biological replicates;
# all replicate BAMs of a group are passed together to exomePeak2 for joint
# peak calling (OpenAc4C unified acRIP-seq pipeline).
#
# Each replicate may be single-end (1 FASTQ) or paired-end (R1 + R2).
# Accepted per-group FASTQ fields:
#   ip_fastq:      [ [R1,R2], "single.fq.gz", ... ]  one entry per sample
#   ip_fastq_r1:   [R1, R1, ...]                     one R1 per sample
#   ip_fastq_r2:   [R2, "", ...]                     optional; "" = single-end
# and the same for input_fastq / input_fastq_r1 / input_fastq_r2.
# A flat 2-element list with a single sample is also treated as R1/R2.
# ---------------------------------------------------------------------------
def as_list(value):
    if value is None or value == "":
        return []
    if isinstance(value, (list, tuple)):
        return [str(v) for v in value]
    return [str(value)]


def fastq_entry(raw):
    """Normalize one sample FASTQ entry to a list of 1 or 2 resolved paths."""
    if raw is None or raw == "":
        return []
    if isinstance(raw, (list, tuple)):
        return [resolve_path(x) for x in raw if x not in (None, "")]
    return [resolve_path(raw)]


def parse_fastq_field(raw, sample_names):
    """Return a list of per-sample FASTQ path lists."""
    n = len(sample_names)
    if raw is None or raw == "":
        return [[] for _ in range(n)]
    if n == 1:
        return [fastq_entry(raw)]
    if not isinstance(raw, (list, tuple)) or len(raw) != n:
        raise ValueError(
            f"Expected {n} FASTQ entries (one per sample), got: {raw!r}"
        )
    return [fastq_entry(x) for x in raw]


def parse_fastq_r1_r2(pair, prefix, sample_names):
    """Parse <prefix>_fastq_r1 / <prefix>_fastq_r2 lists, if present."""
    r1_raw = pair.get(f"{prefix}_fastq_r1")
    if r1_raw is None:
        return None
    r1_list = as_list(r1_raw)
    r2_list = as_list(pair.get(f"{prefix}_fastq_r2"))
    if len(r1_list) != len(sample_names):
        raise ValueError(
            f"pair '{pair['pair_name']}': {len(r1_list)} {prefix}_fastq_r1 "
            f"but {len(sample_names)} {prefix}_samples"
        )
    if r2_list and len(r2_list) != len(sample_names):
        raise ValueError(
            f"pair '{pair['pair_name']}': {len(r2_list)} {prefix}_fastq_r2 "
            f"but {len(sample_names)} {prefix}_samples"
        )
    entries = []
    for i in range(len(sample_names)):
        fqs = [r1_list[i]]
        if r2_list and r2_list[i]:
            fqs.append(r2_list[i])
        entries.append([resolve_path(f) for f in fqs])
    return entries


SAMPLE_PAIRS = config.get("sample_pairs", [])
if not SAMPLE_PAIRS:
    # Fallback to legacy single-pair config for backward compatibility
    SAMPLE_PAIRS = [{
        "pair_name": "default",
        "ip_sample": config.get("ip_sample", ""),
        "input_sample": config.get("input_sample", None),
        "ip_fastq": config.get("ip_fastq", ""),
        "input_fastq": config.get("input_fastq", "")
    }]

SAMPLE_FASTQ = {}
SAMPLE_LAYOUT = {}

for pair in SAMPLE_PAIRS:
    pair.setdefault("pair_name", "pair")

    ip_raw = pair.get("ip_samples")
    if ip_raw is None:
        ip_raw = pair.get("ip_sample")
    input_raw = pair.get("input_samples")
    if input_raw is None:
        input_raw = pair.get("input_sample")

    ip_samples = as_list(ip_raw)
    input_samples = as_list(input_raw)

    ip_entries = parse_fastq_r1_r2(pair, "ip", ip_samples)
    if ip_entries is None:
        ip_entries = parse_fastq_field(pair.get("ip_fastq"), ip_samples)
    input_entries = parse_fastq_r1_r2(pair, "input", input_samples)
    if input_entries is None:
        input_entries = parse_fastq_field(pair.get("input_fastq"), input_samples)

    for name, fqs in zip(ip_samples, ip_entries):
        if len(fqs) not in (1, 2):
            raise ValueError(
                f"sample '{name}': expected 1 or 2 FASTQ files, got {len(fqs)}"
            )
        SAMPLE_FASTQ[name] = fqs
        SAMPLE_LAYOUT[name] = "pe" if len(fqs) == 2 else "se"
    for name, fqs in zip(input_samples, input_entries):
        if len(fqs) not in (1, 2):
            raise ValueError(
                f"sample '{name}': expected 1 or 2 FASTQ files, got {len(fqs)}"
            )
        SAMPLE_FASTQ[name] = fqs
        SAMPLE_LAYOUT[name] = "pe" if len(fqs) == 2 else "se"

    pair["ip_samples"] = ip_samples
    pair["input_samples"] = input_samples
    pair["ip_fastq"] = ip_entries
    pair["input_fastq"] = input_entries
    pair["ip_sample"] = ip_samples[0] if ip_samples else ""
    pair["input_sample"] = input_samples[0] if input_samples else None

PAIR_NAMES = [p["pair_name"] for p in SAMPLE_PAIRS]
PAIR_INFO = {p["pair_name"]: p for p in SAMPLE_PAIRS}
ALL_SAMPLES = list(SAMPLE_FASTQ.keys())

# Trimmed FASTQ paths produced by the SE / PE trim rules
TRIMMED_FASTQ = {}
for _sample in ALL_SAMPLES:
    if SAMPLE_LAYOUT[_sample] == "pe":
        TRIMMED_FASTQ[_sample] = [
            TRIM_DIR / f"trim_{_sample}_pe_val_1.fq.gz",
            TRIM_DIR / f"trim_{_sample}_pe_val_2.fq.gz",
        ]
    else:
        TRIMMED_FASTQ[_sample] = [TRIM_DIR / f"trim_{_sample}_se_trimmed.fq.gz"]


def trim_log_path(sample):
    return TRIM_DIR / f"{sample}.{SAMPLE_LAYOUT[sample]}.trim.log"


def trim_log_targets():
    return [trim_log_path(s) for s in ALL_SAMPLES]


# Global settings
GENOME = config.get("genome", "hg38")
SCAN_MOTIF_PVALUE = config.get("scan_motif_pvalue", 0.01)
STAR_EXTRA_ARGS = config.get("star_extra_args", "")
MAPQ_FILTER = config.get("mapq_filter", 255)
REQUIRE_PROPER_PAIR = config.get("require_proper_pair", True)


def csv_bam_paths(wildcards, role):
    return ",".join(
        str(BAM_DIR / f"{sample}.DeDup.bam")
        for sample in PAIR_INFO[wildcards.pair][f"{role}_samples"]
    )


rule all:
    input:
        # Per-sample QC and preprocessing (single-end or paired-end)
        [QC_DIR / f"{s}.fastqc.log" for s in ALL_SAMPLES],
        [f for s in ALL_SAMPLES for f in TRIMMED_FASTQ[s]],
        # Per-sample alignment and dedup
        [BAM_DIR / f"{s}.STAR.log" for s in ALL_SAMPLES],
        [BAM_DIR / f"{s}.DeDup.log" for s in ALL_SAMPLES],
        # MultiQC report (single, across all samples)
        MULTIQC_DIR / "multiqc_report.html",
        # Per-group joint peak calling
        expand(OUTPUT_DIR / "exomePeak2_output_{pair}" / "peaks.bed", pair=PAIR_NAMES) if not config.get("skip_peak", False) else [],
        # Optional liftOver harmonization
        expand(OUTPUT_DIR / "exomePeak2_output_{pair}" / "peaks.lifted.bed", pair=PAIR_NAMES) if config.get("liftover", False) else [],
        # Merge enriched peaks across all sample groups
        OUTPUT_DIR / "acRIP_all_peaks.bed" if not config.get("skip_peak", False) else [],
        # Per-pair motif analysis (off by default for acRIP-seq)
        expand(MOTIF_DIR / "streme_{pair}.html", pair=PAIR_NAMES) if not config.get("skip_peak", False) and not config.get("skip_motif", False) else [],
        expand(MOTIF_DIR / "peak_{pair}.fa", pair=PAIR_NAMES) if not config.get("skip_peak", False) and not config.get("skip_motif", False) else [],
        # Per-pair scan motif
        expand(SCANMOTIF_DIR / "fimo_{pair}.tsv", pair=PAIR_NAMES) if not config.get("skip_motif", False) else [],
        expand(SCANMOTIF_DIR / "fimo_{pair}.html", pair=PAIR_NAMES) if not config.get("skip_motif", False) else [],
        expand(SCANMOTIF_DIR / "fimo_{pair}_motif.bed", pair=PAIR_NAMES) if not config.get("skip_motif", False) else [],
        # Unified report
        REPORT_DIR / "acRIP_Analysis_Report_CN.html",
        REPORT_DIR / "acRIP_Analysis_Report_EN.html"

rule fastqc:
    input:
        fqs = lambda wildcards: SAMPLE_FASTQ[wildcards.sample]
    output:
        log = QC_DIR / "{sample}.fastqc.log"
    log:
        LOG_DIR / "fastqc" / "{sample}.log"
    params:
        threads = config.get("threads", 8),
        outdir = QC_DIR,
        layout = lambda wildcards: SAMPLE_LAYOUT[wildcards.sample],
        r1 = lambda wildcards: SAMPLE_FASTQ[wildcards.sample][0],
        r2 = lambda wildcards: SAMPLE_FASTQ[wildcards.sample][1] if len(SAMPLE_FASTQ[wildcards.sample]) > 1 else ""
    container: SIF
    shell:
        r"""
        mkdir -p {params.outdir} {LOG_DIR}/fastqc
        TMP=$(mktemp -d)
        mkdir -p "$TMP/r1"
        if [ "{params.layout}" = "pe" ]; then
            mkdir -p "$TMP/r2"
            fastqc -o "$TMP/r2" -t {params.threads} --noextract "{params.r2}" > {output.log} 2> {log}
        fi
        fastqc -o "$TMP/r1" -t {params.threads} --noextract "{params.r1}" >> {output.log} 2>> {log}

        b1=$(basename "{params.r1}"); b1=${{b1%.gz}}; b1=${{b1%.fq}}; b1=${{b1%.fastq}}
        if [ "{params.layout}" = "pe" ]; then
            b2=$(basename "{params.r2}"); b2=${{b2%.gz}}; b2=${{b2%.fq}}; b2=${{b2%.fastq}}
            mv "$TMP/r1/${{b1}}_fastqc.html" "{params.outdir}/{wildcards.sample}_R1_fastqc.html"
            mv "$TMP/r2/${{b2}}_fastqc.html" "{params.outdir}/{wildcards.sample}_R2_fastqc.html"
        else
            mv "$TMP/r1/${{b1}}_fastqc.html" "{params.outdir}/{wildcards.sample}_fastqc.html"
        fi
        rm -rf $TMP
        """

rule trim_galore_se:
    input:
        fq = lambda wildcards: SAMPLE_FASTQ[wildcards.sample][0]
    output:
        trimmed = TRIM_DIR / "trim_{sample}_se_trimmed.fq.gz",
        log = TRIM_DIR / "{sample}.se.trim.log"
    log:
        LOG_DIR / "trim_galore" / "{sample}.se.log"
    params:
        threads = config.get("threads", 8),
        outdir = TRIM_DIR,
        qual = 25,
        length = 36,
        e = 0.1,
        stringency = 3
    container: SIF
    shell:
        """
        mkdir -p {params.outdir} {LOG_DIR}/trim_galore
        trim_galore -q {params.qual} --cores {params.threads} --phred33 --length {params.length} \
            -e {params.e} --stringency {params.stringency} --fastqc --illumina \
            --basename "trim_{wildcards.sample}_se" -o {params.outdir} "{input.fq}" \
            > {output.log} 2> {log}
        """

rule trim_galore_pe:
    input:
        r1 = lambda wildcards: SAMPLE_FASTQ[wildcards.sample][0],
        r2 = lambda wildcards: SAMPLE_FASTQ[wildcards.sample][1]
    output:
        trimmed_1 = TRIM_DIR / "trim_{sample}_pe_val_1.fq.gz",
        trimmed_2 = TRIM_DIR / "trim_{sample}_pe_val_2.fq.gz",
        log = TRIM_DIR / "{sample}.pe.trim.log"
    log:
        LOG_DIR / "trim_galore" / "{sample}.pe.log"
    params:
        threads = config.get("threads", 8),
        outdir = TRIM_DIR,
        qual = 25,
        length = 36,
        e = 0.1,
        stringency = 3
    container: SIF
    shell:
        """
        mkdir -p {params.outdir} {LOG_DIR}/trim_galore
        trim_galore --paired -q {params.qual} --cores {params.threads} --phred33 --length {params.length} \
            -e {params.e} --stringency {params.stringency} --fastqc --illumina \
            --basename "trim_{wildcards.sample}_pe" -o {params.outdir} "{input.r1}" "{input.r2}" \
            > {output.log} 2> {log}
        """

rule star_alignment:
    input:
        trimmed = lambda wildcards: TRIMMED_FASTQ[wildcards.sample]
    output:
        bam = temp(BAM_DIR / "{sample}.Aligned.sortedByCoord.out.bam"),
        log = BAM_DIR / "{sample}.STAR.log"
    log:
        LOG_DIR / "star_alignment" / "{sample}.log"
    # threads MUST be a directive (not params): snakemake sizes scheduling
    # and per-job core allocation from it. With threads hidden in params
    # every star job claims 1 core, so up to 8 STARs run in parallel and
    # each loads a ~35G genome index -> node OOM (exit 137), which is how
    # the 2026-11-16 22:45 run died.
    threads: config.get("threads", 8)
    params:
        ref_dir = REFERENCEDIR,
        outdir = BAM_DIR,
        extra_args = STAR_EXTRA_ARGS,
        read_files = lambda wildcards: " ".join(
            str(f) for f in TRIMMED_FASTQ[wildcards.sample]
        )
    container: SIF
    shell:
        """
        mkdir -p {params.outdir} {LOG_DIR}/star_alignment
        STAR --runThreadN {threads} \
             --genomeDir "{params.ref_dir}" \
             --readFilesIn {params.read_files} \
             --readFilesCommand zcat \
             --outFileNamePrefix "{params.outdir}/{wildcards.sample}." \
             --outSAMtype BAM SortedByCoordinate \
             {params.extra_args} > {output.log} 2> {log}
        """

rule mark_duplicates:
    input:
        bam = BAM_DIR / "{sample}.Aligned.sortedByCoord.out.bam"
    output:
        dedup_bam = BAM_DIR / "{sample}.DeDup.bam",
        flagstat = BAM_DIR / "{sample}.flagstat.txt",
        log = BAM_DIR / "{sample}.DeDup.log"
    log:
        LOG_DIR / "mark_duplicates" / "{sample}.log"
    params:
        threads = config.get("threads", 8),
        layout = lambda wildcards: SAMPLE_LAYOUT[wildcards.sample],
        mapq = MAPQ_FILTER,
        require_proper_pair = REQUIRE_PROPER_PAIR
    container: SIF
    shell:
        r"""
        mkdir -p {LOG_DIR}/mark_duplicates
        if [ "{params.layout}" = "pe" ] && [ "{params.require_proper_pair}" = "True" ]; then
            VIEW_FILTER="-f 0x2 -q {params.mapq}"
        else
            VIEW_FILTER="-q {params.mapq}"
        fi

        samtools sort -n -@ {params.threads} {input.bam} | \
            samtools view $VIEW_FILTER -b | \
            samtools fixmate -m -r - - | \
            samtools sort -@ {params.threads} - | \
            samtools markdup -r -s - {output.dedup_bam} 2> {output.log}

        samtools flagstat {output.dedup_bam} > {output.flagstat}
        samtools index {output.dedup_bam}
        cat {output.log} > {log}
        """

rule peak_calling:
    input:
        ip_bams = lambda wildcards: [BAM_DIR / f"{s}.DeDup.bam" for s in PAIR_INFO[wildcards.pair]["ip_samples"]],
        input_bams = lambda wildcards: [BAM_DIR / f"{s}.DeDup.bam" for s in PAIR_INFO[wildcards.pair]["input_samples"]],
        gtf_file = GTF
    output:
        peak_bed = OUTPUT_DIR / "exomePeak2_output_{pair}" / "peaks.bed"
    log:
        LOG_DIR / "peak_calling" / "{pair}.log"
    params:
        ip_bams = lambda wildcards: csv_bam_paths(wildcards, "ip"),
        input_bams = lambda wildcards: csv_bam_paths(wildcards, "input"),
        genome_version = config.get("genome", "hg38"),
        threads = config.get("threads", 20),
        p_cutoff = config.get("p_cutoff", 1e-10),
        strandness = config.get("strandness", "unstrand"),
        fragment_length = config.get("fragment_length", 100),
        bin_size = config.get("bin_size", 25),
        step_size = config.get("step_size", 25),
        layout = lambda wildcards: (
            "pe"
            if all(SAMPLE_LAYOUT[s] == "pe" for s in PAIR_INFO[wildcards.pair]["ip_samples"] + PAIR_INFO[wildcards.pair]["input_samples"])
            else "se"
            if all(SAMPLE_LAYOUT[s] == "se" for s in PAIR_INFO[wildcards.pair]["ip_samples"] + PAIR_INFO[wildcards.pair]["input_samples"])
            else "mixed"
        )
    container: SIF
    shell:
        """
        mkdir -p {LOG_DIR}/peak_calling
        TMPDIR=$(mktemp -d)
        cd $TMPDIR
        Rscript {workflow.basedir}/scripts/exomePeak2_call_peak.R \
            "{params.ip_bams}" \
            "{params.input_bams}" \
            {input.gtf_file} \
            "{params.genome_version}" \
            {params.threads} \
            {params.p_cutoff} \
            "{params.strandness}" \
            {params.fragment_length} \
            {params.bin_size} \
            {params.step_size} \
            "{params.layout}" 2> {log}

        if [ -f ./exomePeak2_output/peaks.bed ]; then
            mv ./exomePeak2_output/peaks.bed {output.peak_bed}
        else
            : > {output.peak_bed}
        fi
        if [ -f ./exomePeak2_result.rds ]; then
            mv ./exomePeak2_result.rds {OUTPUT_DIR}/exomePeak2_output_{wildcards.pair}/exomePeak2_result.rds
        fi
        rm -rf $TMPDIR
        """

if config.get("liftover", False):
    rule liftover_peaks:
        input:
            peak_bed = OUTPUT_DIR / "exomePeak2_output_{pair}" / "peaks.bed",
            chain = resolve_path(config["liftover_chain"])
        output:
            lifted_bed = OUTPUT_DIR / "exomePeak2_output_{pair}" / "peaks.lifted.bed",
            unmapped_bed = OUTPUT_DIR / "exomePeak2_output_{pair}" / "peaks.unmapped.bed"
        log:
            LOG_DIR / "liftover" / "{pair}.log"
        container: SIF
        shell:
            """
            mkdir -p {LOG_DIR}/liftover
            liftOver {input.peak_bed} {input.chain} {output.lifted_bed} {output.unmapped_bed} > {log} 2>&1
            """

rule merge_peaks:
    input:
        peaks = expand(OUTPUT_DIR / "exomePeak2_output_{pair}" / "peaks.bed", pair=PAIR_NAMES)
    output:
        merged = OUTPUT_DIR / "acRIP_all_peaks.bed"
    log:
        LOG_DIR / "merge_peaks.log"
    container: SIF
    shell:
        """
        mkdir -p {LOG_DIR}
        cat {input.peaks} | sort -k1,1 -k2,2n > {output.merged} 2> {log}
        """

rule motif_analysis:
    input:
        peak_bed = OUTPUT_DIR / "exomePeak2_output_{pair}" / "peaks.bed"
    output:
        streme_html = MOTIF_DIR / "streme_{pair}.html",
        peak_fa = MOTIF_DIR / "peak_{pair}.fa",
        bed6_peak = MOTIF_DIR / "bed6_peak_{pair}.bed"
    log:
        LOG_DIR / "motif_analysis" / "{pair}.log"
    params:
        genome_fasta = GENOME_FASTA,
        outdir = MOTIF_DIR
    container: SIF
    shell:
        r"""
        mkdir -p {params.outdir} {LOG_DIR}/motif_analysis

        if [ ! -s {input.peak_bed} ]; then
            echo "<html><body><h1>No peaks detected</h1><p>Motif analysis skipped because peaks.bed is empty.</p></body></html>" > {output.streme_html}
            touch {output.peak_fa} {output.bed6_peak}
            exit 0
        fi

        bed12ToBed6 -i {input.peak_bed} -n > {output.bed6_peak}

        if [ ! -f "{params.genome_fasta}.fai" ]; then
            samtools faidx {params.genome_fasta}
        fi

        bed2fasta -s -both -o {output.peak_fa} {output.bed6_peak} {params.genome_fasta}

        if [ ! -s {output.peak_fa} ]; then
            echo "<html><body><h1>No motif sequences extracted</h1><p>peak.fa is empty.</p></body></html>" > {output.streme_html}
            exit 0
        fi

        TMPDIR=$(mktemp -d)
        streme --oc $TMPDIR --rna --minw 5 --maxw 5 --thresh 0.1 --align center --p {output.peak_fa} 2> {log}
        mv $TMPDIR/streme.html {output.streme_html}
        rm -rf $TMPDIR
        """

rule scanmotif_analysis:
    input:
        peak_fa = MOTIF_DIR / "peak_{pair}.fa",
        bed6_peak = MOTIF_DIR / "bed6_peak_{pair}.bed",
        motif_meme = SCAN_MOTIF
    output:
        fimo_result = SCANMOTIF_DIR / "fimo_{pair}.tsv",
        fimo_html = SCANMOTIF_DIR / "fimo_{pair}.html",
        fimo_bed = SCANMOTIF_DIR / "fimo_{pair}_motif.bed"
    log:
        LOG_DIR / "scanmotif_analysis" / "{pair}.log"
    params:
        scan_motif_pvalue = config.get("scan_motif_pvalue", 0.01),
        outdir = SCANMOTIF_DIR
    container: SIF
    shell:
        r"""
        mkdir -p {params.outdir} {LOG_DIR}/scanmotif_analysis
        TMPDIR=$(mktemp -d)
        fimo --oc $TMPDIR --verbosity 1 --bgfile --nrdb-- --thresh {params.scan_motif_pvalue} {input.motif_meme} {input.peak_fa} 2> {log}
        mv $TMPDIR/fimo.tsv {output.fimo_result}
        mv $TMPDIR/fimo.html {output.fimo_html}
        rm -rf $TMPDIR

        # Convert FIMO tsv to BED6 and intersect with peaks
        awk 'BEGIN{{FS=OFS="\t"}} /^#/ || $1=="motif_id" || $4+0<=0 {{next}} {{print $3, $4-1, $5, $1, $2, $6}}' {output.fimo_result} > {output.fimo_bed}.tmp
        if [ -s {output.fimo_bed}.tmp ]; then
            intersectBed -s -wa -a {input.bed6_peak} -b {output.fimo_bed}.tmp | sort -k1,1 -k2,2n | uniq > {output.fimo_bed}
        else
            touch {output.fimo_bed}
        fi
        rm -f {output.fimo_bed}.tmp
        """

rule multiqc:
    input:
        fastqc = [QC_DIR / f"{s}.fastqc.log" for s in ALL_SAMPLES],
        trim = trim_log_targets(),
        star = [BAM_DIR / f"{s}.STAR.log" for s in ALL_SAMPLES],
        flagstat = [BAM_DIR / f"{s}.flagstat.txt" for s in ALL_SAMPLES]
    output:
        report = MULTIQC_DIR / "multiqc_report.html"
    log:
        LOG_DIR / "multiqc" / "multiqc.log"
    params:
        outdir = MULTIQC_DIR,
        qc_dir = QC_DIR,
        trim_dir = TRIM_DIR,
        bam_dir = BAM_DIR
    container: SIF
    shell:
        """
        mkdir -p {LOG_DIR}/multiqc
        multiqc {params.qc_dir}/* {params.trim_dir}/* {params.bam_dir}/*.STAR.log \
            {params.bam_dir}/*flagstat.txt -o {params.outdir} --force 2> {log}
        """

rule collect_statistics:
    input:
        flagstat = [BAM_DIR / f"{s}.flagstat.txt" for s in ALL_SAMPLES],
        star_logs = [BAM_DIR / f"{s}.STAR.log" for s in ALL_SAMPLES],
        trim_logs = trim_log_targets(),
        peaks = OUTPUT_DIR / "exomePeak2_output_{pair}" / "peaks.bed" if not config.get("skip_peak", False) else [],
        bed6_peak = MOTIF_DIR / "bed6_peak_{pair}.bed" if not config.get("skip_peak", False) and not config.get("skip_motif", False) else [],
        fimo_result = SCANMOTIF_DIR / "fimo_{pair}.tsv" if not config.get("skip_motif", False) else [],
        fimo_bed = SCANMOTIF_DIR / "fimo_{pair}_motif.bed" if not config.get("skip_motif", False) else []
    output:
        json = REPORT_DATA_DIR / "statistics_{pair}.json"
    log:
        LOG_DIR / "collect_statistics" / "{pair}.log"
    params:
        workdir = OUTPUT_DIR,
        ip_sample = lambda wildcards: ",".join(PAIR_INFO[wildcards.pair]["ip_samples"]),
        input_sample = lambda wildcards: ",".join(PAIR_INFO[wildcards.pair]["input_samples"]),
        peak_bed = lambda wildcards: OUTPUT_DIR / f"exomePeak2_output_{wildcards.pair}" / "peaks.bed",
        bed6_peak = lambda wildcards: MOTIF_DIR / f"bed6_peak_{wildcards.pair}.bed",
        fimo_result = lambda wildcards: SCANMOTIF_DIR / f"fimo_{wildcards.pair}.tsv",
        fimo_bed = lambda wildcards: SCANMOTIF_DIR / f"fimo_{wildcards.pair}_motif.bed"
    container: SIF
    shell:
        """
        mkdir -p {LOG_DIR}/collect_statistics
        python {workflow.basedir}/scripts/collect_stats_multi.py \
            -w {params.workdir} \
            --ip-sample "{params.ip_sample}" \
            --input-sample "{params.input_sample}" \
            --peak-bed {params.peak_bed} \
            --bed6-peak {params.bed6_peak} \
            --fimo-result {params.fimo_result} \
            --fimo-bed {params.fimo_bed} \
            -o {output.json} 2> {log}
        """

rule generate_report_figures:
    input:
        json = REPORT_DATA_DIR / "statistics_{pair}.json"
    output:
        reads_plot = REPORT_FIG_DIR / "{pair}_reads_statistics.png",
        mapping_plot = REPORT_FIG_DIR / "{pair}_mapping_rate.png",
        peak_plot = REPORT_FIG_DIR / "{pair}_peak_counts.png",
        motif_pie = REPORT_FIG_DIR / "{pair}_m6a_motif_pie.png"
    log:
        LOG_DIR / "generate_report_figures" / "{pair}.log"
    params:
        workdir = OUTPUT_DIR,
        ip_sample = lambda wildcards: ",".join(PAIR_INFO[wildcards.pair]["ip_samples"]),
        prefix = lambda wildcards: wildcards.pair
    container: SIF
    shell:
        """
        mkdir -p {LOG_DIR}/generate_report_figures
        python {workflow.basedir}/scripts/generate_figures_multi.py \
            -w {params.workdir} \
            --ip-sample "{params.ip_sample}" \
            --stats-json {input.json} \
            --prefix {params.prefix} \
            -o {REPORT_FIG_DIR} 2> {log}
        """

rule compile_report:
    input:
        json = REPORT_DATA_DIR / "statistics_{pair}.json",
        reads_plot = REPORT_FIG_DIR / "{pair}_reads_statistics.png",
        mapping_plot = REPORT_FIG_DIR / "{pair}_mapping_rate.png",
        peak_plot = REPORT_FIG_DIR / "{pair}_peak_counts.png",
        motif_pie = REPORT_FIG_DIR / "{pair}_m6a_motif_pie.png",
        multiqc = MULTIQC_DIR / "multiqc_report.html"
    output:
        html_report = REPORT_DIR / "acRIP_Analysis_Report_{pair}.html",
        html_report_en = REPORT_DIR / "acRIP_Analysis_Report_{pair}_EN.html"
    log:
        LOG_DIR / "compile_report" / "{pair}.log"
    params:
        workdir = OUTPUT_DIR,
        ip_sample = lambda wildcards: ",".join(PAIR_INFO[wildcards.pair]["ip_samples"]),
        input_sample = lambda wildcards: ",".join(PAIR_INFO[wildcards.pair]["input_samples"]),
        genome = GENOME,
        report_title = "acRIP-seq Analysis Report",
        report_subtitle = "ac4C RNA Immunoprecipitation sequencing Analysis Report",
        company_name = "pi-omicskit Bioinformatics",
        project_id = f"acRIP{datetime.now().strftime('%Y')}",
        scan_motif_pvalue = SCAN_MOTIF_PVALUE
    container: SIF
    shell:
        """
        mkdir -p {LOG_DIR}/compile_report
        python {workflow.basedir}/scripts/compile_merip_report_multi_i18n.py \
            -w {params.workdir} \
            --ip-sample "{params.ip_sample}" \
            --input-sample "{params.input_sample}" \
            -g {params.genome} \
            -t "{params.report_title}" \
            -st "{params.report_subtitle}" \
            -c "{params.company_name}" \
            -p "{params.project_id}" \
            --scan-motif-pvalue {params.scan_motif_pvalue} \
            --stats-json {input.json} \
            -l cn \
            -o {output.html_report} 2> {log}
        python {workflow.basedir}/scripts/compile_merip_report_multi_i18n.py \
            -w {params.workdir} \
            --ip-sample "{params.ip_sample}" \
            --input-sample "{params.input_sample}" \
            -g {params.genome} \
            -t "{params.report_title}" \
            -st "{params.report_subtitle}" \
            -c "{params.company_name}" \
            -p "{params.project_id}" \
            --scan-motif-pvalue {params.scan_motif_pvalue} \
            --stats-json {input.json} \
            -l en \
            -o {output.html_report_en} 2>> {log}
        """

rule collect_all_statistics:
    input:
        jsons = expand(REPORT_DATA_DIR / "statistics_{pair}.json", pair=PAIR_NAMES)
    output:
        json = REPORT_DATA_DIR / "statistics_all.json"
    log:
        LOG_DIR / "collect_all_statistics.log"
    container: SIF
    shell:
        """
        mkdir -p {LOG_DIR}
        python {workflow.basedir}/scripts/collect_all_stats.py \
            -i {input.jsons} \
            -o {output.json} 2> {log}
        """

rule generate_unified_figures:
    input:
        json = REPORT_DATA_DIR / "statistics_all.json"
    output:
        reads_plot = REPORT_FIG_DIR / "all_reads_statistics.png",
        mapping_plot = REPORT_FIG_DIR / "all_mapping_rate.png"
    log:
        LOG_DIR / "generate_unified_figures.log"
    params:
        outdir = REPORT_FIG_DIR
    container: SIF
    shell:
        """
        mkdir -p {LOG_DIR}
        python {workflow.basedir}/scripts/generate_unified_figures.py \
            --stats-json {input.json} \
            -o {params.outdir} 2> {log}
        """

rule compile_unified_report:
    input:
        json = REPORT_DATA_DIR / "statistics_all.json",
        reads_plot = REPORT_FIG_DIR / "all_reads_statistics.png",
        mapping_plot = REPORT_FIG_DIR / "all_mapping_rate.png",
        pair_peak_plots = expand(REPORT_FIG_DIR / "{pair}_peak_counts.png", pair=PAIR_NAMES),
        pair_motif_pies = expand(REPORT_FIG_DIR / "{pair}_m6a_motif_pie.png", pair=PAIR_NAMES),
        multiqc = MULTIQC_DIR / "multiqc_report.html"
    output:
        html_report = REPORT_DIR / "acRIP_Analysis_Report_CN.html",
        html_report_en = REPORT_DIR / "acRIP_Analysis_Report_EN.html"
    log:
        LOG_DIR / "compile_unified_report.log"
    params:
        workdir = OUTPUT_DIR,
        genome = GENOME,
        report_title = "acRIP-seq Analysis Report",
        report_subtitle = "ac4C RNA Immunoprecipitation sequencing Analysis Report",
        company_name = "pi-omicskit Bioinformatics",
        project_id = f"acRIP{datetime.now().strftime('%Y')}",
        scan_motif_pvalue = SCAN_MOTIF_PVALUE,
        pair_names = ",".join(PAIR_NAMES)
    container: SIF
    shell:
        """
        mkdir -p {LOG_DIR}
        python {workflow.basedir}/scripts/compile_unified_report_i18n.py \
            -w {params.workdir} \
            -g {params.genome} \
            -t "{params.report_title}" \
            -st "{params.report_subtitle}" \
            -c "{params.company_name}" \
            -p "{params.project_id}" \
            --scan-motif-pvalue {params.scan_motif_pvalue} \
            --stats-json {input.json} \
            --pair-names {params.pair_names} \
            -l cn \
            -o {output.html_report} 2> {log}
        python {workflow.basedir}/scripts/compile_unified_report_i18n.py \
            -w {params.workdir} \
            -g {params.genome} \
            -t "{params.report_title}" \
            -st "{params.report_subtitle}" \
            -c "{params.company_name}" \
            -p "{params.project_id}" \
            --scan-motif-pvalue {params.scan_motif_pvalue} \
            --stats-json {input.json} \
            --pair-names {params.pair_names} \
            -l en \
            -o {output.html_report_en} 2>> {log}
        """
