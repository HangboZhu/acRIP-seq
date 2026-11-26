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
CONSENSUS_DIR = OUTPUT_DIR / "consensus"
LOG_DIR = OUTPUT_DIR / "logs"

for dir_path in [QC_DIR, TRIM_DIR, BAM_DIR, MOTIF_DIR, SCANMOTIF_DIR,
                 MULTIQC_DIR, REPORT_DIR, REPORT_DATA_DIR, REPORT_FIG_DIR,
                 CONSENSUS_DIR, LOG_DIR]:
    dir_path.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Parse sample groups.  A group may contain multiple biological replicates.
# How those replicates reach exomePeak2 is set by `peak_pairing` further down:
# per_replicate (default) = one call per replicate (1 IP + 1 Input),
# per_pair = one joint call over all replicate BAMs.
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
        # make_config_from_pairs.py always emits a list of per-sample lists,
        # even for single-sample (no-replicate) groups: ip_fastq: [[r1, r2]].
        # Unwrap that one level so the inner [r1, r2] is treated as the
        # sample's files; a flat legacy "ip_fastq: [r1, r2]" still works.
        if (
            isinstance(raw, (list, tuple))
            and len(raw) == 1
            and isinstance(raw[0], (list, tuple))
        ):
            return [fastq_entry(raw[0])]
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

# ------------------------------------------------------------
# Peak-calling granularity  (config: peak_pairing)
# ------------------------------------------------------------
# A "sample_pairs" entry groups a biology: N IP replicates + N Input
# replicates. How many exomePeak2 runs that becomes is controlled here:
#
#   per_replicate (default)  one run per replicate, IP and Input paired by
#                            POSITION inside the group:
#                              ip_samples[0] vs input_samples[0]
#                              ip_samples[1] vs input_samples[1]  ...
#                            This is how pi-OmicsKit / OpenAc4C call peaks
#                            (exactly 1 IP + 1 Input per exomePeak2 call),
#                            so every peak set has a single, unambiguous
#                            IP/Input pair behind it.
#   per_pair                 one joint run over ALL IP reps vs ALL Input
#                            reps (exomePeak2 joint peak calling).
#
# Each run gets a peak-group id used as the {peak} wildcard and as the
# exomePeak2_output_<id> suffix. In per_pair mode the id equals the
# pair_name, so existing output paths are preserved.
PEAK_PAIRING = str(config.get("peak_pairing", "per_replicate")).lower()
if PEAK_PAIRING not in ("per_replicate", "per_pair"):
    raise ValueError(
        f"peak_pairing must be 'per_replicate' or 'per_pair', got {PEAK_PAIRING!r}"
    )

PEAK_GROUPS = []
for pair in SAMPLE_PAIRS:
    ip_list = list(pair["ip_samples"])
    input_list = list(pair["input_samples"])

    if PEAK_PAIRING == "per_pair":
        PEAK_GROUPS.append({
            "peak_id": pair["pair_name"],
            "pair_name": pair["pair_name"],
            "replicate": 0,
            "ip_samples": ip_list,
            "input_samples": input_list,
        })
        continue

    if not ip_list:
        raise ValueError(
            f"pair '{pair['pair_name']}': per_replicate peak calling needs ip_samples"
        )
    if not input_list:
        raise ValueError(
            f"pair '{pair['pair_name']}': per_replicate peak calling needs "
            f"input_samples (each IP is compared against its own Input)"
        )
    if len(ip_list) != len(input_list):
        raise ValueError(
            f"pair '{pair['pair_name']}': per_replicate peak calling pairs IP and "
            f"Input BY POSITION, so both lists must have the same length "
            f"(got {len(ip_list)} IP vs {len(input_list)} Input). "
            f"Fix the config, or set `peak_pairing: per_pair` to call peaks "
            f"jointly across replicates."
        )

    for rep, (ip_name, input_name) in enumerate(zip(ip_list, input_list), start=1):
        PEAK_GROUPS.append({
            "peak_id": f"{pair['pair_name']}_rep{rep}",
            "pair_name": pair["pair_name"],
            "replicate": rep,
            "ip_samples": [ip_name],
            "input_samples": [input_name],
        })

PEAK_IDS = [g["peak_id"] for g in PEAK_GROUPS]
PEAK_INFO = {g["peak_id"]: g for g in PEAK_GROUPS}
if len(set(PEAK_IDS)) != len(PEAK_IDS):
    _dupes = sorted({x for x in PEAK_IDS if PEAK_IDS.count(x) > 1})
    raise ValueError(f"duplicate peak-group ids: {', '.join(_dupes)}")

# Peak groups of each sample_pairs entry, in replicate order. Used by the
# majority-vote consensus step.
PAIR_PEAK_IDS = {p["pair_name"]: [] for p in SAMPLE_PAIRS}
for _g in PEAK_GROUPS:
    PAIR_PEAK_IDS[_g["pair_name"]].append(_g["peak_id"])

# ------------------------------------------------------------
# Cross-replicate consensus peaks (majority vote)
# ------------------------------------------------------------
# per_replicate peak calling produces one list per replicate, so a group has no
# single answer to "how many ac4C peaks does it have". The consensus rules
# merge the replicate lists, count how many replicates support each region
# (>= consensus_min_overlap of it covered), and keep regions supported by at
# least `consensus_min_support` replicates:
#   "all"      strict intersection - every replicate must support the region
#   "majority" floor(N/2)+1 replicates
#   <int>      that many replicates
# Groups with a single replicate are passed through unchanged.
CONSENSUS_PEAKS = config.get("consensus_peaks", True)
CONSENSUS_MIN_SUPPORT = config.get("consensus_min_support", "all")
CONSENSUS_MIN_OVERLAP = float(config.get("consensus_min_overlap", 0.5))


def consensus_min_support_for(pair_name):
    n_reps = len(PAIR_PEAK_IDS[pair_name])
    mode = str(CONSENSUS_MIN_SUPPORT).strip().lower()
    if mode == "all":
        return n_reps
    if mode == "majority":
        return n_reps // 2 + 1
    return int(CONSENSUS_MIN_SUPPORT)

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


def peak_bam_paths(wildcards, role):
    """Comma-separated DeDup BAMs of one peak group (IP or Input side)."""
    return ",".join(
        str(BAM_DIR / f"{sample}.DeDup.bam")
        for sample in PEAK_INFO[wildcards.peak][f"{role}_samples"]
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
        # Peak calling: one exomePeak2 run per peak group (1 IP + 1 Input when
        # peak_pairing == "per_replicate")
        OUTPUT_DIR / "peak_groups.tsv" if not config.get("skip_peak", False) else [],
        expand(OUTPUT_DIR / "exomePeak2_output_{peak}" / "peaks.bed", peak=PEAK_IDS) if not config.get("skip_peak", False) else [],
        # Optional liftOver harmonization
        expand(OUTPUT_DIR / "exomePeak2_output_{peak}" / "peaks.lifted.bed", peak=PEAK_IDS) if config.get("liftover", False) else [],
        # Merge enriched peaks across all peak groups
        OUTPUT_DIR / "acRIP_all_peaks.bed" if not config.get("skip_peak", False) else [],
        # Majority-vote consensus peaks per sample group
        CONSENSUS_DIR / "consensus_summary.tsv" if not config.get("skip_peak", False) and CONSENSUS_PEAKS else [],
        # Per-peak-group motif analysis (off by default for acRIP-seq)
        expand(MOTIF_DIR / "streme_{peak}.html", peak=PEAK_IDS) if not config.get("skip_peak", False) and not config.get("skip_motif", False) else [],
        expand(MOTIF_DIR / "peak_{peak}.fa", peak=PEAK_IDS) if not config.get("skip_peak", False) and not config.get("skip_motif", False) else [],
        # Per-peak-group scan motif
        expand(SCANMOTIF_DIR / "fimo_{peak}.tsv", peak=PEAK_IDS) if not config.get("skip_motif", False) else [],
        expand(SCANMOTIF_DIR / "fimo_{peak}.html", peak=PEAK_IDS) if not config.get("skip_motif", False) else [],
        expand(SCANMOTIF_DIR / "fimo_{peak}_motif.bed", peak=PEAK_IDS) if not config.get("skip_motif", False) else [],
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
    # Genome load (~35G for human/mouse) + BAM sort headroom (--limitBAMsortRAM)
    resources:
        mem_mb = config.get("star_mem_mb", 96000)
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
             --limitBAMsortRAM 40000000000 \
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

rule peak_groups_manifest:
    """Write the peak_id -> IP/Input sample map, so every peak set can be
    traced back to exactly which two samples produced it."""
    output:
        manifest = OUTPUT_DIR / "peak_groups.tsv"
    run:
        import pathlib
        pathlib.Path(str(output.manifest)).parent.mkdir(parents=True, exist_ok=True)
        with open(str(output.manifest), "w") as fh:
            fh.write(f"# peak_pairing: {PEAK_PAIRING}\n")
            fh.write("peak_id\tpair_name\treplicate\tip_sample\tinput_sample\n")
            for g in PEAK_GROUPS:
                fh.write("\t".join([
                    g["peak_id"],
                    g["pair_name"],
                    str(g["replicate"]),
                    ",".join(g["ip_samples"]),
                    ",".join(g["input_samples"]),
                ]) + "\n")


rule peak_calling:
    input:
        ip_bams = lambda wildcards: [BAM_DIR / f"{s}.DeDup.bam" for s in PEAK_INFO[wildcards.peak]["ip_samples"]],
        input_bams = lambda wildcards: [BAM_DIR / f"{s}.DeDup.bam" for s in PEAK_INFO[wildcards.peak]["input_samples"]],
        gtf_file = GTF
    output:
        peak_bed = OUTPUT_DIR / "exomePeak2_output_{peak}" / "peaks.bed"
    log:
        LOG_DIR / "peak_calling" / "{peak}.log"
    params:
        ip_bams = lambda wildcards: peak_bam_paths(wildcards, "ip"),
        input_bams = lambda wildcards: peak_bam_paths(wildcards, "input"),
        genome_version = config.get("genome", "hg38"),
        threads = config.get("threads", 20),
        p_cutoff = config.get("p_cutoff", 1e-10),
        strandness = config.get("strandness", "unstrand"),
        fragment_length = config.get("fragment_length", 100),
        bin_size = config.get("bin_size", 25),
        step_size = config.get("step_size", 25),
        layout = lambda wildcards: (
            "pe"
            if all(SAMPLE_LAYOUT[s] == "pe" for s in PEAK_INFO[wildcards.peak]["ip_samples"] + PEAK_INFO[wildcards.peak]["input_samples"])
            else "se"
            if all(SAMPLE_LAYOUT[s] == "se" for s in PEAK_INFO[wildcards.peak]["ip_samples"] + PEAK_INFO[wildcards.peak]["input_samples"])
            else "mixed"
        )
    resources:
        mem_mb = 48000
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
            mv ./exomePeak2_result.rds {OUTPUT_DIR}/exomePeak2_output_{wildcards.peak}/exomePeak2_result.rds
        fi
        rm -rf $TMPDIR
        """

if config.get("liftover", False):
    rule liftover_peaks:
        input:
            peak_bed = OUTPUT_DIR / "exomePeak2_output_{peak}" / "peaks.bed",
            chain = resolve_path(config["liftover_chain"])
        output:
            lifted_bed = OUTPUT_DIR / "exomePeak2_output_{peak}" / "peaks.lifted.bed",
            unmapped_bed = OUTPUT_DIR / "exomePeak2_output_{peak}" / "peaks.unmapped.bed"
        log:
            LOG_DIR / "liftover" / "{peak}.log"
        container: SIF
        shell:
            """
            mkdir -p {LOG_DIR}/liftover
            liftOver {input.peak_bed} {input.chain} {output.lifted_bed} {output.unmapped_bed} > {log} 2>&1
            """

rule merge_peaks:
    input:
        peaks = expand(OUTPUT_DIR / "exomePeak2_output_{peak}" / "peaks.bed", peak=PEAK_IDS)
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

if CONSENSUS_PEAKS:
    rule peak_consensus:
        """Majority vote across the replicates of one sample group."""
        input:
            peaks = lambda wildcards: [
                OUTPUT_DIR / f"exomePeak2_output_{pid}" / "peaks.bed"
                for pid in PAIR_PEAK_IDS[wildcards.pair]
            ]
        output:
            bed = CONSENSUS_DIR / "{pair}_consensus_peaks.bed",
            stats = CONSENSUS_DIR / "{pair}_peak_support.txt"
        log:
            LOG_DIR / "peak_consensus" / "{pair}.log"
        params:
            rep_names = lambda wildcards: [
                ",".join(PEAK_INFO[pid]["ip_samples"]) for pid in PAIR_PEAK_IDS[wildcards.pair]
            ],
            min_support = lambda wildcards: consensus_min_support_for(wildcards.pair),
            min_overlap = CONSENSUS_MIN_OVERLAP,
            outdir = CONSENSUS_DIR
        container: SIF
        shell:
            """
            mkdir -p {params.outdir} {LOG_DIR}/peak_consensus
            python {workflow.basedir}/scripts/consensus_peaks.py \
                --pair-name "{wildcards.pair}" \
                --peaks {input.peaks} \
                --rep-names {params.rep_names} \
                --min-support {params.min_support} \
                --min-overlap {params.min_overlap} \
                --out-bed {output.bed} \
                --out-stats {output.stats} 2> {log}
            """

    rule consensus_summary:
        input:
            stats = expand(CONSENSUS_DIR / "{pair}_peak_support.txt", pair=PAIR_NAMES)
        output:
            summary = CONSENSUS_DIR / "consensus_summary.tsv"
        log:
            LOG_DIR / "consensus_summary.log"
        container: SIF
        shell:
            """
            mkdir -p {LOG_DIR}
            awk 'FNR==1 && NR!=1 {{next}} 1' {input.stats} > {output.summary} 2> {log}
            """

rule motif_analysis:
    input:
        peak_bed = OUTPUT_DIR / "exomePeak2_output_{peak}" / "peaks.bed"
    output:
        streme_html = MOTIF_DIR / "streme_{peak}.html",
        peak_fa = MOTIF_DIR / "peak_{peak}.fa",
        bed6_peak = MOTIF_DIR / "bed6_peak_{peak}.bed"
    log:
        LOG_DIR / "motif_analysis" / "{peak}.log"
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
        peak_fa = MOTIF_DIR / "peak_{peak}.fa",
        bed6_peak = MOTIF_DIR / "bed6_peak_{peak}.bed",
        motif_meme = SCAN_MOTIF
    output:
        fimo_result = SCANMOTIF_DIR / "fimo_{peak}.tsv",
        fimo_html = SCANMOTIF_DIR / "fimo_{peak}.html",
        fimo_bed = SCANMOTIF_DIR / "fimo_{peak}_motif.bed"
    log:
        LOG_DIR / "scanmotif_analysis" / "{peak}.log"
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
        peaks = OUTPUT_DIR / "exomePeak2_output_{peak}" / "peaks.bed" if not config.get("skip_peak", False) else [],
        bed6_peak = MOTIF_DIR / "bed6_peak_{peak}.bed" if not config.get("skip_peak", False) and not config.get("skip_motif", False) else [],
        fimo_result = SCANMOTIF_DIR / "fimo_{peak}.tsv" if not config.get("skip_motif", False) else [],
        fimo_bed = SCANMOTIF_DIR / "fimo_{peak}_motif.bed" if not config.get("skip_motif", False) else []
    output:
        json = REPORT_DATA_DIR / "statistics_{peak}.json"
    log:
        LOG_DIR / "collect_statistics" / "{peak}.log"
    params:
        workdir = OUTPUT_DIR,
        ip_sample = lambda wildcards: ",".join(PEAK_INFO[wildcards.peak]["ip_samples"]),
        input_sample = lambda wildcards: ",".join(PEAK_INFO[wildcards.peak]["input_samples"]),
        peak_bed = lambda wildcards: OUTPUT_DIR / f"exomePeak2_output_{wildcards.peak}" / "peaks.bed",
        bed6_peak = lambda wildcards: MOTIF_DIR / f"bed6_peak_{wildcards.peak}.bed",
        fimo_result = lambda wildcards: SCANMOTIF_DIR / f"fimo_{wildcards.peak}.tsv",
        fimo_bed = lambda wildcards: SCANMOTIF_DIR / f"fimo_{wildcards.peak}_motif.bed"
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
        json = REPORT_DATA_DIR / "statistics_{peak}.json"
    output:
        reads_plot = REPORT_FIG_DIR / "{peak}_reads_statistics.png",
        mapping_plot = REPORT_FIG_DIR / "{peak}_mapping_rate.png",
        peak_plot = REPORT_FIG_DIR / "{peak}_peak_counts.png",
        motif_pie = REPORT_FIG_DIR / "{peak}_m6a_motif_pie.png"
    log:
        LOG_DIR / "generate_report_figures" / "{peak}.log"
    params:
        workdir = OUTPUT_DIR,
        ip_sample = lambda wildcards: ",".join(PEAK_INFO[wildcards.peak]["ip_samples"]),
        prefix = lambda wildcards: wildcards.peak
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
        json = REPORT_DATA_DIR / "statistics_{peak}.json",
        reads_plot = REPORT_FIG_DIR / "{peak}_reads_statistics.png",
        mapping_plot = REPORT_FIG_DIR / "{peak}_mapping_rate.png",
        peak_plot = REPORT_FIG_DIR / "{peak}_peak_counts.png",
        motif_pie = REPORT_FIG_DIR / "{peak}_m6a_motif_pie.png",
        multiqc = MULTIQC_DIR / "multiqc_report.html"
    output:
        html_report = REPORT_DIR / "acRIP_Analysis_Report_{peak}.html",
        html_report_en = REPORT_DIR / "acRIP_Analysis_Report_{peak}_EN.html"
    log:
        LOG_DIR / "compile_report" / "{peak}.log"
    params:
        workdir = OUTPUT_DIR,
        ip_sample = lambda wildcards: ",".join(PEAK_INFO[wildcards.peak]["ip_samples"]),
        input_sample = lambda wildcards: ",".join(PEAK_INFO[wildcards.peak]["input_samples"]),
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
            -g "{params.genome}" \
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
            -g "{params.genome}" \
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
        jsons = expand(REPORT_DATA_DIR / "statistics_{peak}.json", peak=PEAK_IDS)
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
        pair_peak_plots = expand(REPORT_FIG_DIR / "{peak}_peak_counts.png", peak=PEAK_IDS),
        pair_motif_pies = expand(REPORT_FIG_DIR / "{peak}_m6a_motif_pie.png", peak=PEAK_IDS),
        multiqc = MULTIQC_DIR / "multiqc_report.html"
    output:
        html_report = REPORT_DIR / "acRIP_Analysis_Report_CN.html",
        html_report_en = REPORT_DIR / "acRIP_Analysis_Report_EN.html"
    log:
        LOG_DIR / "compile_unified_report.log"
    params:
        workdir = OUTPUT_DIR,
        # -g is a plain HTML label here (no BSgenome lookup), but argparse
        # rejects an empty value ("expected one argument", job 1872919/1872923),
        # so substitute an honest label when GC correction is disabled.
        genome = GENOME if GENOME else "na (no BSgenome; GC correction disabled)",
        report_title = "acRIP-seq Analysis Report",
        report_subtitle = "ac4C RNA Immunoprecipitation sequencing Analysis Report",
        company_name = "pi-omicskit Bioinformatics",
        project_id = f"acRIP{datetime.now().strftime('%Y')}",
        scan_motif_pvalue = SCAN_MOTIF_PVALUE,
        pair_names = ",".join(PEAK_IDS)
    container: SIF
    shell:
        """
        mkdir -p {LOG_DIR}
        python {workflow.basedir}/scripts/compile_unified_report_i18n.py \
            -w {params.workdir} \
            -g "{params.genome}" \
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
            -g "{params.genome}" \
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
