#!/usr/bin/env Rscript

# ---------------------------------------------------------------------------
# exomePeak2 peak calling for acRIP-seq
#
# OpenAc4C unified acRIP-seq pipeline:
#   - joint peak calling when biological replicates are present
#     (all replicate BAMs are passed together as vectors)
#   - GC content bias correction (enabled whenever `genome` is provided)
#   - unified peak calling cut-off P < 1e-10
#
# Usage:
#   Rscript exomePeak2_call_peak.R <ip_bams_csv> <input_bams_csv> <gtf> \
#       <genome> <threads> [p_cutoff] [strandness] [fragment_length] \
#       [bin_size] [step_size] [layout]
#
# ip_bams_csv / input_bams_csv : comma-separated BAM paths.
# genome : UCSC genome name (e.g. hg38, mm10) or "" to disable GC correction.
# layout : pe | se | mixed (used only by exomePeak2 versions that expose
#          a `paired_end` argument; newer versions auto-handle both).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(exomePeak2))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop(
    "Usage: Rscript exomePeak2_call_peak.R <ip_bams_csv> <input_bams_csv> ",
    "<gtf> <genome> <threads> [p_cutoff] [strandness] ",
    "[fragment_length] [bin_size] [step_size] [layout]"
  )
}

split_bams <- function(x) {
  if (is.null(x) || !nzchar(x)) {
    return(NULL)
  }
  bams <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  bams <- bams[nzchar(bams)]
  if (length(bams) == 0) NULL else bams
}

ip_bams <- split_bams(args[1])
input_bams <- split_bams(args[2])
gene_anno_gtf <- args[3]
genome_version <- args[4]
threads <- as.numeric(args[5])
p_cutoff <- if (length(args) >= 6 && nzchar(args[6])) as.numeric(args[6]) else 1e-10
strandness <- if (length(args) >= 7 && nzchar(args[7])) args[7] else "unstrand"
fragment_length <- if (length(args) >= 8 && nzchar(args[8])) as.numeric(args[8]) else 100
bin_size <- if (length(args) >= 9 && nzchar(args[9])) as.numeric(args[9]) else 25
step_size <- if (length(args) >= 10 && nzchar(args[10])) as.numeric(args[10]) else 25
layout <- if (length(args) >= 11 && nzchar(args[11])) args[11] else "auto"

if (is.null(ip_bams)) {
  stop("At least one IP BAM file is required.")
}

# Empty genome string disables GC correction (for species without a BSgenome).
genome_arg <- if (nzchar(genome_version)) genome_version else NULL

cat("Running exomePeak2 with parameters:\n")
cat("IP BAM(s):", paste(ip_bams, collapse = ", "), "\n")
cat("Input BAM(s):", if (is.null(input_bams)) "none" else paste(input_bams, collapse = ", "), "\n")
cat("GTF:", gene_anno_gtf, "\n")
cat("Genome:", if (is.null(genome_arg)) "none (GC correction disabled)" else genome_arg, "\n")
cat("P cutoff:", p_cutoff, "\n")
cat("Strandness:", strandness, "\n")
cat("Fragment length:", fragment_length, "\n")
cat("Bin size:", bin_size, "\n")
cat("Step size:", step_size, "\n")
cat("Threads:", threads, "\n")
cat("Library layout:", layout, "\n")
cat("Replicates -> joint peak calling:", length(ip_bams) > 1, "\n")

ep_args <- list(
  bam_ip = ip_bams,
  bam_input = input_bams,
  gff = gene_anno_gtf,
  genome = genome_arg,
  strandness = strandness,
  fragment_length = fragment_length,
  bin_size = bin_size,
  step_size = step_size,
  p_cutoff = p_cutoff,
  mode = "exon",
  parallel = threads,
  save_output = TRUE,
  save_dir = getwd(),
  experiment_name = "exomePeak2_output"
)

# exomePeak2 v1.14+ auto-handles fragments/singletons and has no
# `paired_end` argument; older versions expose it as a global switch.
if ("paired_end" %in% names(formals(exomePeak2))) {
  if (layout == "mixed") {
    warning(
      "Mixed SE/PE samples in one group; using paired_end = TRUE. ",
      "Consider splitting mixed groups for older exomePeak2 versions."
    )
    ep_args$paired_end <- TRUE
  } else if (layout %in% c("pe", "se")) {
    ep_args$paired_end <- (layout == "pe")
  }
}

result <- do.call(exomePeak2, ep_args)

saveRDS(result, file = "exomePeak2_result.rds")
cat("Peak calling completed successfully!\n")
