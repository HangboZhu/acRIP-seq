#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Majority-vote consensus peaks across the replicates of one IP/Input group.

Why this exists
---------------
With ``peak_pairing: per_replicate`` every replicate gets its own exomePeak2
run, so a group of N replicates yields N independent peak lists. Those lists
overlap but never coincide, so the group has no single "this is our peak set"
answer. This script produces one, by majority vote.

Method
------
1. Collapse the peaks of all replicates into candidate regions: intervals that
   overlap or touch are merged (bedtools-merge semantics).
2. For every candidate region count the replicates that have a peak covering at
   least ``--min-overlap`` (default 0.5) of that region.
3. Keep regions whose support is >= ``--min-support``.

A replicate with a single data point is passed through unchanged (min-support
1), so un-replicated groups behave exactly as before.

Outputs
-------
``<pair>_consensus_peaks.bed``
    BED6 plus the support count in column 5. Strand is the common strand of the
    supporting peaks, or ``.`` when they disagree.
``<pair>_peak_support.txt``
    One TSV row summarising the vote for this group.
"""

import argparse
import sys
from collections import defaultdict


def read_bed(path):
    """chrom -> sorted list of (start, end, strand)."""
    intervals = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                continue
            if end <= start:
                continue
            strand = fields[5] if len(fields) > 5 and fields[5] in "+-" else "."
            intervals[fields[0]].append((start, end, strand))
    for chrom in intervals:
        intervals[chrom].sort()
    return intervals


def merge_intervals(records):
    """Merge overlapping/touching intervals, keeping the set of strands."""
    merged = []
    for start, end, strand in records:
        if merged and start <= merged[-1][1]:
            if end > merged[-1][1]:
                merged[-1][1] = end
            merged[-1][2].add(strand)
        else:
            merged.append([start, end, {strand}])
    return merged


def supported(regions, rep_intervals, min_overlap):
    """Boolean per region: does this replicate cover >= min_overlap of it?"""
    result = [False] * len(regions)
    n = len(rep_intervals)
    j = 0
    for i, (region_start, region_end, _strands) in enumerate(regions):
        while j < n and rep_intervals[j][1] <= region_start:
            j += 1
        length = region_end - region_start
        k = j
        while k < n and rep_intervals[k][0] < region_end:
            overlap = min(rep_intervals[k][1], region_end) - max(rep_intervals[k][0], region_start)
            if overlap > 0 and overlap / length >= min_overlap:
                result[i] = True
                break
            k += 1
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pair-name", required=True)
    parser.add_argument("--peaks", nargs="+", required=True, help="One peaks.bed per replicate")
    parser.add_argument("--rep-names", nargs="+", required=True, help="Replicate labels, same order as --peaks")
    parser.add_argument("--min-support", required=True, help="Integer, or 'majority'")
    parser.add_argument("--min-overlap", type=float, default=0.5)
    parser.add_argument("--out-bed", required=True)
    parser.add_argument("--out-stats", required=True)
    args = parser.parse_args()

    if len(args.peaks) != len(args.rep_names):
        sys.exit("--peaks and --rep-names must have the same length")

    rep_beds = [read_bed(p) for p in args.peaks]
    rep_counts = [sum(len(v) for v in bed.values()) for bed in rep_beds]
    n_reps = len(rep_beds)

    if str(args.min_support).lower() == "majority":
        min_support = n_reps // 2 + 1
    else:
        min_support = int(args.min_support)

    # Single replicate (or a group where voting is meaningless): pass through.
    if n_reps == 1:
        kept = []
        for chrom, records in sorted(rep_beds[0].items()):
            for start, end, strand in records:
                kept.append((chrom, start, end, strand, 1))
        support_hist = {1: len(kept)}
    else:
        candidate = {}
        for bed in rep_beds:
            for chrom, records in bed.items():
                candidate.setdefault(chrom, []).extend(records)

        support_hist = defaultdict(int)
        kept = []
        for chrom in sorted(candidate):
            regions = merge_intervals(sorted(candidate[chrom]))
            masks = [supported(regions, rep_beds[r].get(chrom, []), args.min_overlap) for r in range(n_reps)]
            for idx, (start, end, strands) in enumerate(regions):
                support = sum(1 for mask in masks if mask[idx])
                support_hist[support] += 1
                if support < min_support:
                    continue
                # Strand: the common one, or '.' when the replicates disagree.
                contributing = {s for r in range(n_reps) if masks[r][idx] for s in strands}
                strand = next(iter(contributing)) if len(contributing) == 1 else "."
                kept.append((chrom, start, end, strand, support))

    with open(args.out_bed, "w") as fh:
        for i, (chrom, start, end, strand, support) in enumerate(kept, start=1):
            fh.write(f"{chrom}\t{start}\t{end}\t{args.pair_name}_consensus_{i}\t{support}\t{strand}\n")

    union_regions = sum(support_hist.values())
    hist_text = ",".join(f"{k}:{support_hist.get(k, 0)}" for k in range(1, n_reps + 1))
    fraction = len(kept) / union_regions if union_regions else 0.0

    with open(args.out_stats, "w") as fh:
        fh.write(
            "pair_name\tn_reps\tmin_support\tmin_overlap\trep_names\trep_peak_counts\t"
            "union_regions\tsupport_hist\tconsensus_peaks\tfraction_of_union\n"
        )
        fh.write("\t".join([
            args.pair_name,
            str(n_reps),
            str(min_support),
            f"{args.min_overlap:g}",
            ",".join(args.rep_names),
            ",".join(str(c) for c in rep_counts),
            str(union_regions),
            hist_text,
            str(len(kept)),
            f"{fraction:.3f}",
        ]) + "\n")

    print(f"{args.pair_name}: {n_reps} replicates -> {union_regions} union regions, "
          f"{len(kept)} consensus peaks (support >= {min_support})")


if __name__ == "__main__":
    main()
