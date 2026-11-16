#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate unified acRIP-seq report figures across all samples/pairs
"""

import os
import sys
import json
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
plt.rcParams['axes.unicode_minus'] = False


def load_stats(stats_file):
    with open(stats_file, 'r', encoding='utf-8') as f:
        return json.load(f)


def generate_reads_barplot(stats, output_dir, prefix="all"):
    output_dir = Path(output_dir)
    samples = list(stats.get("samples", {}).keys())

    if not samples:
        return None

    sample_names = []
    total_reads = []
    trimmed_reads = []
    mapped_reads = []

    for sample_name in samples:
        sample_data = stats["samples"][sample_name]
        sample_names.append(sample_name)
        trim_data = sample_data.get("trimming", {})
        align_data = sample_data.get("alignment", {})
        flagstat_data = sample_data.get("flagstat", {})

        total = trim_data.get("total_reads", 0)
        trimmed = align_data.get("input_reads", trim_data.get("trimmed_reads", total))
        mapped = align_data.get("uniquely_mapped", flagstat_data.get("mapped_reads", 0))

        total_reads.append(total / 1e6 if total else 0)
        trimmed_reads.append(trimmed / 1e6 if trimmed else 0)
        mapped_reads.append(mapped / 1e6 if mapped else 0)

    x = np.arange(len(sample_names))
    width = 0.25

    fig, ax = plt.subplots(figsize=(max(12, len(sample_names) * 1.2), 6))

    bars1 = ax.bar(x - width, total_reads, width, label='Total Reads', color='#4299e1')
    bars2 = ax.bar(x, trimmed_reads, width, label='Trimmed Reads', color='#48bb78')
    bars3 = ax.bar(x + width, mapped_reads, width, label='Mapped Reads', color='#ed8936')

    ax.set_xlabel('Sample', fontsize=12, fontweight='bold')
    ax.set_ylabel('Reads (Million)', fontsize=12, fontweight='bold')
    ax.set_title('Reads Statistics by Sample (All Pairs)', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(sample_names, rotation=30, ha='right')
    ax.legend(loc='upper right')
    ax.grid(axis='y', alpha=0.3, linestyle='--')

    for bars in [bars1, bars2, bars3]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax.annotate(f'{height:.1f}',
                            xy=(bar.get_x() + bar.get_width() / 2, height),
                            xytext=(0, 3),
                            textcoords="offset points",
                            ha='center', va='bottom', fontsize=8)

    plt.tight_layout()
    output_file = output_dir / f"{prefix}_reads_statistics.png"
    plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()

    return output_file


def generate_mapping_rate_plot(stats, output_dir, prefix="all"):
    output_dir = Path(output_dir)
    samples = list(stats.get("samples", {}).keys())

    if not samples:
        return None

    sample_names = []
    mapping_rates = []
    colors = []

    for sample_name in samples:
        sample_data = stats["samples"][sample_name]
        sample_names.append(sample_name)
        flagstat_data = sample_data.get("flagstat", {})
        align_data = sample_data.get("alignment", {})
        rate = align_data.get("mapped_rate", flagstat_data.get("mapping_rate", 0))
        mapping_rates.append(rate)
        if rate >= 70:
            colors.append('#48bb78')
        elif rate >= 50:
            colors.append('#ecc94b')
        else:
            colors.append('#f56565')

    fig, ax = plt.subplots(figsize=(max(10, len(sample_names) * 1.0), 6))
    bars = ax.bar(sample_names, mapping_rates, color=colors, edgecolor='black', linewidth=0.5)

    ax.set_xlabel('Sample', fontsize=12, fontweight='bold')
    ax.set_ylabel('Mapping Rate (%)', fontsize=12, fontweight='bold')
    ax.set_title('Mapping Rate by Sample (All Pairs)', fontsize=14, fontweight='bold')
    ax.set_ylim(0, 100)
    ax.set_xticklabels(sample_names, rotation=30, ha='right')
    ax.grid(axis='y', alpha=0.3, linestyle='--')

    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height:.1f}%',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=10, fontweight='bold')

    ax.axhline(y=70, color='#3182ce', linestyle='--', alpha=0.5, label='Good (>70%)')

    plt.tight_layout()
    output_file = output_dir / f"{prefix}_mapping_rate.png"
    plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()

    return output_file


def main():
    parser = argparse.ArgumentParser(description='Generate Unified acRIP-seq Report Figures')
    parser.add_argument('--stats-json', required=True, help='Input unified statistics JSON file')
    parser.add_argument('-o', '--outdir', default=None, help='Output directory for figures')
    args = parser.parse_args()

    stats_file = Path(args.stats_json)
    outdir = Path(args.outdir) if args.outdir else Path("result/report/figures")

    if not stats_file.exists():
        print(f"Error: Statistics file not found: {stats_file}")
        sys.exit(1)

    stats = load_stats(stats_file)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Generating unified figures in {outdir}...")

    generate_reads_barplot(stats, outdir, prefix="all")
    print("  - all_reads_statistics.png")

    generate_mapping_rate_plot(stats, outdir, prefix="all")
    print("  - all_mapping_rate.png")

    print("Unified figure generation completed!")


if __name__ == "__main__":
    main()
