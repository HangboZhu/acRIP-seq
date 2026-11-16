#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
acRIP-seq Analysis Figures Generator (Multi-pair version)
生成acRIP-seq分析的图表用于报告
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
    with open(stats_file, 'r') as f:
        return json.load(f)


def generate_reads_barplot(stats, output_dir, prefix):
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
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    bars1 = ax.bar(x - width, total_reads, width, label='Total Reads', color='#4299e1')
    bars2 = ax.bar(x, trimmed_reads, width, label='Trimmed Reads', color='#48bb78')
    bars3 = ax.bar(x + width, mapped_reads, width, label='Mapped Reads', color='#ed8936')
    
    ax.set_xlabel('Sample', fontsize=12, fontweight='bold')
    ax.set_ylabel('Reads (Million)', fontsize=12, fontweight='bold')
    ax.set_title('Reads Statistics by Sample', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(sample_names, rotation=15, ha='right')
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


def generate_mapping_rate_plot(stats, output_dir, prefix):
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
    
    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.bar(sample_names, mapping_rates, color=colors, edgecolor='black', linewidth=0.5)
    
    ax.set_xlabel('Sample', fontsize=12, fontweight='bold')
    ax.set_ylabel('Mapping Rate (%)', fontsize=12, fontweight='bold')
    ax.set_title('Mapping Rate by Sample', fontsize=14, fontweight='bold')
    ax.set_ylim(0, 100)
    ax.set_xticklabels(sample_names, rotation=15, ha='right')
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


def generate_peak_count_plot(stats, output_dir, ip_sample, prefix):
    output_dir = Path(output_dir)
    
    peak_count_exome = stats.get("peak_count_exomepeak2", 0)
    peak_count_bed6 = stats.get("peak_count_bed6", 0)
    
    if peak_count_exome == 0 and peak_count_bed6 == 0:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, 'No Peaks Detected',
                ha='center', va='center', fontsize=14, fontweight='bold')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        plt.tight_layout()
        output_file = output_dir / f"{prefix}_peak_counts.png"
        plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close()
        return output_file
    
    categories = ['ExomePeak2\nPeaks', 'BED6 Format\nPeaks']
    values = [peak_count_exome, peak_count_bed6]
    colors = ['#4299e1', '#48bb78']
    
    fig, ax = plt.subplots(figsize=(8, 6))
    bars = ax.bar(categories, values, color=colors, edgecolor='black', linewidth=0.5)
    
    ax.set_ylabel('Peak Count', fontsize=12, fontweight='bold')
    ax.set_title(f'Peak Calling Results\n({ip_sample})', fontsize=14, fontweight='bold')
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{int(height):,}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    output_file = output_dir / f"{prefix}_peak_counts.png"
    plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    
    return output_file


def generate_m6a_motif_pie(stats, output_dir, prefix):
    output_dir = Path(output_dir)
    
    m6a_peaks = stats.get("m6a_motif_peaks", 0)
    total_peaks = stats.get("peak_count_bed6", 0)
    
    if total_peaks == 0:
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.text(0.5, 0.5, 'No Motif Data Available',
                ha='center', va='center', fontsize=14, fontweight='bold')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        plt.tight_layout()
        output_file = output_dir / f"{prefix}_m6a_motif_pie.png"
        plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close()
        return output_file
    
    non_m6a_peaks = total_peaks - m6a_peaks
    labels = [f'ac4C Motif Peaks\n({m6a_peaks:,})', f'Other Peaks\n({non_m6a_peaks:,})']
    sizes = [m6a_peaks, non_m6a_peaks]
    colors = ['#4299e1', '#a0aec0']
    explode = (0.05, 0)
    
    fig, ax = plt.subplots(figsize=(8, 8))
    wedges, texts, autotexts = ax.pie(sizes, explode=explode, labels=labels, colors=colors,
                                       autopct='%1.1f%%', shadow=False, startangle=90,
                                       textprops={'fontsize': 11, 'fontweight': 'bold'})
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontsize(12)
        autotext.set_fontweight('bold')
    
    ax.set_title(f'ac4C Motif Distribution\n(Total: {total_peaks:,} peaks)',
                 fontsize=14, fontweight='bold', pad=20)
    
    plt.tight_layout()
    output_file = output_dir / f"{prefix}_m6a_motif_pie.png"
    plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    
    return output_file


def main():
    parser = argparse.ArgumentParser(description='Generate acRIP-seq Report Figures')
    parser.add_argument('-w', '--workdir', default='result', help='Working directory')
    parser.add_argument('--ip-sample', required=True, help='IP sample name')
    parser.add_argument('--stats-json', required=True, help='Input statistics JSON file')
    parser.add_argument('--prefix', required=True, help='Output file prefix')
    parser.add_argument('-o', '--outdir', default=None, help='Output directory for figures')
    args = parser.parse_args()
    
    workdir = Path(args.workdir)
    ip_sample = args.ip_sample
    stats_file = Path(args.stats_json)
    prefix = args.prefix
    outdir = Path(args.outdir) if args.outdir else workdir / "report" / "figures"
    
    if not stats_file.exists():
        print(f"Error: Statistics file not found: {stats_file}")
        sys.exit(1)
    
    stats = load_stats(stats_file)
    outdir.mkdir(parents=True, exist_ok=True)
    
    print(f"Generating figures in {outdir} with prefix {prefix}...")
    
    generate_reads_barplot(stats, outdir, prefix)
    print(f"  - {prefix}_reads_statistics.png")
    
    generate_mapping_rate_plot(stats, outdir, prefix)
    print(f"  - {prefix}_mapping_rate.png")
    
    generate_peak_count_plot(stats, outdir, ip_sample, prefix)
    print(f"  - {prefix}_peak_counts.png")
    
    generate_m6a_motif_pie(stats, outdir, prefix)
    print(f"  - {prefix}_m6a_motif_pie.png")
    
    print("Figure generation completed!")


if __name__ == "__main__":
    main()
