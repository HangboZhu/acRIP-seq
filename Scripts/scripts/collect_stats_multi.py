#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
acRIP-seq Analysis Statistics Collector (Multi-pair version)
收集acRIP-seq分析的统计数据用于报告生成
"""

import os
import sys
import json
import re
import argparse
from pathlib import Path
from datetime import datetime


def parse_flagstat(flagstat_file):
    """解析samtools flagstat输出文件"""
    stats = {
        "total_reads": 0,
        "mapped_reads": 0,
        "mapping_rate": 0.0
    }
    
    if not os.path.exists(flagstat_file):
        return stats
    
    with open(flagstat_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if "in total" in line or "QC-passed reads" in line:
                match = re.search(r'(\d+)', line)
                if match:
                    stats["total_reads"] = int(match.group(1))
            elif "mapped (" in line:
                match = re.search(r'(\d+)\s+\+\s+\d+\s+mapped', line)
                if match:
                    stats["mapped_reads"] = int(match.group(1))
                rate_match = re.search(r'mapped \((.+?)%:', line)
                if rate_match:
                    stats["mapping_rate"] = float(rate_match.group(1))
    
    return stats


def parse_star_log(log_file):
    """解析STAR比对日志"""
    stats = {
        "input_reads": 0,
        "uniquely_mapped": 0,
        "mapped_rate": 0.0
    }
    
    if not os.path.exists(log_file):
        return stats
    
    with open(log_file, 'r') as f:
        content = f.read()
        
        match = re.search(r'Number of input reads \|\s+(\d+)', content)
        if match:
            stats["input_reads"] = int(match.group(1))
        
        match = re.search(r'Uniquely mapped reads number \|\s+(\d+)', content)
        if match:
            stats["uniquely_mapped"] = int(match.group(1))
        
        match = re.search(r'Uniquely mapped reads % \|\s+([\d.]+)%', content)
        if match:
            stats["mapped_rate"] = float(match.group(1))
    
    return stats


def parse_trim_log(log_file):
    """解析trim_galore日志"""
    stats = {
        "total_reads": 0,
        "trimmed_reads": 0
    }
    
    if not os.path.exists(log_file):
        return stats
    
    with open(log_file, 'r') as f:
        content = f.read()
        lines = content.split('\n')
        for line in lines:
            if 'reads processed' in line.lower():
                match = re.search(r'(\d+)', line.replace(',', ''))
                if match:
                    stats["total_reads"] = int(match.group(1))
            elif 'reads written' in line.lower():
                match = re.search(r'(\d+)', line.replace(',', ''))
                if match:
                    stats["trimmed_reads"] = int(match.group(1))
    
    if stats["trimmed_reads"] == 0 and stats["total_reads"] > 0:
        stats["trimmed_reads"] = stats["total_reads"]
    
    return stats


def parse_trimming_report(report_file):
    """解析trim_galore trimming report文件"""
    stats = {
        "total_reads": 0,
        "trimmed_reads": 0
    }
    
    if not os.path.exists(report_file):
        return stats
    
    with open(report_file, 'r') as f:
        content = f.read()
        match = re.search(r'Total reads processed:\s+([\d,]+)', content)
        if match:
            stats["total_reads"] = int(match.group(1).replace(',', ''))
        match = re.search(r'Reads written \(passing filters\):\s+([\d,]+)', content)
        if match:
            stats["trimmed_reads"] = int(match.group(1).replace(',', ''))
    
    return stats


def parse_fastqc_html(html_file):
    """解析FastQC HTML报告中的Total Sequences"""
    if not os.path.exists(html_file):
        return 0
    
    with open(html_file, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
    # FastQC HTML中 Total Sequences 通常形如: <td>Total Sequences</td><td>50000000</td>
    # 先去掉内容中的换行符以便跨行匹配
    content_oneline = re.sub(r'\s+', ' ', content)
    match = re.search(r'Total Sequences\s*</td>\s*<td>(\d+)</td>', content_oneline, re.IGNORECASE)
    if match:
        return int(match.group(1))
    
    # 备用：直接搜索 Total Sequences 后面的数字
    match = re.search(r'Total Sequences.*?<td[^>]*>(\d+)</td>', content_oneline, re.IGNORECASE)
    if match:
        return int(match.group(1))
    
    return 0


def count_bed_lines(bed_file):
    """统计BED文件的行数"""
    if not os.path.exists(bed_file):
        return 0
    
    count = 0
    with open(bed_file, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                count += 1
    return count


def split_sample_names(value):
    """Split a comma-separated sample list (supports single names)."""
    if not value:
        return []
    return [name.strip() for name in str(value).split(',') if name.strip()]


def collect_statistics(workdir, ip_sample, input_sample, peak_bed, bed6_peak, fimo_result, fimo_bed=None, ip_trim_report=None, input_trim_report=None):
    """收集所有统计数据（支持逗号分隔的多个重复样本）"""
    workdir = Path(workdir)
    
    qc_dir = workdir / "qc"
    trim_dir = workdir / "trim"
    bam_dir = workdir / "bam"
    motif_dir = workdir / "motif_analysis"
    scanmotif_dir = workdir / "scan_motif"
    
    ip_samples = split_sample_names(ip_sample)
    input_samples = split_sample_names(input_sample)
    
    stats = {
        "ip_sample": ip_sample,
        "input_sample": input_sample,
        "collection_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "samples": {}
    }
    
    def get_sample_trimming_stats(sample_name, trim_report=None):
        """获取单个样本的trimming统计，优先从fastqc html读取total reads"""
        trimming_stats = {"total_reads": 0, "trimmed_reads": 0}
        
        # 优先从fastqc html读取Total Sequences作为total_reads
        # paired-end: <sample>_R1_fastqc.html; single-end: <sample>_fastqc.html
        fastqc_candidates = [
            qc_dir / f"{sample_name}_R1_fastqc.html",
            qc_dir / f"{sample_name}_fastqc.html",
            qc_dir / f"{sample_name}_R2_fastqc.html",
        ]
        for fastqc_html in fastqc_candidates:
            fastqc_reads = parse_fastqc_html(fastqc_html)
            if fastqc_reads > 0:
                trimming_stats["total_reads"] = fastqc_reads
                break
        
        # 再尝试从trim report/log读取
        if trim_report and os.path.exists(trim_report):
            report_stats = parse_trimming_report(trim_report)
            if trimming_stats["total_reads"] == 0:
                trimming_stats["total_reads"] = report_stats.get("total_reads", 0)
            trimming_stats["trimmed_reads"] = report_stats.get("trimmed_reads", 0)
        else:
            trim_candidates = [
                trim_dir / f"{sample_name}.se.trim.log",
                trim_dir / f"{sample_name}.pe.trim.log",
                trim_dir / f"{sample_name}.trim.log",
            ]
            trim_log = next((p for p in trim_candidates if p.exists()), None)
            if trim_log is not None:
                log_stats = parse_trim_log(trim_log)
                if trimming_stats["total_reads"] == 0:
                    trimming_stats["total_reads"] = log_stats.get("total_reads", 0)
                trimming_stats["trimmed_reads"] = log_stats.get("trimmed_reads", 0)
        
        # 如果没有trimmed_reads，用total_reads兜底
        if trimming_stats["trimmed_reads"] == 0 and trimming_stats["total_reads"] > 0:
            trimming_stats["trimmed_reads"] = trimming_stats["total_reads"]
        
        return trimming_stats

    def add_sample(sample_name, sample_type, trim_report=None):
        sample_stats = {
            "type": sample_type,
            "fastqc": {},
            "trimming": get_sample_trimming_stats(sample_name, trim_report),
            "alignment": {}
        }

        star_log = bam_dir / f"{sample_name}.Log.final.out"
        if star_log.exists():
            sample_stats["alignment"] = parse_star_log(star_log)

        flagstat_file = bam_dir / f"{sample_name}.flagstat.txt"
        if flagstat_file.exists():
            sample_stats["flagstat"] = parse_flagstat(flagstat_file)

        stats["samples"][sample_name] = sample_stats

    for name in ip_samples:
        add_sample(name, "IP", ip_trim_report)
    for name in input_samples:
        add_sample(name, "Input", input_trim_report)
    
    stats["peak_count_exomepeak2"] = count_bed_lines(peak_bed)
    stats["peak_count_bed6"] = count_bed_lines(bed6_peak)
    
    # 使用fimo_bed文件行数作为motif peaks数量
    stats["m6a_motif_peaks"] = count_bed_lines(fimo_bed) if fimo_bed else 0
    
    if stats["peak_count_bed6"] > 0:
        stats["m6a_motif_ratio"] = stats["m6a_motif_peaks"] / stats["peak_count_bed6"] * 100
    else:
        stats["m6a_motif_ratio"] = 0.0
    
    return stats


def main():
    parser = argparse.ArgumentParser(description='Collect acRIP-seq Statistics')
    parser.add_argument('-w', '--workdir', default='result', help='Working directory')
    parser.add_argument('--ip-sample', required=True, help='IP sample name')
    parser.add_argument('--input-sample', default='', help='Input sample name (optional)')
    parser.add_argument('--ip-trim-report', default='', help='IP sample trimming report path')
    parser.add_argument('--input-trim-report', default='', help='Input sample trimming report path')
    parser.add_argument('--peak-bed', required=True, help='Peak BED file path')
    parser.add_argument('--bed6-peak', required=True, help='BED6 peak file path')
    parser.add_argument('--fimo-result', required=True, help='FIMO result file path')
    parser.add_argument('--fimo-bed', default='', help='FIMO motif BED file path')
    parser.add_argument('-o', '--output', default='result/report/data/statistics.json', help='Output JSON file')
    args = parser.parse_args()
    
    workdir = args.workdir
    ip_sample = args.ip_sample
    input_sample = args.input_sample if args.input_sample else None
    ip_trim_report = args.ip_trim_report if args.ip_trim_report else None
    input_trim_report = args.input_trim_report if args.input_trim_report else None
    peak_bed = args.peak_bed
    bed6_peak = args.bed6_peak
    fimo_result = args.fimo_result
    fimo_bed = args.fimo_bed if args.fimo_bed else None
    output_json = Path(args.output)
    
    print(f"Collecting statistics from {workdir}...")
    print(f"  IP sample: {ip_sample}")
    print(f"  Input sample: {input_sample}")
    print(f"  Peak BED: {peak_bed}")
    
    stats = collect_statistics(workdir, ip_sample, input_sample, peak_bed, bed6_peak, fimo_result, fimo_bed, ip_trim_report, input_trim_report)
    
    output_json.parent.mkdir(parents=True, exist_ok=True)
    with open(output_json, 'w', encoding='utf-8') as f:
        json.dump(stats, f, indent=2, ensure_ascii=False)
    
    print(f"Statistics saved to: {output_json}")


if __name__ == "__main__":
    main()
