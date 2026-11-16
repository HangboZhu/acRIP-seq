#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unified acRIP-seq Analysis HTML Report Compiler
"""

import os
import sys
import json
import base64
from datetime import datetime
from pathlib import Path
import argparse


def image_to_base64(image_path):
    try:
        with open(image_path, 'rb') as f:
            return base64.b64encode(f.read()).decode('utf-8')
    except Exception as e:
        print(f"Warning: Could not convert image {image_path}: {e}")
        return None


def format_number(num):
    if num is None:
        return "0"
    try:
        return f"{int(num):,}"
    except Exception:
        return str(num)


def generate_html_report(workdir, genome, report_title, report_subtitle,
                         company_name, project_id, scan_motif_pvalue,
                         stats_file_override=None, pair_names=None):
    workdir = Path(workdir)
    reportdir = workdir / "report"
    datadir = reportdir / "data"
    figuredir = reportdir / "figures"

    stats_file = Path(stats_file_override) if stats_file_override else datadir / "statistics_all.json"
    if stats_file.exists():
        with open(stats_file, 'r', encoding='utf-8') as f:
            stats = json.load(f)
    else:
        stats = {}

    report_date = datetime.now().strftime("%Y-%m-%d")
    current_year = datetime.now().strftime("%Y")
    sample_count = len(stats.get("samples", {}))
    pairs = stats.get("pairs", [])
    pair_names = pair_names or [p["pair_name"] for p in pairs]

    # Build sample summary table rows
    table_rows = ""
    for sample_name, sample_data in stats.get("samples", {}).items():
        sample_type = sample_data.get("type", "Unknown")
        trim_data = sample_data.get("trimming", {})
        align_data = sample_data.get("alignment", {})
        flagstat_data = sample_data.get("flagstat", {})
        total_reads = trim_data.get("total_reads", 0)
        trimmed_reads = align_data.get("input_reads", trim_data.get("trimmed_reads", total_reads))
        mapped_reads = align_data.get("uniquely_mapped", flagstat_data.get("mapped_reads", 0))
        mapping_rate = align_data.get("mapped_rate", flagstat_data.get("mapping_rate", 0))

        table_rows += f"""
                <tr>
                    <td class="sample-cell">{sample_name}</td>
                    <td class="read-cell">{sample_type}</td>
                    <td>{format_number(total_reads)}</td>
                    <td>{format_number(trimmed_reads)}</td>
                    <td>{format_number(mapped_reads)}</td>
                    <td>{mapping_rate:.2f}%</td>
                </tr>"""

    # Unified figures
    figures = {}
    for fig_name in ["all_reads_statistics", "all_mapping_rate"]:
        fig_path = figuredir / f"{fig_name}.png"
        if fig_path.exists():
            figures[fig_name] = image_to_base64(fig_path)

    total_peaks_all = sum(p.get("peak_count_exomepeak2", 0) for p in pairs)

    # Build per-pair peak/motif sections
    pair_sections = ""
    for p in pairs:
        pair_name = p.get("pair_name", "")
        ip_sample = p.get("ip_sample", "")
        input_sample = p.get("input_sample", "")
        peak_count_exome = p.get("peak_count_exomepeak2", 0)
        peak_count_bed6 = p.get("peak_count_bed6", 0)
        m6a_motif_peaks = p.get("m6a_motif_peaks", 0)
        m6a_motif_ratio = p.get("m6a_motif_ratio", 0.0)

        peak_counts_img = image_to_base64(figuredir / f"{pair_name}_peak_counts.png")
        motif_pie_img = image_to_base64(figuredir / f"{pair_name}_m6a_motif_pie.png")

        pair_sections += f"""
                <div class="info-card mt-4">
                    <h5><i class="fas fa-layer-group me-2"></i>Pair: {pair_name}</h5>
                    <table class="table param-table mt-3">
                        <tr><td>IP样本</td><td>{ip_sample}</td></tr>
                        <tr><td>Input样本</td><td>{input_sample if input_sample else 'N/A'}</td></tr>
                        <tr><td>exomePeak2原始Peaks</td><td>{format_number(peak_count_exome)}</td></tr>
                        <tr><td>BED6格式转换后Peaks</td><td>{format_number(peak_count_bed6)}</td></tr>
                    </table>
                </div>

                <div class="row">
                    <div class="col-md-12">
                        <div class="figure-container">
                            <h6><i class="fas fa-mountain me-2"></i>Peak Counts ({pair_name})</h6>
                            {f'<img src="data:image/png;base64,{peak_counts_img}" alt="Peak Counts">' if peak_counts_img else '<p class="text-muted">图表生成中...</p>'}
                            <p class="figure-caption">Peak calling结果统计 ({pair_name})</p>
                        </div>
                    </div>
                </div>

                <div class="info-card">
                    <h5><i class="fas fa-chart-pie me-2"></i>ac4C Motif Peaks比例 ({pair_name})</h5>
                    <table class="table param-table mt-3">
                        <tr><td>总Peaks数 (BED6格式)</td><td>{format_number(peak_count_bed6)}</td></tr>
                        <tr><td>含ac4C motif的Peaks</td><td>{format_number(m6a_motif_peaks)}</td></tr>
                        <tr><td>ac4C motif占比</td><td>{m6a_motif_ratio:.2f}%</td></tr>
                    </table>
                </div>

                <div class="row">
                    <div class="col-md-12">
                        <div class="figure-container">
                            <h6><i class="fas fa-chart-pie me-2"></i>ac4C Motif Distribution ({pair_name})</h6>
                            {f'<img src="data:image/png;base64,{motif_pie_img}" alt="ac4C Motif Distribution">' if motif_pie_img else '<p class="text-muted">图表生成中...</p>'}
                            <p class="figure-caption">ac4C motif在Peaks中的分布比例 ({pair_name})</p>
                        </div>
                    </div>
                </div>
        """

    html_content = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{report_title} - {project_id}</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
    <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css" rel="stylesheet">
    <style>
        :root {{
            --primary-color: #2c5282;
            --secondary-color: #4a5568;
            --accent-color: #3182ce;
            --light-bg: #f7fafc;
            --border-color: #e2e8f0;
        }}
        body {{
            font-family: 'Segoe UI', -apple-system, BlinkMacSystemFont, sans-serif;
            line-height: 1.6;
            color: #2d3748;
        }}
        .report-header {{
            background: linear-gradient(135deg, var(--primary-color) 0%, var(--accent-color) 100%);
            color: white;
            padding: 3rem 0;
            margin-bottom: 2rem;
        }}
        .report-header h1 {{
            font-size: 2.5rem;
            font-weight: 700;
            margin-bottom: 0.5rem;
        }}
        .report-header .subtitle {{
            font-size: 1.2rem;
            opacity: 0.9;
        }}
        .company-badge {{
            background: rgba(255,255,255,0.15);
            padding: 0.5rem 1rem;
            border-radius: 20px;
            font-size: 0.9rem;
            display: inline-block;
            margin-top: 1rem;
        }}
        .three-line-table {{
            width: 100%;
            border-collapse: collapse;
            margin: 1.5rem 0;
            font-size: 0.95rem;
        }}
        .three-line-table thead {{
            border-top: 2px solid #2d3748;
            border-bottom: 1px solid #2d3748;
        }}
        .three-line-table thead th {{
            padding: 12px 8px;
            font-weight: 600;
            text-align: center;
            background: var(--light-bg);
        }}
        .three-line-table tbody tr {{
            border-bottom: 1px solid var(--border-color);
        }}
        .three-line-table tbody tr:last-child {{
            border-bottom: 2px solid #2d3748;
        }}
        .three-line-table tbody td {{
            padding: 10px 8px;
            text-align: center;
        }}
        .three-line-table tbody tr:hover {{
            background: #f7fafc;
        }}
        .three-line-table .sample-cell {{
            font-weight: 600;
            color: var(--primary-color);
            vertical-align: middle;
        }}
        .three-line-table .read-cell {{
            font-weight: 500;
            color: var(--secondary-color);
        }}
        .section {{
            margin-bottom: 3rem;
        }}
        .section-header {{
            border-left: 4px solid var(--accent-color);
            padding-left: 1rem;
            margin-bottom: 1.5rem;
        }}
        .section-header h2 {{
            font-size: 1.75rem;
            font-weight: 600;
            color: var(--primary-color);
            margin-bottom: 0.25rem;
        }}
        .section-header .section-desc {{
            color: var(--secondary-color);
            font-size: 0.95rem;
        }}
        .info-card {{
            background: var(--light-bg);
            border: 1px solid var(--border-color);
            border-radius: 8px;
            padding: 1.5rem;
            margin-bottom: 1.5rem;
        }}
        .info-card h5 {{
            color: var(--primary-color);
            font-weight: 600;
            margin-bottom: 1rem;
        }}
        .stat-card {{
            background: white;
            border: 1px solid var(--border-color);
            border-radius: 8px;
            padding: 1.5rem;
            text-align: center;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05);
            transition: transform 0.2s, box-shadow 0.2s;
        }}
        .stat-card:hover {{
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
        }}
        .stat-card .stat-icon {{
            font-size: 2rem;
            color: var(--accent-color);
            margin-bottom: 0.5rem;
        }}
        .stat-card .stat-value {{
            font-size: 2rem;
            font-weight: 700;
            color: var(--primary-color);
        }}
        .stat-card .stat-label {{
            color: var(--secondary-color);
            font-size: 0.9rem;
        }}
        .figure-container {{
            background: white;
            border: 1px solid var(--border-color);
            border-radius: 8px;
            padding: 1rem;
            margin-bottom: 1.5rem;
            text-align: center;
        }}
        .figure-container img {{
            max-width: 100%;
            height: auto;
        }}
        .figure-caption {{
            color: var(--secondary-color);
            font-size: 0.9rem;
            margin-top: 0.5rem;
            font-style: italic;
        }}
        .report-footer {{
            background: var(--light-bg);
            border-top: 1px solid var(--border-color);
            padding: 2rem 0;
            margin-top: 3rem;
            text-align: center;
            color: var(--secondary-color);
            font-size: 0.9rem;
        }}
        .toc {{
            position: sticky;
            top: 20px;
        }}
        .toc .list-group-item {{
            border: none;
            padding: 0.5rem 1rem;
            color: var(--secondary-color);
        }}
        .toc .list-group-item:hover {{
            background: var(--light-bg);
            color: var(--accent-color);
        }}
        .toc .list-group-item i {{
            margin-right: 0.5rem;
            width: 20px;
        }}
        .param-table {{
            font-size: 0.9rem;
        }}
        .param-table td:first-child {{
            font-weight: 600;
            color: var(--primary-color);
            width: 30%;
        }}
        .motif-link-card {{
            background: white;
            border: 1px solid var(--border-color);
            border-radius: 8px;
            padding: 1.5rem;
            margin-bottom: 1rem;
            transition: all 0.3s;
        }}
        .motif-link-card:hover {{
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            border-color: var(--accent-color);
        }}
        .motif-link-card a {{
            color: var(--accent-color);
            text-decoration: none;
            font-weight: 500;
        }}
        .motif-link-card a:hover {{
            text-decoration: underline;
        }}
    </style>
</head>
<body>
    <header class="report-header">
        <div class="container">
            <div class="row align-items-center">
                <div class="col-md-8">
                    <h1><i class="fas fa-dna me-3"></i>{report_title}</h1>
                    <p class="subtitle">{report_subtitle}</p>
                    <span class="company-badge">
                        <i class="fas fa-building me-2"></i>{company_name}
                    </span>
                </div>
                <div class="col-md-4 text-md-end">
                    <div class="mt-3">
                        <p class="mb-1"><strong>Project ID:</strong> {project_id}</p>
                        <p class="mb-1"><strong>Genome:</strong> {genome}</p>
                        <p class="mb-0"><strong>Report Date:</strong> {report_date}</p>
                    </div>
                </div>
            </div>
        </div>
    </header>

    <div class="container">
        <div class="row">
            <div class="col-md-3 d-none d-md-block">
                <div class="toc">
                    <div class="card">
                        <div class="card-header bg-primary text-white">
                            <i class="fas fa-list me-2"></i>目录
                        </div>
                        <div class="list-group list-group-flush">
                            <a href="#overview" class="list-group-item">
                                <i class="fas fa-info-circle"></i>项目概述
                            </a>
                            <a href="#summary" class="list-group-item">
                                <i class="fas fa-table"></i>样本统计汇总
                            </a>
                            <a href="#peaks" class="list-group-item">
                                <i class="fas fa-mountain"></i>Peaks统计
                            </a>
                            <a href="#motif" class="list-group-item">
                                <i class="fas fa-search"></i>Motif分析
                            </a>
                            <a href="#methods" class="list-group-item">
                                <i class="fas fa-cogs"></i>方法学说明
                            </a>
                        </div>
                    </div>
                </div>
            </div>
            
            <div class="col-md-9">
                <section id="overview" class="section">
                    <div class="section-header">
                        <h2><i class="fas fa-info-circle me-2"></i>项目概述</h2>
                        <p class="section-desc">Project Overview</p>
                    </div>
                    
                    <div class="info-card">
                        <h5><i class="fas fa-project-diagram me-2"></i>A. 项目信息</h5>
                        <table class="table param-table">
                            <tr><td>项目名称</td><td>{report_title}</td></tr>
                            <tr><td>项目编号</td><td>{project_id}</td></tr>
                            <tr><td>参考基因组</td><td>{genome}</td></tr>
                            <tr><td>样本数量</td><td>{sample_count} 个样本</td></tr>
                            <tr><td>分析Pair数</td><td>{len(pairs)} 个pair</td></tr>
                            <tr><td>分析平台</td><td>pi-omicskit Bioinformatics Pipeline</td></tr>
                            <tr><td>分析日期</td><td>{report_date}</td></tr>
                        </table>
                    </div>
                    
                    <div class="info-card">
                        <h5><i class="fas fa-microscope me-2"></i>B. acRIP-seq 技术简介</h5>
                        <p>
                            acRIP-seq (acetylated RNA Immunoprecipitation Sequencing) 是一种用于检测RNA分子上表观遗传修饰的高通量测序技术。
                            该技术通过特异性识别RNA修饰（如ac4C乙酰化）的抗体进行免疫沉淀，富集带有修饰的RNA片段，然后进行高通量测序分析。
                            acRIP-seq能够全基因组范围内鉴定RNA修饰位点，对于研究RNA表观遗传学调控机制具有重要意义。
                        </p>
                    </div>
                    
                    <div class="info-card">
                        <h5><i class="fas fa-bullseye me-2"></i>C. 分析目标</h5>
                        <p>本项目的主要分析目标包括：</p>
                        <ol>
                            <li>对acRIP-seq原始数据进行质量控制和预处理</li>
                            <li>将测序reads比对到参考基因组</li>
                            <li>使用exomePeak2进行Peak calling，鉴定RNA修饰富集区域</li>
                            <li>进行Motif分析，识别可能的RNA修饰结合序列模式</li>
                            <li>对已知ac4C motif进行扫描和定量分析</li>
                        </ol>
                    </div>
                    
                    <div class="row">
                        <div class="col-md-3">
                            <div class="stat-card">
                                <div class="stat-icon"><i class="fas fa-vial"></i></div>
                                <div class="stat-value">{sample_count}</div>
                                <div class="stat-label">样本数量</div>
                            </div>
                        </div>
                        <div class="col-md-3">
                            <div class="stat-card">
                                <div class="stat-icon"><i class="fas fa-mountain"></i></div>
                                <div class="stat-value">{format_number(total_peaks_all)}</div>
                                <div class="stat-label">总Peaks</div>
                            </div>
                        </div>
                        <div class="col-md-3">
                            <div class="stat-card">
                                <div class="stat-icon"><i class="fas fa-search"></i></div>
                                <div class="stat-value">STREME</div>
                                <div class="stat-label">Motif分析</div>
                            </div>
                        </div>
                        <div class="col-md-3">
                            <div class="stat-card">
                                <div class="stat-icon"><i class="fas fa-dna"></i></div>
                                <div class="stat-value">ac4C</div>
                                <div class="stat-label">Motif扫描</div>
                            </div>
                        </div>
                    </div>
                </section>
                
                <section id="summary" class="section">
                    <div class="section-header">
                        <h2><i class="fas fa-table me-2"></i>样本统计汇总</h2>
                        <p class="section-desc">Sample Statistics Summary</p>
                    </div>
                    
                    <div class="info-card">
                        <h5><i class="fas fa-chart-bar me-2"></i>A. 数据质量汇总表</h5>
                        <p class="text-muted mb-3">下表展示了各样本的数据处理和质量控制统计信息。</p>
                        
                        <div class="table-responsive">
                            <table class="three-line-table">
                                <thead>
                                    <tr>
                                        <th>Sample</th>
                                        <th>Type</th>
                                        <th>Total Reads</th>
                                        <th>Trimmed Reads</th>
                                        <th>Mapped Reads</th>
                                        <th>Mapping Rate</th>
                                    </tr>
                                </thead>
                                <tbody>
                                    {table_rows}
                                </tbody>
                            </table>
                        </div>
                        
                        <p class="text-muted small mt-2">
                            <strong>注：</strong> Total Reads为原始测序reads数；Trimmed Reads为质控过滤后的reads数；
                            Mapped Reads为成功比对到参考基因组的reads数；Mapping Rate为比对成功率。
                        </p>
                    </div>
                    
                    <div class="row">
                        <div class="col-md-12">
                            <div class="figure-container">
                                <h6><i class="fas fa-chart-bar me-2"></i>Reads Statistics (All Samples)</h6>
                                {f'<img src="data:image/png;base64,{figures.get("all_reads_statistics", "")}" alt="Reads Statistics">' if figures.get("all_reads_statistics") else '<p class="text-muted">图表生成中...</p>'}
                                <p class="figure-caption">所有样本reads数量统计（单位：百万）</p>
                            </div>
                        </div>
                    </div>
                    
                    <div class="row">
                        <div class="col-md-12">
                            <div class="figure-container">
                                <h6><i class="fas fa-map-marker-alt me-2"></i>Mapping Rate (All Samples)</h6>
                                {f'<img src="data:image/png;base64,{figures.get("all_mapping_rate", "")}" alt="Mapping Rate">' if figures.get("all_mapping_rate") else '<p class="text-muted">图表生成中...</p>'}
                                <p class="figure-caption">所有样本比对成功率统计</p>
                            </div>
                        </div>
                    </div>
                </section>
                
                <section id="peaks" class="section">
                    <div class="section-header">
                        <h2><i class="fas fa-mountain me-2"></i>Peaks统计</h2>
                        <p class="section-desc">Peak Calling Results</p>
                    </div>
                    
                    <div class="info-card">
                        <h5><i class="fas fa-chart-pie me-2"></i>Peak calling数量统计（所有Pairs汇总）</h5>
                        <p>
                            使用exomePeak2软件对各IP样本进行Peak calling分析，鉴定RNA修饰富集区域。
                            以下是各Pair的Peak calling统计结果：
                        </p>
                    </div>
                    
                    {pair_sections}
                </section>
                
                <section id="motif" class="section">
                    <div class="section-header">
                        <h2><i class="fas fa-search me-2"></i>Motif分析</h2>
                        <p class="section-desc">Motif Analysis Results</p>
                    </div>
                    
                    <div class="info-card">
                        <h5><i class="fas fa-dna me-2"></i>A. 无偏Motif寻找 (STREME)</h5>
                        <p>使用STREME软件对所有Peaks进行无偏Motif寻找，发现可能的RNA修饰相关序列模式。</p>
                        <div class="motif-link-card mt-3">
                            <i class="fas fa-external-link-alt me-2 text-primary"></i>
                            <a href="../motif_analysis/" target="_blank">
                                查看各Pair的STREME分析结果目录
                            </a>
                            <p class="text-muted small mt-2 mb-0">
                                <i class="fas fa-folder-open me-1"></i>文件位置: <code>{workdir}/motif_analysis/</code>
                            </p>
                        </div>
                    </div>
                    
                    <div class="info-card">
                        <h5><i class="fas fa-crosshairs me-2"></i>B. ac4C Motif扫描 (FIMO)</h5>
                        <p>
                            使用FIMO软件对已知ac4C motif (C-rich)在所有Peaks中进行扫描，
                            p值阈值设置为 {scan_motif_pvalue}。
                        </p>
                        <div class="motif-link-card mt-3">
                            <i class="fas fa-external-link-alt me-2 text-primary"></i>
                            <a href="../scan_motif/" target="_blank">
                                查看各Pair的FIMO扫描结果目录
                            </a>
                            <p class="text-muted small mt-2 mb-0">
                                <i class="fas fa-folder-open me-1"></i>文件位置: <code>{workdir}/scan_motif/</code>
                            </p>
                        </div>
                    </div>
                </section>
                
                <section id="methods" class="section">
                    <div class="section-header">
                        <h2><i class="fas fa-cogs me-2"></i>方法学说明</h2>
                        <p class="section-desc">Methodology</p>
                    </div>
                    
                    <div class="info-card">
                        <h5><i class="fas fa-stream me-2"></i>A. 分析流程</h5>
                        <ol>
                            <li><strong>质量控制与预处理：</strong>使用 FastQC 进行原始数据质控，使用 Trim Galore 进行接头序列切除和质量过滤</li>
                            <li><strong>序列比对：</strong>使用 STAR 将 reads 比对到 {genome} 参考基因组</li>
                            <li><strong>去重处理：</strong>使用 samtools 进行PCR重复序列去除</li>
                            <li><strong>Peak calling：</strong>使用 exomePeak2 进行 RNA 修饰富集区域鉴定</li>
                            <li><strong>Motif分析：</strong>使用 STREME 进行无偏 Motif 发现</li>
                            <li><strong>Motif扫描：</strong>使用 FIMO 对已知 ac4C motif 进行扫描</li>
                        </ol>
                    </div>
                    
                    <div class="info-card">
                        <h5><i class="fas fa-sliders-h me-2"></i>B. 主要软件版本</h5>
                        <table class="table param-table">
                            <tr><td>FastQC</td><td>v0.12.1</td></tr>
                            <tr><td>Trim Galore</td><td>v0.6.10</td></tr>
                            <tr><td>STAR</td><td>v2.7.11b</td></tr>
                            <tr><td>samtools</td><td>v1.18</td></tr>
                            <tr><td>exomePeak2</td><td>v1.10</td></tr>
                            <tr><td>bedtools</td><td>v2.31.0</td></tr>
                            <tr><td>MEME Suite (STREME)</td><td>v5.5.5</td></tr>
                        </table>
                    </div>
                    
                    <div class="info-card">
                        <h5><i class="fas fa-book me-2"></i>C. 方法学补充</h5>
                        <p>
                            <strong>Peak calling参数：</strong>使用 exomePeak2 进行peak calling，
                            基于负二项分布模型进行统计检验，鉴定IP样本相对于Input样本的显著富集区域。
                        </p>
                        <p class="mb-0">
                            <strong>Motif分析参数：</strong>STREME用于发现de novo motif，
                            设置motif宽度为5bp，显著性阈值0.1；FIMO用于扫描已知motif，
                            p值阈值设置为 {scan_motif_pvalue}。
                        </p>
                    </div>
                    
                    <div class="info-card">
                        <h5><i class="fas fa-forward me-2"></i>D. 后续工作与交互</h5>
                        <p>本报告提供了acRIP-seq分析的主要结果。如需进一步分析，建议考虑以下方向：</p>
                        <ul>
                            <li>Peak注释：将Peaks定位到基因区域（5'UTR, CDS, 3'UTR, intron等）</li>
                            <li>差异修饰分析：如有多个样本组，可进行差异RNA修饰分析</li>
                            <li>Motif富集分析：进一步分析特定基因区域的motif分布</li>
                            <li>功能富集分析：对含有Peaks的基因进行GO/KEGG富集分析</li>
                        </ul>
                    </div>
                </section>
            </div>
        </div>
    </div>
    
    <footer class="report-footer">
        <div class="container">
            <p class="mb-1">
                <i class="fas fa-flask me-2"></i>
                <strong>pi-omicskit Bioinformatics</strong>
            </p>
            <p class="mb-1">Professional Bioinformatics Service Provider</p>
            <p class="text-muted small mb-0">
                Report generated on {report_date} | Project: {project_id} | 
                <i class="fas fa-envelope me-1"></i>support@pi-omicskit.com
            </p>
        </div>
    </footer>

    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js"></script>
</body>
</html>"""

    return html_content


def main():
    parser = argparse.ArgumentParser(description='Compile Unified acRIP-seq HTML Report')
    parser.add_argument('-w', '--workdir', default='result', help='Working directory')
    parser.add_argument('-g', '--genome', default='hg38', help='Genome version')
    parser.add_argument('-t', '--title', default='acRIP-seq Analysis Report', help='Report title')
    parser.add_argument('-st', '--subtitle', default='ac4C RNA Immunoprecipitation sequencing Analysis Report', help='Report subtitle')
    parser.add_argument('-c', '--company', default='pi-omicskit Bioinformatics', help='Company name')
    parser.add_argument('-p', '--project', default='acRIP2024', help='Project ID')
    parser.add_argument('--scan-motif-pvalue', type=float, default=0.01, help='Scan motif p-value threshold')
    parser.add_argument('--stats-json', required=True, help='Unified statistics JSON file path')
    parser.add_argument('--pair-names', default='', help='Comma-separated pair names')
    parser.add_argument('-o', '--output', required=True, help='Output HTML file')
    args = parser.parse_args()

    workdir = Path(args.workdir)
    genome = args.genome
    report_title = args.title
    report_subtitle = args.subtitle
    company_name = args.company
    project_id = args.project
    scan_motif_pvalue = args.scan_motif_pvalue
    stats_json = args.stats_json
    output_html = Path(args.output)
    pair_names = [x.strip() for x in args.pair_names.split(',') if x.strip()] if args.pair_names else None

    print("Generating unified HTML report...")
    html_content = generate_html_report(
        workdir, genome, report_title, report_subtitle,
        company_name, project_id, scan_motif_pvalue,
        stats_json, pair_names
    )

    output_html.parent.mkdir(parents=True, exist_ok=True)
    with open(output_html, 'w', encoding='utf-8') as f:
        f.write(html_content)

    print(f"Unified HTML report saved to: {output_html}")


if __name__ == "__main__":
    main()
