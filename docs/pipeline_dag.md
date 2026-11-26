# acRIP-seq pipeline DAG

生成方式（可信来源，非手绘）：

```bash
ACRIP_CONFIG=Scripts/zhuhb_run/config.rat.yaml \
snakemake -s Scripts/pipeline.smk --directory .zhuhb_run/dryrun --nolock --dag \
  | dot -Tsvg > docs/dag.svg
```

规则级边关系如下（由 `snakemake --dag` 输出的 DOT 解析得到）：

```
collect_all_statistics  -> compile_unified_report, generate_unified_figures
collect_statistics      -> collect_all_statistics, generate_report_figures
compile_unified_report  -> all
fastqc                  -> all, multiqc
generate_report_figures -> compile_unified_report
generate_unified_figures-> compile_unified_report
mark_duplicates         -> all, collect_statistics, multiqc, peak_calling
merge_peaks             -> all
multiqc                 -> all, compile_unified_report
peak_calling            -> all, collect_statistics, merge_peaks
peak_groups_manifest    -> all
star_alignment          -> all, collect_statistics, mark_duplicates, multiqc
trim_galore_pe          -> all, collect_statistics, multiqc, star_alignment
trim_galore_se          -> all, collect_statistics, multiqc, star_alignment
```

## 三个扇出层级

| 层级 | 规则 | 粒度 | rat | fly | human | mouse |
|---|---|---|---|---|---|---|
| A 每样本 | fastqc | 1 / 样本 | 12 | 18 | 120 | 150 |
| A | trim_galore_se / pe | 1 / 样本 | 12 | 18 | 120 | 150 |
| A | star_alignment | 1 / 样本 | 12 | 18 | 120 | 150 |
| A | mark_duplicates | 1 / 样本 | 12 | 18 | 120 | 150 |
| — | multiqc | 1 / 物种 | 1 | 1 | 1 | 1 |
| **B 每 peak group** | **peak_calling** | **1 / (1 IP + 1 Input)** | **6** | **9** | **60** | **75** |
| B | collect_statistics | 1 / peak group | 6 | 9 | 60 | 75 |
| B | generate_report_figures | 1 / peak group | 6 | 9 | 60 | 75 |
| — | merge_peaks / peak_groups_manifest | 1 / 物种 | 2 | 2 | 2 | 2 |
| C 全局 | collect_all_statistics | 1 / 物种 | 1 | 1 | 1 | 1 |
| C | generate_unified_figures | 1 / 物种 | 1 | 1 | 1 | 1 |
| C | compile_unified_report | 1 / 物种 | 1 | 1 | 1 | 1 |
| | **合计 job** | | **72** | **105** | **666** | **831** |

## 完整结构（ASCII）

```
╔══════════════════════════════════════════════════════════════════════════════╗
║ Zone A · 每样本层   fan-out = S 个样本（human 120 / mouse 150 / rat 12 / fly 18）║
╚══════════════════════════════════════════════════════════════════════════════╝
  <sample>_R1.fastq.gz [ + R2 ]
       │
       ├─────────────► fastqc ─────────────► qc/<sample>_R1(_R2)_fastqc.html
       │                    │
       └─► trim_galore_se   │
           trim_galore_pe   │      （SE / PE 由 SAMPLE_LAYOUT 自动选择）
                │           │
                │           ├───────────────► trim/<sample>.{se,pe}.trim.log
                ▼           │
          star_alignment ───┼───────────────► bam/<sample>.STAR.log
                │           │                 bam/<sample>.Aligned.sortedByCoord.out.bam
                ▼           │
          mark_duplicates ──┼───────────────► bam/<sample>.DeDup.bam (+ .bai)
                            │                 bam/<sample>.flagstat.txt
                            │                 bam/<sample>.DeDup.log
                            │
                            ├──► multiqc ──► multiqc/multiqc_report.html
                            │
╔═══════════════════════════╪══════════════════════════════════════════════════╗
║ Zone B · 每 peak group 层  │  G = human 60 / mouse 75 / rat 6 / fly 9          ║
║   每个 group = 1 个 IP + 1 个 Input（peak_pairing: per_replicate）             ║
╚═══════════════════════════╪══════════════════════════════════════════════════╝
   <pair>.rep1.IP.DeDup.bam ┐
   <pair>.rep1.Input.DeDup.bam ┘
        ▼
   peak_calling ──► exomePeak2_output_<pair>_rep1/peaks.bed
        │           exomePeak2_output_<pair>_rep1/exomePeak2_result.rds
        │           （P<1e-10，GC 校正，exon 模式，1 IP vs 1 Input）
        │
        │      rep2 / rep3 …… 各自独立一条支线，互不依赖
        │
        ├──► [liftover_peaks]        可选，默认关 ──► peaks.lifted.bed
        ├──► [motif_analysis]        可选，默认关 ──► bed6_peak_<pair>_repN.bed
        │         └──► [scanmotif_analysis] 可选，默认关 ──► fimo_*.tsv
        ├──► merge_peaks ──► acRIP_all_peaks.bed（全部 group 的 cat+sort 并集）
        │
        ├──► peak_consensus ──► consensus/<pair>_consensus_peaks.bed
        │        （同一 pair 的 rep1+rep2+rep3 取交集：每个 rep 都必须支持）
        │    └─► consensus/<pair>_peak_support.txt
        │         └─► consensus/consensus_summary.tsv（全部 pair 汇总）
        │
        ▼
   collect_statistics ──► report/data/statistics_<pair>_repN.json
        ▼
   generate_report_figures ──► report/figures/<pair>_repN_*.png
        │                （reads_statistics / mapping_rate / peak_counts / m6a_motif_pie）
╔═══════╪══════════════════════════════════════════════════════════════════════╗
║ Zone C · 全局汇总层   fan-in = G 个 json + G 组图 + multiqc                    ║
╚═══════╪══════════════════════════════════════════════════════════════════════╝
        ▼
   collect_all_statistics ──► report/data/statistics_all.json
        ├──► generate_unified_figures ──► report/figures/all_{reads_statistics,mapping_rate}.png
        └──► compile_unified_report ──► report/acRIP_Analysis_Report_{CN,EN}.html

   旁支：peak_groups_manifest ──► result_*/peak_groups.tsv（peak_id ↔ IP/Input 对照）
```

## 并行特性

- **Zone A 完全并行**：样本之间零依赖，STAR 是主要耗时点。
- **Zone B 在样本层之上再展开**：`<pair>_rep1` 只要 rep1 的两个 BAM 就绪即可开跑，不必等整组。
  （旧版联合 calling 必须等一组 6 个样本全部 STAR 完才能开始。）
- **Zone C 是串行收口点**：全部 json + 图 + multiqc 齐备后才出总报告。

## 与上游 π-OmicsKit MeRIP-seq kit 的结构差异

| | 上游 kit | 本地 |
|---|---|---|
| pair 粒度 | 1 IP + 1 Input，无重复概念 | 1 IP + 1 Input（`per_replicate`），重复在 config 层组织 |
| peak 数 | 每条 config 条目 1 次 | 每条 config 条目展开为 N 次（N = 该组重复数） |
| SE / PE | 仅 SE | SE + PE 自动识别 |
| 多样本汇总 | 有（multi-pair 报告） | 有（i18n 中英文报告） |
| merge_peaks | 无 | 有（并集，非 consensus） |
| 跨重复 consensus | 无（上游无重复概念） | 有，严格取交集 `peak_consensus` |

## consensus 的三种模式

`consensus_min_support` 控制一条区域需要几个重复支持才保留：

| 取值 | 含义 | 3 重复时 | 风险 |
|---|---|---|---|
| `"all"`（默认） | 严格交集，每个重复都要支持 | ≥3 | 一个弱重复会把整组 peak 清空 |
| `"majority"` | 多数投票 | ≥2 | 允许单个重复不同 |
| 整数（如 `2`） | 指定支持数 | ≥2 | 同 majority |

`consensus_min_overlap`（默认 0.5）规定一个重复的 peak 要覆盖该区域至少多少比例才算支持。
单重复分组一律直通（support 全为 1），不会因为取交集而丢数据。
| liftOver | 无 | 有（默认关） |
| bamCoverage / bigWig | 有 | 无 |
| motif（STREME/FIMO） | 有，默认开（DRACH） | 有，默认关（DRACH 是 m6A motif） |
| 执行后端 | 本地 | `--executor slurm`，按物种独立工作目录 |
