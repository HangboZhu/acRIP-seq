suppressPackageStartupMessages(library(GenomicRanges))
r <- readRDS("/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/result_full/exomePeak2_output_GSE288115.HepAD38.WT/exomePeak2_result.rds")
g <- unlist(r, use.names=FALSE)
cat("peak 总数:", length(g), "\n")
cat("元数据列:", paste(colnames(mcols(g)), collapse=" | "), "\n")
mc <- as.data.frame(mcols(g))
pcol <- grep("pvalue|padj|p\\.value", colnames(mc), value=TRUE)[1]
fcol <- grep("log2FC|log2fc|fold|score", colnames(mc), value=TRUE)[1]
cat("排序用 P 列:", pcol, " FC 列:", fcol, "\n")
ord <- order(mc[[pcol]])[1:8]
tab <- data.frame(chr=as.character(seqnames(g))[ord], start=start(g)[ord], end=end(g)[ord],
                  width=width(g)[ord], mc[ord, c(fcol, pcol), drop=FALSE], check.names=FALSE)
cat("\n== top 8 peaks（按 P 值）==\n")
print(tab, row.names=FALSE, digits=3)
cat("\n== peaks.bed 第 4/5 列语义对照（bed 前两行）==\n")
p <- read.delim("/hpcdisk1/apod/shijiachen/yuliao/acRIPseq_pipeline/result_full/exomePeak2_output_GSE288115.HepAD38.WT/peaks.bed", header=FALSE, nrows=2)
cat("bed V4(峰ID) V5(score):", p$V4[1], p$V5[1], "|", p$V4[2], p$V5[2], "\n")
