library(tximport)
library(DESeq2)
library(data.table)

files = c("~/Desktop/sequencing_class/pset8/WT_1/abundance.tsv","~/Desktop/sequencing_class/pset8/WT_2/abundance.tsv","~/Desktop/sequencing_class/pset8/WT_3/abundance.tsv","~/Desktop/sequencing_class/pset8/WT_4/abundance.tsv","~/Desktop/sequencing_class/pset8/WT_5/abundance.tsv","~/Desktop/sequencing_class/pset8/SNF2_1/abundance.tsv","~/Desktop/sequencing_class/pset8/SNF2_2/abundance.tsv","~/Desktop/sequencing_class/pset8/SNF2_3/abundance.tsv","~/Desktop/sequencing_class/pset8/SNF2_4/abundance.tsv","~/Desktop/sequencing_class/pset8/SNF2_5/abundance.tsv")
names(files) =c("WT1","WT2","WT3","WT4","WT5","SNF2_1","SNF2_2","SNF2_3","SNF2_4","SNF2_5")

txdat = tximport(files, type = "kallisto", txOut=TRUE)

coldata = data.frame(condition=c("WT","WT","WT","WT","WT","SNF2","SNF2","SNF2","SNF2","SNF2"))
rownames(coldata) = names(files)
coldata

dds = DESeqDataSetFromTximport(txdat, colData = coldata, design=~ condition)
dds = DESeq(dds)
res = results(dds)

plotMA(res)
plotDispEsts(dds)

res$padj
res$pvalue

sum(res$pvalue < .05, na.rm = TRUE)
sum(res$padj < .05, na.rm = TRUE)
