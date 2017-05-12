library(DESeq2)
library(tximport)
library(data.table)
library(GenomeInfoDb)
library(rhdf5)

RNA_Files = c("~/Desktop/sequencing_class/final/ScermRNAOut0/abundance.tsv","~/Desktop/sequencing_class/final/ScermRNAOut1/abundance.tsv","~/Desktop/sequencing_class/final/SparmRNAOut0/abundance.tsv","~/Desktop/sequencing_class/final/SparmRNAOut1/abundance.tsv")
names(RNA_Files) = c("Scer1","Scer2","Spar1","Spar2")

RNA_txdat = tximport(RNA_Files, type = "kallisto", txOut = TRUE)

RNA_coldata = data.frame(condition = c("Scer","Scer","Spar","Spar"))
rownames(RNA_coldata) = names(RNA_Files)
RNA_coldata

RNA_dds = DESeqDataSetFromTximport(RNA_txdat, colData = RNA_coldata, design =~ condition)
RNA_dds = DESeq(RNA_dds)
RNA_res = results(RNA_dds)

RNA_res_FDR10 = subset(RNA_res, padj <0.1)
write.csv(as.data.frame(RNA_res_FDR10[order(RNA_res_FDR10$log2FoldChange),]), file = "~/Desktop/sequencing_class/final/mRNA_diffexp.csv", col.name=FALSE, row.name=FALSE)

RNA_res_table = as.data.table(as.data.frame(RNA_res), keep.rownames="id")
RNA_res_raw = RNA_res_table[,.(id,log2FoldChange)]
write.table(RNA_res_raw, file ="~/Desktop/sequencing_class/final/mRNA_raw.txt")

RNA_res_sig = RNA_res_table[padj<0.1,.(id,log2FoldChange)]
write.table(RNA_rew_sig, file = "~/Desktop/sequencing_class/final/mRNA_diffexpsig.txt", sol.name = FALSE, row.name = FALSE)

plotMA(RNA_res)
plotDispEsts(RNA_dds)

RNA_pvalue_sb = sum(RNA_res$pvalue <0.1, na.rm = TRUE)
RNA_pvalue_sb

RNA_false_discovery = sum(RNA_res$padj <0.1, na.rm = TRUE)
RNA_false_discovery 

Ribo_files = c("~/Desktop/sequencing_class/final/ScerRiboOut0/abundance.tsv","~/Desktop/sequencing_class/final/ScerRiboOut1/abundance.tsv","~/Desktop/sequencing_class/final/SparRiboOut0/abundance.tsv","~/Desktop/sequencing_class/final/SparRiboOut1/abundance.tsv")
names(Ribo_files) = c("Scer1","Scer2","Spar1","Spar2")

Ribo_txdat = tximport(Ribo_files, type="kallisto", txOut = TRUE)

Ribo_coldata = data.frame(condition = c("Scer","Scer","Spar","Spar"))
rownames(Ribo_coldata) = names(Ribo_files)
Ribo_coldata

Ribo_dds = DESeqDataSetFromTximport(Ribo_tsdat, colData = Ribo_coldata, design=~ conditon)
Ribo_dds = DESeq(Ribo_dds)

Ribo_res = results(Ribo_dds)
Ribo_res_FDR10 = subset(Ribo_res, padj <0.1)
write.csv(as.data.frame(Ribo_res_FDR10[order(Ribo_res_FDR10$log2FoldChange),]), file = "~/Desktop/sequencing_class/final/Ribo_diffoccusig.txt")

plotMA(Ribo_res)
plotDisEsts(Ribo_dds)

Ribo_pvalue_sb = sum(Ribo_res$pvalue <0.1, na.rm = TRUE)
Ribo_pvalue_sb

Ribo_false_discovery = sum(Ribo_res$padj <0.1, na.rm = TRUE)
Ribo_false_discovery

TE_files = c("~/Desktop/sequencing_class/final/ScermRNAOut0/abundance.tsv","~/Desktop/sequencing_class/final/ScermRNAOut1/abundance.tsv","~/Desktop/sequencing_class/final/SparmRNAOut0/abundance.tsv","~/Desktop/sequencing_class/final/SparmRNAOut1/abundance.tsv","~/Desktop/sequencing_class/final/ScerRiboOut0/abundance.tsv","~/Desktop/sequencing_class/final/ScerRiboOut1/abundance.tsv","~/Desktop/sequencing_class/final/SparRiboOut0/abundance.tsv","~/Desktop/sequencing_class/final/SparRiboOut1/abundance.tsv")
names(TE_files) = c("ScerRNA1","ScerRNA2","SparRNA1","SparRNA2","ScerRibo1","ScerRibo2","SparRibo1","SparRibo2")

TE_txdat = tximport(TE_files, type = "kallisto", txOut = TRUE)

TE_coldata = data.frame(condition = c("Scer","Scer","Spar","Spar","Scer","Scer","Spar","Spar"), assay = c("RNA","RNA","RNA","RNA","Ribo","Ribo","Ribo","Ribo"))
rownames(TE_coldata) = names(TE_files)
TE_coldata

TE_dds = DESeqDataSetFromTximport(TE_txdat, colData = TE_coldata, design =~ assay + condition + assay:condition)
TE_res = results(TE_dds)
TE_res_FDR10 = subset(TE_res, padj <0.1)
write.csv(as.data.frame(TE_res_FDR10[order(TE_res_FDR10$log2FoldChange),]), file = "~/Desktop/sequencing_class/final/TE_diffsig.csv")

TE_res_table = as.data.table(as.data.frame(TE_res), keep.rownames="id")
TE_res_raw = TE_res_table[,.(id,log2FoldChange)]
write.table(TE_res_raw, file = "~/Desktop/sequencing_class/final/TE_raw.txt", col.name = FALSE)

TE_res_sig = TE_res_table[padj <0.1,.(id,log2FoldChange)]
write.table(TE_res_sig, file = "~/Desktop/sequencing_class/final/TE_sig.txt", col.name = FALSE)

plotMA(TE_res)
plotDispEsts(TE_dds)

TE_pvalue_sb = sum(TE_res$pvalue <0.1, na.rm = TRUE)
TE_pvalue_sb

TE_false_discovery = sum(TE_res$padj <0.1, na.rm = TRUE)
TE_false_dicsovery

plot(RNA_res$log2FoldChange, TE_res$log2FoldChange, main = "Log fold change in mRNA abundance vs. translation efficiency", xlab = "LFC in mRNA abundance", ylab = "LFC in translation efficiency")

RNA_raw = read.table("~/Desktop/sequencing_class/final/mRNA_raw.txt")
TE_raw = read.table("~/Desktop/sequencing_class/final/TE_raw.txt")
table = as.data.table(merge(RNA_raw, TE_raw, by "V1"))

RNAup_TEdown <- length(which(table$V2.x>0&table$V2.y<0))
RNAup_TEdown

RNAdown_TEup <- length(which(table$V2.x<0&table$V2.y>0))
RNAdown_TEup

RNAup_TEup <- length(which(table$V2.x>0&table$V2.y>0))
RNAup_TEup

RNAdown_TEdown <- length(which(table$V2.x<0&table$V2.y<0))
RNAdown_TEdown

RNA_sig = read.table("~/Desktop/sequencing_class/final/mRNA_diffexpsig.txt")
TE_sig = read.table("~/Desktop/sequencing_class/final/TE_sig.txt")
table_sig = as.data.table(merge(RNA_sig, TE_sig, by="V1", all = FALSE))
table_sig_num <- transform(table_sig, V2.x = as.numeric(V2.x), V2.y = as.numeric(V2.y))
sapply(table_sig_num, mode)

RNAup_TEdown_sig=length(which(table_sig_num$V2.x>0&table_sig_num$V2.y<0))
RNAup_TEdown_sig

RNAdown_TEup_sig=length(which(table_sig_num$V2.x<0&table_sig_num$V2.y>0))
RNAdown_TEup_sig

RNAup_TEup_sig=length(which(table_sig_num$V2.x>0&table_sig_num$V2.y>0))
RNAup_TEup_sig

RNAdown_TEdown_sig=length(which(table_sig_num$V2.x<0&table_sig_num$V2.y<0))
RNAdown_TEdown_sig