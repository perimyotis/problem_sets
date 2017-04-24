library(data.table)
counts <- fread ("~/Desktop/sequencing_class/pset7/7.4.txt")
counts
counts.log2FC.twofold = counts[abs(counts$log2fold) >= 2.0,]