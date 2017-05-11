library(data.table)
counts <- fread ("~/Desktop/sequencing_class/pset7/7.4.txt")
counts

#7.4 D
countslog <- counts[,.(Name, (log2(BY_counts/RM_counts)))]
countslog

#7.4 E
countslog2 <- counts[BY_counts >0 & RM_counts > 0, .(Name, (log2(BY_counts/RM_counts)))]
countslog2

#7.4 F

pseudocounts <- counts[,.(Name, BY_counts +1, RM_counts +1)]
pseudocounts

#7.4 G

pseudocountslog <- pseudocounts [,. (Name, (log2(V2/V3)))]
pseudocountslog