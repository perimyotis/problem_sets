import pysam
import matplotlib.pyplot as mpl
samfile = pysam.AlignmentFile ("mappedreads.sorted.bam", "rb")
mapquallist= []
for read in samfile.fetch():
	mapqual = read.mapping_quality
	mapquallist.append(mapqual)
	mapstrand = read.is_reverse
mpl.hist(mapquallist)
mpl.savefig("mapqualhist.pdf")

plusnum = 0
for strand in strands:
	if strand ==False:
		plusnum +=1
plusfrac = plusnum/len(strands)
print "%i of reads are on the plus strand"%(plusfrac)

