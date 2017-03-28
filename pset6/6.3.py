import gffutils as gu
import pysam

db = gu.FeatureDB("yeast.db")

bamfile = pysam.AlignmentFile("set6readssort.bam", "rb")


Genecount = []
Genelen = []
mappedreads = 0.0
for mRNA in db.features_of_type("mRNA"):
	mappedreads +=1

	mRNAlength = 0.0
	ref = mRNA.chrom
	start =  mRNA.start
	stop = mRNA.stop
	
	try:
		genereads = bamfile.count(reference = ref, start = start,end = stop)
	except ValueError:
		continue
	
	Genecount.append(genereads)
	
	for CDS in db.children(mRNA, featuretype = "CDS"):
		CDSlength = CDS.stop - CDS.start +1
		mRNAlength += CDSlength
	Genelen.append(mRNAlength)
FPKM = []
for i in range(len(Genelen)):
	FPKM.append((Genecount[i]/(Genelen[i]*mappedreads))*10**9)

for elem in FPKM:
	print mRNA["Name"], elem


