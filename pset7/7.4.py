import gffutils as gu
import pysam
import scipy.stats as st
import numpy as np 

db = gu.FeatureDB("yeast.db")

bamfile_BY = pysam.AlignmentFile("BYReads.sorted.bam", "rb")
bamfile_RM = pysam.AlignmentFile("RMReads.sorted.bam","rb")

out = open("7.4.txt","w")

GenecountBY = []
GenecountRM = []
Genelen = []
GeneName = []
mappedreadsBY = 0.0
mappedreadsRM = 0.0
for mRNA in db.features_of_type("mRNA"):
	GeneName.append(mRNA["Name"][0])

	mRNAlength = 0.0
	ref = mRNA.chrom
	start = mRNA.start
	stop = mRNA.stop

	try:
		genereadsBY = bamfile_BY.count(reference = ref, start = start, end = stop)
		genereadsRM = bamfile_RM.count(reference = ref, start = start, end = stop)
	except ValueError:
		continue
	GenecountBY.append(genereadsBY)
	GenecountRM.append(genereadsRM)

	mappedreadsBY += genereadsBY
	mappedreadsRM += genereadsRM

	for CDS in db.children(mRNA, featuretype = "CDS"):
		CDSlength = CDS.stop - CDS.start +1
		mRNAlength += CDSlength
	Genelen.append(mRNAlength)

out.write("Name\tLength\tBY_counts\tRM_counts\n")
for i in range(len(Genelen)):
	k11 = GenecountBY[i]
	N1 = mappedreadsBY

	k21 = GenecountRM[i]
	N2 = mappedreadsRM

	lambda1 = ((k11)/N1+(k21)/N2)/2

	lambda11 = (k11)/N1
	lambda21 = (k21)/N2

	prob_data_null = st.poisson.pmf(k11,N1*lambda1)*st.poisson.pmf(k21,N2*lambda1)
	prob_data_alt = st.poisson.pmf(k11,N1*lambda11)*st.poisson.pmf(k21,N2*lambda21)

	LRT = 2*(np.log(prob_data_alt) - np.log(prob_data_null))

	p_value = st.chi2.sf(LRT,df=1)


	out.write("%s\t%d\t%d\t%d\n"%(GeneName[i],Genelen[i],GenecountBY[i],GenecountRM[i]))

print "Done"
