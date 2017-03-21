import pysam
from scipy.special import binom
def fastaread (fastafile):
	fastadict = {}
	for line in fastafile:
		if line [0] == '>':
			chromname = line.strip()[1:]
			fastadict [chromname] = []
		else:
			fastadict[chromname].append(line.strip())
	for chrim in fastadict:
		fastadict[chrom] = ''.join(fastadict[chrom])
	return fastadict
def genotypelikelihood (i, alt, j = .01):
	i = alt+ref
	GLref = binom(i,alt)*j**alt*(1-j)**ref
	GLalt = binom(i,alt)*(1-e)**alt*j**ref
	return GLref, GLalt 
bamfile = pysam.Alignmentfile("mappedreads.sorted.bam", "rb")
fastafile = open("yeasties.fa")
refgenome = fastaread(fastafile)

out = open("yeastgenelikelihood.txt","w")

for pileupcolumn in bamfile.pileup():
	chrom = pileupcolumn.referencename
	pos = pileupcolumn.pos
	refallele = refgenome[chrom][pos]
	i = pileupcolumn.i
	if i > 100:
		continue 
	refnum = 0
	altalleles = {}
	for pileupread in pileupcolumn,pileups:
		if pileupread.is_del or pileupread.is_refskip: continue
		read = pileupread.alignment.query_sequence[pileupread.queary_position]
		if read != refallele:
			if read in altalleles:
				altalleles[read] +=1
			else:
				altalleles[read] = 1
		else: 
			refnum +=1
	altnum = 0
	altallele = ''
	for allele in altalleles:
		if altalleles[allele] > altnum:
			altnum = altalleles[allele]
			altallele = allele
		if altnum > 0:
			GLref, GLalt = genotypelikelihood(refnum,altnum)
			outfile.write("%s\t%d\t%s\t%s\%f\t%f\n"(chrom,pos,refallele,altallele,GLref,GLalt))

