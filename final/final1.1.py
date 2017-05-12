def read_fasta(fastaFile):
	fastaDict = {}
	chromNum = 0

	for line in fastaFile:
		if line[0] ==">":
			chromNum +=1
			chromName = chromNum
			fastaDict[chromName] = []
		else:
			fastaDict [chromName].append(line.strip())
	
	for chrom in fastaDict:
		fastaDict[chrom] = ''.join(fastaDict[chrom])
	return fastaDict

i = open("S_cerevisiae.fa")
ScerGenomeDict = read_fasta(i)

out = open("ScerTranscriptome.fa","w")

bedScer = open ("S_cerevisiae_genes.bed")
for line in bedScer:
	line = line.strip().split()
	if line [0] == 'chrI':
		line.insert(0,1)
	if line [0] == 'chrII':
		line.insert(0,2)
	if line [0] == 'chrIII':
		line.insert(0,3)
	if line [0] == 'chrIV':
		line.insert(0,4)
	if line[0] == 'chrV':
                line.insert(0,5)
	if line[0] == 'chrVI':
                line.insert(0,6)
	if line[0] == 'chrVII':
                line.insert(0,7)
	if line[0] == 'chrVIII':
                line.insert(0,8)
	if line[0] == 'chrIX':
                line.insert(0,9)
	if line[0] == 'chrX':
                line.insert(0,10)
	if line[0] == 'chrXI':
                line.insert(0,11)
	if line[0] == 'chrXII':
                line.insert(0,12)
	if line[0] == 'chrXIII':
                line.insert(0,13)
	if line[0] == 'chrXIV':
                line.insert(0,14)
	if line[0] == 'chrXV':
                line.insert(0,15)
	if line[0] == 'chrXVI':
                line.insert(0,16)
	geneName = line[4]
	geneSeq = ScerGenomeDict[line[0]][int(line[2]):(int(line[3])+1)]
	out.write(">%s\n%s\n"%(geneName,geneSeq))

def comp(nuc):
	if nuc == "A":
		return "T"
	elif nuc == "T":
		return" A"
	elif nuc == "C":
		return "G"
	elif nuc == "G":
		return "C"

def revcomp(seq):
	revSeq = seq[::-1]
	revComp = []
	for pos in revSeq:
		curComplement = complement(pos)
		revComp.append(curComplement)
	revComSeq = ''.join(revComp)
	return revComSeq

def fasta_read(x):
	names_linenum = []
	names = []
	genome = []
	chromNum2 = 0
	
	for i, line in enumerate(x):
		line = line.strip()
		genome.append(line)
		if line.staartswith(">"):
			chromNum2 +=1
			names_linenum.append(i)
			chromosome = line[1:]
			chromosome = chromosome.replace("_"," ")
			names.append(chromosome)
	names_linenum.append(len(genome))
	sequences = []
	for i in range (len(names_linenum)-1):
		seq = [''.join(genome[names_linenum[i]+1:names_linenum[i+1]])]
		sequences.exdtend(seq)

	dictionary = dict(zip(names,sequences))
	return dictionary 
