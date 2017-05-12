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

i = open ("S_paradoxus.fa")
SparGenomeDict = read_fasta(i)

out = open("SparTranscriptome.fa","w")

bedSpar = open ("S_paradoxus_genes.bed")
for line in bedSpar:
	line = line.strip().split()
	if line[0] == 'Spar_1':
		line.insert(0,1)
	if line[0] == 'Spar_2':
                line.insert(0,2)
	if line[0] == 'Spar_3':
                line.insert(0,3)
	if line[0] == 'Spar_4':
                line.insert(0,4)
	if line[0] == 'Spar_5':
                line.insert(0,5)
	if line[0] == 'Spar_6':
                line.insert(0,6)
	if line[0] == 'Spar_7':
                line.insert(0,7)
	if line[0] == 'Spar_8':
                line.insert(0,8)
	if line[0] == 'Spar_9':
                line.insert(0,9)
	if line[0] == 'Spar_10':
                line.insert(0,10)
	if line[0] == 'Spar_11':
                line.insert(0,11)
	if line[0] == 'Spar_12':
                line.insert(0,12)
	if line[0] == 'Spar_13':
                line.insert(0,13)
	if line[0] == 'Spar_14':
                line.insert(0,14)
	if line[0] == 'Spar_15':
                line.insert(0,15)
	if line[0] == 'Spar_16':
                line.insert(0,16)
	stop = line[4].find('.')
	if stop == -1:
		geneName = line[4]
	else:
		geneName = line[4][:stop]
	geneSeq = SparGenomeDict[line[0]][int(line[2]):(int(line[3])+1)]
	out.write(">%s\n%s\n"%(geneName,geneSeq))

i.close()
out.close()	
