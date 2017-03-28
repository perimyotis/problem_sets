
echo "# a_thing" >> README.md
git init
git add README.md
git commit -m "first commit"
git remote add origin https://github.com/perimyotis/a_thing.git
git push -u origin master

class FastaRecord(object):
	def _init_(self,title,sequence):
		self.title = title
		self.sequence = sequence

def read_fasta_records(infile):
	fastadict = {}
	for line in infile:

		if line[0] == ">":
			name = line.strip()[1:]
			fastadict[name] = []
			#I'm at a new sequence
			#so do stuff for that case
		
		else:
			fastadict[name].append(line.strip())
			#I'm reading lines from a sequence
			#append lines and stuff like that
	for sequence in fastafict:
		fastadict[sequenceq] = "".join(fastadict[sequence])
	return fastadict
	
		"title = line[1:].rstrip() #useful

		sequence_lines = []
		while 1:
			line =infile.readline().rstrip()
			if line[0] ==">" or line=="":
				break
			sequence_lines.append(line) #useful
		sequence = "".join(sequence_lines) #useful
		
		return FastaRecord(title,sequence)"
#in quotes just to save previous attempts

def comp(nuc):
	if nuc == "A":
		return "T"
	elif nuc == "T":
		return "A"
	elif nuc == "C":
		return "G"
	elif nuc == "G":
		return "C"
def revcomp(seq):
	revSeq = seq[::-1]
	rev = []
	for pos in revSeq:
		curComplement = complement(pos)
		revCom.append(curComplement)
	rev_com_seq = "".join(rev)
	return rev_com_seq

def stringsearch (substring, string):
	for i in range(len(string)-len(substring)):
		for x in range (len(substring)):
			if substring[x] != string [i+x]:
				break
			if x == len(substring)-1:
				return i
	return None

opened_file = open("testfasta.fa")

fasta_dictionary = read_fasta_records(infile = opened_file)

compdict = {}
for chrom in compdict:
	compdict[chrom] = revcomp(compdict[chrom])

fastq = open("my_reads.fastq")
Output = open("3.2output.txt","w")

read = 0
while 1:
	if read % 500 == 0:
		print read
	read +=1
	name2 = fastq.readline().strip()
	if name2 =="":
		break
	seq2 = fastq.readline().strip()
	name3 = fastq.readline().strip()
	qual = fastq.readline().strip()

	for chrom in compdict:
		map = stringsearch(seq2,compdict[chrom])
		if map is not None:
			Output.write("%s\t%d\t%s\n"%(chrom,map,name2))
			break
print "Done."
