import os

ScermRNA = ["Scer_RNA_seq_1.fastq.gz","Scer_RNA_seq_2.fastq.gz"]

for sample in ScermRNA:
	command = "kallisto quant -i ScerTranscriptome.idx -o ScermRNAOut%s --single -l 180 -s 20 %s" %(ScermRNA.index(sample), sample)
	print "Currently running: %s"% command
	os.system(command)

ScerRibo = ["Scer_ribo_seq_1.fastq.gz","Scer_ribo_seq_2.fastq.gz"]

for sample in ScerRibo:
	command = "kallisto quant -i ScerTranscriptome.idx -o ScerRiboOut%s --single -l 180 -s 20 %s" %(ScerRibo.index(sample), sample)
	print "Currently running: %s"% command
	os.system(command)

SparmRNA = ["Spar_RNA_seq_1.fastq.gz", "Spar_RNA_seq_2.fastq.gz"]

for sample in SparmRNA:
	command = "kallisto quant -i SparTranscriptome.idx -o SparmRNAOut%s --single -l 180 -s 20 %s" %(SparmRNA.index(sample),sample)
	print"Currently running %s"% command
	os.system(command)

SparRibo = ["Spar_ribo_seq_1.fastq.gz","Spar_ribo_seq_2.fastq.gz"]

for sample in SparRibo:
	command = "kallisto quant -i SparTranscriptome.idx -o SparRiboOut%s --single -l 180 -s 20 %s" %(SparRibo.index(sample),sample)
	print"Currently running %s"% command
	os.system(command)
