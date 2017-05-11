import os

samples = ["WT_2","WT_3","WT_4","WT_5","SNF2_1","SNF2_2","SNF2_3","SNF2_4","SNF2_5"]

for sample in samples:
	command = "kallisto quant -i Scere.idx -o ~/Desktop/sequencing_class/pset8/%s --single -l 180 -s 20 %s.fastq.gz"%(sample, sample)
	print "Currently running: %s"%command
	os.system(command)

