import gffutils
# to open this database in any scripts use the code "db = gffutils.FeatureDB("yeast.db")
db = gffutils.FeatureDB("yeast.db")

mRNAs = []
allmRNA = 0.0
for mRNA in db.features_of_type("mRNA"):
	allmRNA += 1 
	for intron in db.children(mRNA, featuretype = "intron"):
		mRNAs.append( mRNA["Name"])

intronmRNA =  len(mRNAs)
		
print intronmRNA/allmRNA


for mRNA in db.features_of_type("mRNA"):
	mRNAlength = 0.0
	for CDS in db.children(mRNA, featuretype = "CDS"):
		CDSlength = CDS.stop - CDS.start +1
		mRNAlength += CDSlength
	#print mRNA["Name"], mRNAlength
	#print "%s\t%d"%(mRNA["Name"][0], mRNAlength)
	print "Gene %s has a length of %d" % (mRNA["Name"][0], mRNAlength)
