#!/usr/bin/env python

#Makes a list of CDS/intergenic regions from GFF of CDSs and a list of SNP locations 



#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------

import string, re
import os, sys, commands, getopt, gzip


def Usage():
	print 'getcdssnpsfromgff.py Usage:'
	print 'getcdssnpsfromgff.py -i [input alignment] -o [output file name] {-h}'
	print 'or'
	print 'getcdssnpsfromgff.py --in [input alignment] --out [output file name] --help'
	print 'Copyright Simon R Harris, Wellcome Trust Sanger Institute, UK. 2008'



#------------------------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------------------------

def getOptions(arg):


	try:
		opts, args = getopt.getopt(argv, "hg:s:o:", ["help", "gff=", "out=", "snp="])
	except getopt.GetoptError:
		Usage()
		sys.exit(2)

	gfffile=''
	outfile=''
	snpfile=''

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-g", "--gff"):
			gfffile=arg
		elif opt in ("-s", "--snp"):
			snpfile=arg
                elif opt in ("-o", "--out"):
			outfile=arg


	if gfffile=='':
		print 'no input gff file selected'
		Usage()
		sys.exit()
	if snpfile=='':
		print 'no input snp file selected'
		Usage()
		sys.exit()
	if outfile=='':
		outfile=gfffile.replace(".gff","")+"_cdssnps.out"
	
	return gfffile, snpfile, outfile


class progressBar:
	def __init__(self, minValue = 0, maxValue = 10, totalWidth=12):
		self.progBar = "[]"   # This holds the progress bar string
		self.min = minValue
		self.max = maxValue
		self.span = maxValue - minValue
		self.width = totalWidth
		self.amount = 0       # When amount == max, we are 100% done 
		self.updateAmount(0)  # Build progress bar string

	def updateAmount(self, newAmount = 0):
		if newAmount < self.min: newAmount = self.min
		if newAmount > self.max: newAmount = self.max
		self.amount = newAmount

		# Figure out the new percent done, round to an integer
		diffFromMin = float(self.amount - self.min)
		percentDone = (diffFromMin / float(self.span)) * 100.0
		percentDone = round(percentDone)
		percentDone = int(percentDone)

		# Figure out how many hash bars the percentage should be
		allFull = self.width - 2
		numHashes = (percentDone / 100.0) * allFull
		numHashes = int(round(numHashes))

		# build a progress bar with hashes and spaces
		self.progBar = "[" + '='*numHashes + '-'*(allFull-numHashes) + "]"

		# figure out where to put the percentage, roughly centered
		percentPlace = (len(self.progBar) / 2) - len(str(percentDone)) 
		percentString = str(percentDone) + "%"

		# slice the percentage into the bar
		self.progBar = self.progBar[0:percentPlace] + percentString + self.progBar[percentPlace+len(percentString):]

	def __str__(self):
		return str(self.progBar)




if __name__ == "__main__":
	argv=sys.argv[1:]
	gfffile, snpfile, outfile=getOptions(argv)

#Open the input file

lines=open(gfffile, "rU").readlines()

print "Reading gff file"

gffs=[]
for line in lines:
	for i in range(int(line.split()[3]), int(line.split()[4])+1):
		gffs.append(i)


lines=open(snpfile, "rU").readlines()

print "Reading SNP file"

output=open(outfile,"w")
count=0
total=0
prog = progressBar(0, len(lines), 77)
for line in lines:
	if int(line.strip()) in gffs:
		print >> output, 'g'
		gffs=gffs[gffs.index(int(line.strip())):]
	else:
		print >> output, 'i'
	output.flush()
	count=count+1
	if count==100:
		total=total+count
		prog.updateAmount(total)
		print prog, "\r",
		count=0


output.close()
print "\nDone."

