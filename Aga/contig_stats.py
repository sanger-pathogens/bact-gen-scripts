#!/usr/bin/env python
#/usr/bin/env python
import string
import os, sys
import pysam
from numpy import std, mean, median, max, min
from optparse import OptionParser


#################################
# Simple Error Printing Funtion #
#################################

def DoError(ErrorString):
	print "!!!Error:", ErrorString,"!!!"
	sys.exit()


#############################################
# Function to reverse complement a sequence #
#############################################

def revcomp(sequence):
	rev=sequence[::-1]
	revcomp=''
	d={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'g':'c', 'c':'g', "n":"n", "N":"N"}
	for i in rev:
		if d.has_key(i):
			revcomp=revcomp+d[i]
		else:
			revcomp=revcomp+i
	
	return revcomp


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-b", "--bam", action="store", dest="bam", help="input file in sam or bam format (must have suffix .sam or .bam)", default="", metavar="FILE")
	parser.add_option("-H", "--header", action="store_false", dest="header", help="Do not show header line", default=True, metavar="FILE")
	
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):

	if options.bam=='':
		DoError('No bam/sam input file selected')
	elif options.bam.split('.')[-1] not in ['sam', 'bam']:
		DoError('Input file '+options.bam+' must have the suffix .sam or .bam')
	elif not os.path.isfile(options.bam):
		DoError('Cannot find file '+options.bam)

	
	return








########
# Main #
########

if __name__ == "__main__":
	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)
	
	#open the bam/sam file
	
	if options.bam.split(".")[-1]=="bam":
		try:
			samfile = pysam.Samfile( options.bam, "rb" )
		except StandardError:
			DoError('Failed to open '+options.bam+'. Is it in bam format?')
	elif options.bam.split(".")[-1]=="sam":
		try:
			samfile = pysam.Samfile( options.bam, "r" )
		except StandardError:
			DoError('Failed to open '+options.bam+'. Is it in sam format?')
	
	#get reference names and lengths from the sam file header
	refs=samfile.references
	lengths=samfile.lengths
	
	if options.header:
		print '\t'.join(["contig", "length", "mean", "median", "std", "min", "max", "mapped%"])

	#get pileup depths for each reference
	refpileups={}
	for x, ref in enumerate(refs):
		depths=[]
		for pileupcolumn in samfile.pileup(ref):
			depths.append(pileupcolumn.n)
		refpileups[ref]={}
		refpileups[ref]["mean"]=mean(depths)
		refpileups[ref]["median"]=median(depths)
		refpileups[ref]["std"]=std(depths)
		refpileups[ref]["min"]=min(depths)
		refpileups[ref]["max"]=max(depths)
		refpileups[ref]["mapped%"]=(float(len(depths))/lengths[x])*100
		if refpileups[ref]["mapped%"]!=100:
			refpileups[ref]["min"]=0
		else:
			refpileups[ref]["min"]=min(depths)
		refpileups[ref]["length"]=lengths[x]
		print '\t'.join(map(str,[ref, refpileups[ref]["length"], refpileups[ref]["mean"], refpileups[ref]["median"], refpileups[ref]["std"], refpileups[ref]["min"], refpileups[ref]["max"], refpileups[ref]["mapped%"]]))
	
	sys.exit()
			
