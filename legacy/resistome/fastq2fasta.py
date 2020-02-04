#!/usr/bin/env python


import string
import os, sys
from optparse import OptionParser
from Bio import SeqIO

##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options] <list of mfa files>"
	parser = OptionParser(usage=usage)

	parser.add_option("-f", "--fastq", action="store", dest="fastq", help="Input fastq file", default="")
	
	
	return parser.parse_args()



def DoError(errorstring):
	print "Error:", errorstring
	sys.exit()

################
# Main program #
################		

if __name__ == "__main__":

	#Get command line arguments

	(options, args) = main()
	
	if options.fastq=="":
		DoError("Fastq file (-f) required")
	
	if not os.path.isfile(options.fastq):
		DoError("Cannot find file "+options.fastq)
	
	

	
	
	count=SeqIO.convert(options.fastq, "fastq", "test.fasta", "fasta")

	print "converted", count, "reads"
