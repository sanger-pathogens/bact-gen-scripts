#!/usr/bin/env python

#takes an input fasta file (alignment?) and outputs a tab file showing conserved blocks (between Ns) of a size range between a max and min specified

import string
import os, sys
from optparse import OptionParser
from Bio import SeqIO

#################################
# Simple Error Printing Funtion #
#################################

def DoError(ErrorString):
	print "!!!Error:", ErrorString,"!!!"
	sys.exit()


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-f", "--fasta", action="store", dest="fasta", help="fasta input file (containing Ns)", default="", metavar="FILE")
	parser.add_option("-m", "--min", action="store", dest="minimum", help="minimum block size to output", type=int, default=0, metavar="INT")
	parser.add_option("-M", "--max", action="store", dest="maximum", help="maximum block size to output", type=int, default=0, metavar="INT")
	parser.add_option("-o", "--output", action="store", dest="output", help="output tab file name", default="", metavar="FILE")
	
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):

	if options.fasta=='':
		DoError('No fasta input file selected')
	elif not os.path.isfile(options.fasta):
		DoError('Cannot find file '+options.fasta)
	if options.output=="":
		DoError('No output file name (-o) selected')
	if options.minimum<0:
		DoError('Minimum block size must be greater than 0')
	if options.maximum<=options.minimum and options.maximum!=0:
		DoError('Maximum block size must be greater than minimum blocks size (or 0 for no maximum)')
	
	return
	
########
# Main #
########

if __name__ == "__main__":
	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)
	
	try:
		sequences = SeqIO.parse(open(options.fasta, "rU"), "fasta")
	except StandardError:
		DoError('Failed to open fasta file')
	
	
	
	
	output=open(options.output, "w")
	
	total=0
	for sequence in sequences:
		seq=sequence.seq
		runlen=0
		runstart=total
		for base in seq.upper():
			total+=1
			if base=="N":
				if runlen>=options.minimum and (runlen<=options.maximum or options.maximum<=options.minimum):
					print >> output, "FT   misc_feature    "+str(runstart)+".."+str(total-1)

				runlen=0
				runstart=total+1
				
				
			else:
				runlen+=1
			
	output.close()
			
			
	
