#!/usr/bin/env python

import os, sys
from optparse import OptionParser, OptionGroup
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_SeqIO import *

##########################
# Error message function #
##########################

def DoError(errorstring):
	print "\nError:", errorstring
	print "\nFor help use -h or --help\n"
	sys.exit()


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	
	group = OptionGroup(parser, "General options")
	group.add_option("-i", "--input", action="store", dest="input", help="Input fasta file", default="", metavar="FILE")
	group.add_option("-x", "--xmfa", action="store", dest="xmfa", help="Name of output xmfa file", default="", metavar="FILE")
	parser.add_option_group(group)


	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):
	
	options.formatdb=True
	
	if options.input=='':
		DoError('No input file selected')
	elif not os.path.isfile(options.input):
		DoError('Cannot find file '+options.input)
	if options.xmfa=='':
		DoError('No output file name selected')
		
	return

################
# Main program #
################		

if __name__ == "__main__":
	
	
	#Get command line arguments

	(options, args) = main()
	
	
	try:
		seqs=read_seq_file(options.input)
		if len(seqs)==0:
			DoError("Found empy input file")	
	except StandardError:
		DoError("Cannot read input file")
	
	output=open(options.xmfa, "w")
	for seq in seqs:
		if ">" in str(seq.seq):
			DoError("Found > in sequence")
		
		seq.name=seq.name+":0-"+str(len(seq.seq))+" +"
		print >> output, ">"+seq.name
		print >> output, seq.seq
	print >> output, "="
	output.close()
	