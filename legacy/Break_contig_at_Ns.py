#!/usr/bin/env python


SMALT_DIR=""
SAMTOOLS_DIR=""

##################
# Import modules #
##################

import os, sys
from optparse import OptionParser
from random import randrange

##########################
# Error message function #
##########################

def DoError(errorstring):
	print "\nError:", errorstring
	print "\nFor help use -h or --help\n"
	sys.exit()


##############################
# Get command line arguments #
##############################

def get_user_options():
	usage = "usage: %prog [options]"
	version="%prog 1.0. Written by Simon Harris, Wellcome Trust Sanger Institute, 2012"
	parser = OptionParser(usage=usage, version=version)
	
	#do not allow arguments to be interspersed. e.g. -a -b arg1 agr2. MUST be -a arg1 -b arg2.
	parser.disable_interspersed_args()
	

	parser.add_option("-c", "--contigs", action="store", dest="contigs", help="multifasta containing contigs to search in", default="", metavar="FILE")
	parser.add_option("-o", "--output", action="store", dest="output", help="output prefix", default="", metavar="FILE")


	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):

	if options.output=='':
		DoError("No output file selected")

	elif options.contigs=='':
		DoError("No contigs file selected")

	
	return


########
# Main #
########


if __name__ == "__main__":

	(options, args)=get_user_options()
	#print options, args
	check_input_validity(options, args)
	
	
	fastalines=open(options.contigs).read().split(">")
	fastaout=open(options.output, "w")
	for line in fastalines[1:]:
		name=line.split("\n")[0].split()[0]
		sequence="".join(line.split("\n")[1:])
		nons=sequence.replace("n"," ").replace("N"," ").split()
		if len(nons)>1:
			for x, seq in enumerate(nons):
				print >> fastaout, ">"+name+"_part"+str(x+1)
				print >> fastaout, nons[x]
				#print name+"_part"+str(x+1), len(nons[x])
		else:
			print >> fastaout, ">"+name
			print >> fastaout, nons[0]
			#print name, len(nons[0])
	fastaout.close()