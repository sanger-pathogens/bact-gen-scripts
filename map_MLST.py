#!/usr/bin/env python

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
	parser.add_option("-g", "--genes", action="store", dest="genes", help="multifasta containing MLST genes to search for", default="", metavar="FILE")
	parser.add_option("-o", "--output", action="store", dest="output", help="output file prefix", default="", metavar="FILE")
	parser.add_option("-s", "--sts", action="store", dest="sts", help="st file", default="", metavar="FILE")


	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):

	if options.output=='':
		DoError("No output file selected")

	elif options.contigs=='':
		DoError("No contigs file selected")
	elif options.genes=='':
		DoError("No genes file selected")
	elif options.sts=='':
		DoError("No sts file selected")

	
	return


########
# Main #
########


if __name__ == "__main__":



	(options, args)=get_user_options()
	#print options, args
	check_input_validity(options, args)

	
	os.system("samtools faidx "+options.contigs)
	os.system("smalt index -k 13 -s 1 "+options.contigs+".index "+options.contigs)
	os.system("smalt map -y 0.9 -r "+str(randrange(1,99999))+" -f samsoft -o "+options.output+".sam "+options.contigs+".index "+options.genes)
	os.system("samtools view -b -S "+options.output+".sam -t "+options.contigs+".fai > "+options.output+".1.bam")
	os.system("samtools sort "+options.output+".1.bam "+options.output)
	os.system("samtools index "+options.output+".bam")
	os.system("rm -f "+options.output+".1.bam "+options.output+".sam "+options.contigs+".fai "+options.contigs+".index.*")
	os.system("~sh16/scripts/filter_MLST.py -c "+options.contigs+" -b "+options.output+".bam -g "+options.genes+" -o "+options.output+" -s "+options.sts)
