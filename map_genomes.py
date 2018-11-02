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
	

	parser.add_option("-r", "--ref", action="store", dest="ref", help="mutifasta containing reference contigs", default="", metavar="FILE")
	parser.add_option("-q", "--query", action="store", dest="query", help="multifasta containing query contigs", default="", metavar="FILE")
	parser.add_option("-o", "--output", action="store", dest="output", help="output file prefix", default="", metavar="FILE")


	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):

	if options.output=='':
		DoError("No output file selected")

	elif options.ref=='':
		DoError("No reference file selected")
	elif options.query=='':
		DoError("No query contigs file selected")

	
	return


########
# Main #
########


if __name__ == "__main__":



	(options, args)=get_user_options()
	#print options, args
	check_input_validity(options, args)

	
	os.system("samtools faidx "+options.ref)
	os.system("smalt index -k 13 -s 1 "+options.ref+".index "+options.ref)
	os.system("/nfs/pathogen/sh16_scripts/fasta2fastq_shredder.py "+options.query+" "+options.output+" 1000 200 l 1000")
	os.system("smalt map -i 5000 -j 1000 -y 0.7 -r "+str(randrange(1,99999))+" -f samsoft -o "+options.output+".sam "+options.ref+".index "+options.output+"_1.fastq "+options.output+"_2.fastq ")
	os.system("samtools view -b -S "+options.output+".sam -t "+options.ref+".fai > "+options.output+".1.bam")
	os.system("samtools sort "+options.output+".1.bam "+options.output)
	os.system("samtools index "+options.output+".bam")
	os.system("rm -f "+options.output+".1.bam "+options.output+".sam "+options.ref+".fai "+options.ref+".index.*")
	#os.system("/nfs/pathogen/sh16_scripts/filter_MLST.py -c "+options.ref+" -b "+options.output+".bam -g "+options.genes+" -o "+options.output)
