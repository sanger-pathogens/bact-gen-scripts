#!/usr/bin/env python

SMALT_DIR=""
SAMTOOLS_DIR=""
BAM_FILTER_LOCATION="~sh16/scripts/Aga/bam_filter.py"
HUMAN_GENOME_LOCATION="/lustre/scratch101/blastdb/Ensembl/Human/NCBI36/genome/softmasked_dusted.fa"

##################
# Import modules #
##################

import os, sys, string
from random import randint, choice
from optparse import OptionParser


##############################
# Get command line arguments #
##############################

def get_user_options():
	usage = "usage: %prog [options]"
	version="%prog 1.0. Written by Simon Harris, Wellcome Trust Sanger Institute, 2012"
	parser = OptionParser(usage=usage, version=version)
	
	#do not allow arguments to be interspersed. e.g. -a -b arg1 agr2. MUST be -a arg1 -b arg2.
	parser.disable_interspersed_args()
	
	parser.add_option("-o", "--output", action="store", dest="output", help="prefix for output files", default="", metavar="FILE")
	parser.add_option("-f", "--forward", action="store", dest="forward", help="forward fastq file (may be zipped, but must end .gz)", default="", metavar="FILE")
	parser.add_option("-r", "--reverse", action="store", dest="reverse", help="reverse fastq file (may be zipped, but must end .gz)", default="", metavar="FILE")
	parser.add_option("-i", "--id", action="store", dest="id", help="minimum id for mapping (excluding clipping due to contig breaks) [default = %default]", default=0.9, type="float", metavar="float")

	

	return parser.parse_args()

(options, args)=get_user_options()

if not os.path.isfile(options.forward):
	print "Could not find forward fastq file"
	sys.exit()	
if not os.path.isfile(options.reverse):
	print "Could not find reverse fastq file"
	sys.exit()
if options.output=="":
	options.output="Filtered"

if options.id<0 or options.id>1:
	print "percent id (-i) must be between 0 and 1"
	sys.exit()


chars = string.ascii_letters + string.digits
tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))

os.system("cp -s "+HUMAN_GENOME_LOCATION+ " "+tmpname+".ref.fasta")
os.system(SAMTOOLS_DIR+"samtools faidx "+tmpname+".ref.fasta")
os.system(SMALT_DIR+"smalt index -k 20 -s 13 "+tmpname+".ref.fasta.index "+tmpname+".ref.fasta")
os.system(SMALT_DIR+"smalt map -y "+str(options.id)+" -r 12345 -f samsoft -o "+tmpname+".sam "+tmpname+".ref.fasta.index "+options.forward+" "+options.reverse)
os.system(BAM_FILTER_LOCATION+" -t bothunmapped -f pairedfastq -o "+options.output+" -b "+tmpname+".sam")
os.system("rm -rf "+tmpname+"*")

