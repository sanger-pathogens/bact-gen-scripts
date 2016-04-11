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
	parser.add_option("-f", "--forward", action="store", dest="forward", help="forward fastq file (may be zipped, but must end .gz)", default="", metavar="FILE")
	parser.add_option("-r", "--reverse", action="store", dest="reverse", help="reverse fastq file (may be zipped, but must end .gz)", default="", metavar="FILE")
	parser.add_option("-g", "--genes", action="store", dest="genes", help="multifasta containing genes to search for", default="", metavar="FILE")
	parser.add_option("-o", "--output", action="store", dest="output", help="output prefix", default="", metavar="FILE")
	parser.add_option("-i", "--id", action="store", dest="id", help="minimum id to report match (excluding clipping due to contig breaks). Must be between 0 and 1. [default = %default]", default=0.9, type="float", metavar="float")


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
	if options.id<0 or options.id>1:
		print "minimum id (-i) must be between 0 and 1"
		sys.exit()

	
	return


########
# Main #
########


if __name__ == "__main__":



	(options, args)=get_user_options()
	#print options, args
	check_input_validity(options, args)
	
	try:
		genesfile=open(options.genes, "rU").read()
	except StandardError:
		print "Could not open genes file"
		sys.exit()
		
	genes=set([])
	for line in genesfile.split(">")[1:]:
		bits=line.split("\n")
		if bits[0].split()[0] in genes:
			print "Error: your genes file contains more than one gene with the same name:", bits[0].split()[0]
			sys.exit()
		genes.add(bits[0].split()[0])
	
	
	fastalines=open(options.contigs).read().split(">")
	fastaout=open(options.output+"_no_Ns.fasta", "w")
	for line in fastalines[1:]:
		name=line.split("\n")[0].split()[0]
		sequence="".join(line.split("\n")[1:])
		nons=sequence.replace("n"," ").replace("N"," ").split()
		if len(nons)>1:
			for x, seq in enumerate(nons):
				if len(nons[x])>13:#remove contigs smaller than has length - will dies in smalt otherwise
					print >> fastaout, ">"+name+"_part"+str(x+1)
					print >> fastaout, nons[x]
					#print name+"_part"+str(x+1), len(nons[x])
		else:
			print >> fastaout, ">"+name
			print >> fastaout, nons[0]
			#print name, len(nons[0])
	fastaout.close()

	os.system(SAMTOOLS_DIR+"samtools faidx "+options.output+"_no_Ns.fasta")
	os.system(SMALT_DIR+"smalt index -k 13 -s 1 "+options.output+"_no_Ns.fasta.index "+options.output+"_no_Ns.fasta")
	os.system(SMALT_DIR+"smalt map -d -1 -f samsoft -o "+options.output+".sam "+options.output+"_no_Ns.fasta.index "+options.genes)
	os.system(SAMTOOLS_DIR+"samtools view -b -S "+options.output+".sam -t "+options.output+"_no_Ns.fasta.fai > "+options.output+".1.bam")
	os.system(SAMTOOLS_DIR+"samtools sort "+options.output+".1.bam "+options.output)
	os.system(SAMTOOLS_DIR+"samtools index "+options.output+".bam")
	os.system("rm -f "+options.output+".1.bam "+options.output+".sam "+options.output+"_no_Ns.fasta.fai "+options.output+"_no_Ns.fasta.index.sma "+options.output+"_no_Ns.fasta.index.smi ")#+options.output+"_no_Ns.fasta.fai "+options.output+"_no_Ns.fasta.index.*")
	os.system("~sh16/scripts/filter_resistome.py -i "+str(options.id)+" -c "+options.output+"_no_Ns.fasta -b "+options.output+".bam -g "+options.genes+" -o "+options.output+" -f "+options.forward+" -r "+options.reverse)
