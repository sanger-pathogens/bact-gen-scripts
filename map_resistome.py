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
	parser.add_option("-f", "--forward", action="store", dest="forward", help="forward fasta file", default="", metavar="FILE")
	parser.add_option("-r", "--reverse", action="store", dest="reverse", help="reverse fasta file", default="", metavar="FILE")
	parser.add_option("-g", "--genes", action="store", dest="genes", help="multifasta containing genes to search for", default="", metavar="FILE")
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
	elif options.genes=='':
		DoError("No genes file selected")

	
	return


########
# Main #
########


if __name__ == "__main__":



	(options, args)=get_user_options()
	#print options, args
	check_input_validity(options, args)


#	os.system("samtools faidx "+options.genes)
#	os.system("smalt index -k 13 -s 1 "+options.genes+".index "+options.genes)
#	os.system("smalt map -y 0.8 -r "+str(randrange(1,99999))+" -f samsoft -o "+options.output+".sam "+options.genes+".index "+options.forward+" "+options.reverse)
#	
#	os.system("samtools view -S -h -F 4 "+options.output+".sam -t "+options.genes+".fai > "+options.output+".1.sam")
#	os.system("samtools view -S -f 4 -F 8 "+options.output+".sam -t "+options.genes+".fai >> "+options.output+".1.sam")
#	os.system("samtools view -b -h -S "+options.output+".1.sam -t "+options.genes+".fai > "+options.output+".1.bam")
#	os.system("samtools sort -n "+options.output+".1.bam "+options.output)
#
#	try:
#		genesfile=open(options.genes, "rU").read()
#	except StandardError:
#		print "Could not open contigs file"
#		sys.exit()
#	
#	
#	geneorder=[]
#	for line in genesfile.split(">")[1:]:
#		bits=line.split("\n")
#		geneorder.append(bits[0].split()[0])
#	os.system("rm -f tmpall.mfa")
#	count=0
#	for gene in geneorder:
#		print gene
#		os.system("~sh16/scripts/bam_filter.py -b test.1.bam -t contigs -c "+gene+" -o tmp")
#		os.system("~sh16/scripts/shufflefastqSequences.pl tmp_1.fastq tmp_2.fastq tmp_shuffled.fastq")
#		os.system("velveth velvet 21 -shortPaired -fastq tmp_shuffled.fastq")
#		os.system("velvetg velvet -cov_cutoff auto -exp_cov auto -min_contig_lgth 200")
#		os.system("cat velvet/contigs.fa >> tmpall.mfa")
#		os.system("rm -rf velvet")
#		count+=1
#		if gene=="ermC":
#			sys.exit()
#	
#	options.contigs="tmpall.mfa"

	
	os.system("samtools faidx "+options.contigs)
	os.system("smalt index -k 13 -s 1 "+options.contigs+".index "+options.contigs)
	os.system("smalt map -y 0.9 -r "+str(randrange(1,99999))+" -f samsoft -o "+options.output+".sam "+options.contigs+".index "+options.genes)
	os.system("samtools view -b -S "+options.output+".sam -t "+options.contigs+".fai > "+options.output+".1.bam")
	os.system("samtools sort "+options.output+".1.bam "+options.output)
	os.system("samtools index "+options.output+".bam")
	os.system("rm -f "+options.output+".1.bam "+options.output+".sam "+options.contigs+".fai "+options.contigs+".index.*")
	os.system("~sh16/scripts/filter_resistome.py -c "+options.contigs+" -b "+options.output+".bam -g "+options.genes+" -o "+options.output)
