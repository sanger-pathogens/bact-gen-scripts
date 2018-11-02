#!/usr/bin/env python
import string, re
import os, sys
from optparse import OptionParser, OptionGroup
#probably should add this to make sequences in real fasta format?
#from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord



####################
# Set some globals #
####################

SAMTOOLS_DIR=""
BCFTOOLS_DIR=""
SSAHA_DIR="/nfs/users/nfs_s/sh16/ssaha2_v2.5.1_x86_64/"
BWA_DIR=""
MY_SCRIPTS_DIR="/nfs/pathogen/sh16_scripts/"


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
	
	group = OptionGroup(parser, "IO Options")
	group.add_option("-b", "--bcf", action="store", dest="bcf", help="bcf/vcf file", default="", metavar="file")
	group.add_option("-o", "--output", action="store", dest="output", help="output file name", default="", metavar="prefix")
	parser.add_option_group(group)
	
	
	parser.add_option_group(group)
	
	
	
	
	return parser.parse_args()


################
# Main program #
################		

if __name__ == "__main__":


	#Get command line arguments

	(options, args) = main()


	#Do some checking of the input files
	
	if options.bcf=="":
		DoError("No input file specified")
	elif not os.path.isfile(options.bcf):
		DoError("Cannot find input file")
	
	if options.output=="":
		DoError("No output prefix specified")
		
	
	try:
		bcffile=os.popen(BCFTOOLS_DIR+"bcftools view "+options.bcf)
	except StandardError:
		DoError("Cannot open bcf file")
	
			
	mapplot=open(options.output, "w")

	prevbase=0
	for line in bcffile:

		words=line.split()
			
		if words[0][0]=="#":
			if words[0][1]!="#":

				headings=words
				headings[0]=headings[0][1:]
			continue

		
		if len(words)!=len(headings):
			print "words not equal to headings"
			print headings
			print words
			sys.exit()
		
		coverage=0.0
		base=0
		indel=False
		
		for x, heading in enumerate(headings):
			if heading=="POS":
				base=int(words[x])
			elif heading=="INFO":
				
				try: info=words[x].split(";")
				except StandardError:
					print "Cannot split info string", words[x]
					sys.exit()
				for i in info:
					if i=="INDEL":
						indel=True
					infotype=i.split("=")[0]
					if infotype=="DP":
						coverage=int(i.split("=")[1])
		if indel:
			continue
		prevbase+=1
		while prevbase<base:
			print >> mapplot, 0
			prevbase+=1
		print >> mapplot, coverage
		
		
	mapplot.close()


