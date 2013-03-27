#!/usr/bin/env python


SMALT_DIR=""
SAMTOOLS_DIR=""

##################
# Import modules #
##################

import os, sys, string
from optparse import OptionParser
from random import *

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
	

	parser.add_option("-l", "--length", action="store", dest="length", help="length of minimum repeat region to search for (e.g. read length)", default=76, metavar="INT")
	parser.add_option("-r", "--reference", action="store", dest="ref", help="reference fasta file", default="", metavar="FILE")
	parser.add_option("-o", "--output", action="store", dest="output", help="output prefix", default="", metavar="FILE")
	parser.add_option("-c", "--circular", action="store_true", dest="circular", help="contigs are circular", default=False, metavar="FILE")


	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):

	if options.output=='':
		DoError("No output file selected")

	elif options.ref=='':
		DoError("No reference file selected")

	
	return


########
# Main #
########


if __name__ == "__main__":

	step=1

	(options, args)=get_user_options()
	#print options, args
	check_input_validity(options, args)
	
	chars = string.ascii_letters + string.digits
	tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))
	
	if options.ref.split('.')[-1]=='gz':
		refcontigs=gzip.open(options.ref, 'r').read().split(">")[1:]
	else:
		refcontigs=open(options.ref, 'rU').read().split(">")[1:]
	
	
	refs=[]
	for refcontig in refcontigs:
		refs.append(''.join(refcontig.split("\n")[1:]))
	
	outputf=open(tmpname+".fastq", "w")
	for refcount, ref in enumerate(refs):
		count=0
		for x in range(0-int(options.length),len(ref),step):
			if x+int(options.length)>len(ref):
				continue
			count=count+1
			
			if x<0:
				print >> outputf, "@IL0_0000:1:1:"+str(refcount+1)+":"+str(count)+":"+str(len(ref)+x+1)
			else:
				print >> outputf, "@IL0_0000:1:1:"+str(refcount+1)+":"+str(count)+":"+str(x+1)
				
			if not options.circular:
				print >> outputf, "N"*int(options.length)
				print >> outputf, "+"
				print >> outputf, "&"*int(options.length)
					
			else:
				if (x+int(options.length))<0:
					print >> outputf, ref[(x+len(ref)):(x+int(options.length))]
				elif x<0:
					print >> outputf, ref[(x+len(ref)):]+ref[0:x+int(options.length)]
				else:
					print >> outputf, ref[x:x+int(options.length)]
				print >> outputf, "+"
				print >> outputf, "B"*int(options.length)
		
			
	outputf.close()


	os.system(SAMTOOLS_DIR+"samtools faidx "+options.ref)
	os.system(SMALT_DIR+"smalt index -k 13 -s 1 "+options.ref+".index "+options.ref)
	os.system(SMALT_DIR+"smalt map -r -1 -f samsoft -o "+tmpname+".sam "+options.ref+".index "+tmpname+".fastq")
	os.system(SAMTOOLS_DIR+"samtools view -f 4 -b -S "+tmpname+".sam -t "+tmpname+".fai > "+tmpname+".1.bam")
	os.system(SAMTOOLS_DIR+"samtools sort "+tmpname+".1.bam "+options.output)
	os.system(SAMTOOLS_DIR+"samtools index "+options.output+".bam")
	#os.system("rm -f "+tmpname+"*")
	
