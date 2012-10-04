#!/usr/bin/env python


##################
# Import modules #
##################

import string, re, os, sys

from optparse import OptionParser, OptionGroup
from random import *
from numpy import mean, max, min, median, std, sum
import pysam

###########
# Globals #
###########

SAMTOOLS_DIR=""
BCFTOOLS_DIR=""


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
	usage = "usage: %prog [options] <list of bam files>"
	version="%prog 1.0. Written by Simon Harris, Wellcome Trust Sanger Institute, 2011"
	parser = OptionParser(usage=usage, version=version)
	
	#do not allow arguments to be interspersed. e.g. -a -b arg1 agr2. MUST be -a arg1 -b arg2.
	parser.disable_interspersed_args()
#	
	parser.add_option("-r", "--reference", action="store", dest="ref", help="Reference DNA sequence (in fasta or multi-fasta format)", default="", metavar="FILE")
#	parser.add_option("-b", "--bam", action="store", dest="bam", help="Input bam file", default="")
#	parser.add_option("-B", "--bcf", action="store", dest="bcf", help="Input bcf/vcf file", default="")
#	parser.add_option("-v", "--vcf", action="store_true", dest="vcf", help="variant file is in vcf format", default=False)
#	parser.add_option("-o", "--output", action="store", dest="output", help="Output file prefix", default="")
#	parser.add_option("-t", "--tempname", action="store", dest="tmpname", help="Prefix for temporary files", default="")
#	parser.add_option("-p", "--proportion", action="store", dest="proportion", help="Proportion of reads required to support alternate for it to be accepted", default=0.75, type=float)
#	parser.add_option("-d", "--depth", action="store", dest="depth", help="Minimum number of reads for an indel to be accepted", default=4, type=int)
#	parser.add_option("-P", "--pindel", action="store_false", dest="runpindel", help="Do not run pindel", default=True)
#	
#	parser.add_option("-D", "--directory", action="store", dest="wd", help="Working directory for analysis", default="")
	
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):

#	if options.output=='':
#		options.output=options.ref.split("/")[-1].split(".")[0]+"_realigned"
			

	if options.ref=='':
		DoError('No reference dna file (-r) selected!')
	if not os.path.isfile(options.ref):
		DoError('Cannot find file '+options.ref)
#	if options.bam=='':
#		DoError('No bam file (-b) selected!')
#	if not os.path.isfile(options.bam):
#		DoError('Cannot find file '+options.bam)
#	if options.bcf!='' and not os.path.isfile(options.bcf):
#		DoError('Cannot find file '+options.bcf)
#	if options.proportion>1 or options.proportion<0:
#		DoError('Proportion parameter must be between 0 and 1')
#	if options.wd!="" and not os.path.isdir(options.wd):
#		DoError('Working directory '+options.wd+" does not exist")
#	
#	try: samfile = pysam.Samfile( options.bam, "rb" )
#	except StandardError:
#		print options.bam+" not a bam file"
#		sys.exit() 
#	
#	refs=samfile.references
#	lengths=samfile.lengths
#	samfile.close()
#	
#	if not os.path.isfile(options.ref):
#		DoError('Cannot find file '+options.ref)
		
	
#	return refs, lengths




########
# Main #
########


if __name__ == "__main__":



	(options, args)=get_user_options()
	#print options, args
	check_input_validity(options, args)
	
	#make random name for files
#	chars = string.ascii_letters + string.digits
#	if options.tmpname=="":
#		tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))
#	else:
#		tmpname=options.tmpname
#	
#	if options.wd!="":
#		tmpname=options.wd+"/"+tmpname
	
	print "Reading reference contigs"
	
	refdata={}
	contigorder=[]
	reffile=open(options.ref, "rU")
	lines=reffile.read().split('>')[1:]
	for line in lines:
		words=line.strip().split('\n')
		refdata[words[0].split()[0]]=map(list,[[]]*len(''.join(words[1:]).upper()))
		contigorder.append(words[0].split()[0])
	reffile.close()
	
	
	for filename in args:
		print "Reading "+filename
		try:
			bcffile=os.popen(BCFTOOLS_DIR+"bcftools view "+filename)
		except StandardError:
			DoError("Cannot open bcf file")
		lastbase=0
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
			
			BASEINFO={}
			
			for x, heading in enumerate(headings):
			
				if heading=="INFO":
					
					BASEINFO[heading]={}
					
					try: info=words[x].split(";")
					except StandardError:
						print "Cannot split info string", words[x]
						sys.exit()
					for i in info:
						
						infotype=i.split("=")[0]
						
						if len(i.split("="))<2:
							if infotype=="INDEL":
								BASEINFO[heading][infotype]=True
						else:
							infodata=i.split("=")[1]
							try: BASEINFO[heading][infotype]=float(infodata)
							except StandardError:
								try: BASEINFO[heading][infotype]=map(float,infodata.split(","))
								except StandardError:
									BASEINFO[heading][infotype]=infodata
					
					
						
				else:
					try: BASEINFO[heading]=float(words[x])
					except StandardError:
						BASEINFO[heading]=words[x]
			
			if "INDEL" in BASEINFO['INFO']:
				continue
			if "DP4" in BASEINFO['INFO']:
				if (BASEINFO['INFO']["DP4"][0]+BASEINFO['INFO']["DP4"][1])>(BASEINFO['INFO']["DP4"][2]+BASEINFO['INFO']["DP4"][3]):
					try:
						refdata[BASEINFO["CHROM"]][int(BASEINFO["POS"])].append((BASEINFO['INFO']["DP4"][2]+BASEINFO['INFO']["DP4"][3])/(BASEINFO['INFO']["DP4"][0]+BASEINFO['INFO']["DP4"][1]))
					except IndexError:
						print BASEINFO["POS"], len(refdata[BASEINFO["CHROM"]])
				else:
					try:
						refdata[BASEINFO["CHROM"]][int(BASEINFO["POS"])].append((BASEINFO['INFO']["DP4"][0]+BASEINFO['INFO']["DP4"][1])/(BASEINFO['INFO']["DP4"][2]+BASEINFO['INFO']["DP4"][3]))
					except IndexError:
						print BASEINFO["POS"], len(refdata[BASEINFO["CHROM"]])
			else:
				refdata[BASEINFO["CHROM"]][int(BASEINFO["POS"])].append(0)
				
#			print lastbase, int(BASEINFO["POS"])
			x=lastbase
			while x<int(BASEINFO["POS"])+1:
				refdata[BASEINFO["CHROM"]][x].append(0)
				x+=1
			
			lastbase=int(BASEINFO["POS"])
			if lastbase>20000:
				break
				
		bcffile.close()
	
	output=open("test.plot","w")
	print >> output, "#BASE MEAN MEDIAN MIN MAX"
	for contig in contigorder:
		for x, base in enumerate(refdata[contig]):
			if len(base)>0 and max(base)>0:
				print >> output, x+1, mean(base), median(base), min(base), max(base)
	output.close()
