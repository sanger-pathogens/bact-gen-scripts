#!/usr/bin/env python
import string, re
import os, sys
from optparse import OptionParser



####################
# Set some globals #
####################

SAMTOOLS_DIR=""
MY_SCRIPTS_DIR=""


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

	parser.add_option("-p", "--pileup", action="store", dest="pileup", help="pileup file", default="")
	parser.add_option("-b", "--bam", action="store", dest="bam", help="bam file", default="")
	parser.add_option("-s", "--sam", action="store", dest="sam", help="sam file", default="")
	parser.add_option("-q", "--quality", action="store", dest="quality", help="Mapping Quality", default=30, type="int")
	parser.add_option("-r", "--ratio", action="store", dest="ratio", help="SNP/Mapping quality ratio cutoff", default=0.75, type="float")
	
	parser.add_option("-o", "--output", action="store", dest="output", help="output file prefix", default="")
	
	return parser.parse_args()


################
# Main program #
################		

if __name__ == "__main__":


	#Get command line arguments

	(options, args) = main()


	#Do some checking of the input files
	
	if options.pileup=="":
		DoError("No pileup file specified")
	
	if options.output=="":
		DoError("No output prefix specified")
		
	if options.sam=="" and options.bam=="":
		DoError("sam or bam file from which pileup was made must be specified")
	elif options.sam!="":
		header=os.popen(SAMTOOLS_DIR+"samtools view -S -H "+options.sam).readlines()
	else:
		header=os.popen(SAMTOOLS_DIR+"samtools view -H "+options.bam).readlines()
	
	contigsizes={}
	
	for line in header:
		words=line.split()
		name=""
		length=0
		for word in words:
			if word.split(":")[0]=="SN":
				name=word.split(":")[1]
			elif word.split(":")[0]=="LN":
				length=int(word.split(":")[1])
		if name!="" and length!=0:
			contigsizes[name]=length
				
		
	if len(contigsizes)==0:
		DoError("No contigs found. Perhaps your sam/bam has no header?")
		
	
	contigs={}
	
	for contig in contigsizes.keys():
		contigs[contig]=["-"]*contigsizes[contig]
	

	pileupfile=open(options.pileup, "rU")
	
	for line in pileupfile:
		words=line.split()
		if words[3].upper()in ["A", "C", "G", "T"] and int(words[4])>options.quality and (words[2].upper()==words[3].upper() or (float(words[5])/float(words[4]))>=options.ratio) and words[2]!="*":
			contigs[words[0]][int(words[1])-1]=words[3]
	
			
	oneseq=""
	if len(contigs)>1:
		out=open(options.output+".mfa","w")
		for contig in contigs.keys():
			print >> out, ">"+''.join(contig)
			print >> out, ''.join(contigs[contig])
			oneseq=oneseq+''.join(contigs[contig])
		out.close()
	else:
		oneseq=oneseq+''.join(contigs[contigs.keys()[0]])
	
	
	
	out=open(options.output+".dna","w")	
	print >> out, ">"+options.output.split("/")[-1].split(".")[0]
	print >> out, oneseq



