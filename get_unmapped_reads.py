#!/usr/bin/env python
#/usr/bin/env python
import string
import os, sys
import pysam
from optparse import OptionParser
from Bio import SeqIO


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

	parser.add_option("-b", "--bam", action="store", dest="bam", help="mapping input file in sam or bam format (must have suffix .sam or .bam)", default="", metavar="FILE")
	parser.add_option("-f", "--forward", action="store", dest="forward", help="input forward (or single end) fastq file", default="", metavar="FILE")
	parser.add_option("-r", "--reverse", action="store", dest="reverse", help="input reverse fastq file (optional)", default="", metavar="FILE")
	parser.add_option("-o", "--output", action="store", dest="output", help="output prefix for file name(s)", default="", metavar="FILE")
	parser.add_option("-t", "--type", action="store", type="choice", dest="outputtype", choices=["unmapped","unpaired"], help="reads to output (unmapped or unpaired). Unmapped = both reads unmapped. Unpaired = either read is unmapped.", default="unmapped", metavar="FILE")
	
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):

	if options.bam=='':
		DoError('No bam/sam input file selected')
	elif options.bam.split('.')[-1] not in ['sam', 'bam']:
		DoError('Input file '+options.bam+' must have the suffix .sam or .bam')
	elif not os.path.isfile(options.bam):
		DoError('Cannot find file '+options.bam)
	if options.forward=='':
		DoError('No forward fastq input file selected')
	elif not os.path.isfile(options.forward):
		DoError('Cannot find file '+options.forward)
	if options.reverse!="" and not os.path.isfile(options.reverse):
		DoError('Cannot find file '+options.reverse)
	if options.output=="":
		DoError('No output file name (-o) selected')
	
	return




########
# Main #
########

if __name__ == "__main__":
	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)
	
	#open the bam/sam file
	
	if options.bam.split(".")[-1]=="bam":
		print "Reading bam file"
		try:
			samfile = pysam.Samfile( options.bam, "rb" )
		except StandardError:
			DoError('Failed to open '+options.bam+'. Is it in bam format?')
	elif options.bam.split(".")[-1]=="sam":
		print "Reading sam file"
		try:
			samfile = pysam.Samfile( options.bam, "r" )
		except StandardError:
			DoError('Failed to open '+options.bam+'. Is it in sam format?')


				
	print "Identifying reads that match criteria"
	toprint=set()
	add=toprint.add
	#Iterate through input file and print appropriate lines
	for read in samfile:
		
		if options.outputtype=="unmapped" and (not read.is_unmapped and not read.mate_is_unmapped):
			if not "/".join(read.qname.split("/")[:-1]) in toprint:
				add("/".join(read.qname.split("/")[:-1]))
		elif options.outputtype=="unpaired"  and not read.is_unmapped:
			if not "/".join(read.qname.split("/")[:-1]) in toprint:
				add("/".join(read.qname.split("/")[:-1]))
		
		
		
			   
	
	#clean up
	samfile.close
	
	print len(toprint), "reads found"
	
	print "Printing output files"
	
	fastqs = SeqIO.parse(open(options.forward, "rU"), "fastq")
	rfastqs = SeqIO.parse(open(options.reverse, "rU"), "fastq")
	outfastq=[]
	outrfastq=[]
	append=outfastq.append
	appendb=outrfastq.append
	for fastq in fastqs:
		rfastq=rfastqs.next()
		if not "/".join(fastq.id.split("/")[:-1]) in toprint:
			append(fastq)
			appendb(rfastq)
		else:
			toprint.remove("/".join(fastq.id.split("/")[:-1]))
	
	output_handle = open(options.output+"_1.fastq", "w")
	SeqIO.write(outfastq, output_handle, "fastq")
	output_handle.close()
	
	output_handle = open(options.output+"_2.fastq", "w")
	SeqIO.write(outrfastq, output_handle, "fastq")
	output_handle.close()
	
	

	    	
			