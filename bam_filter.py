#!/usr/bin/env python
#/usr/bin/env python
import string
import os, sys
import pysam
from optparse import OptionParser


#################################
# Simple Error Printing Funtion #
#################################

def DoError(ErrorString):
	print "!!!Error:", ErrorString,"!!!"
	sys.exit()


#############################################
# Function to reverse complement a sequence #
#############################################

def revcomp(sequence):
	rev=sequence[::-1]
	revcomp=''
	d={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'g':'c', 'c':'g', "n":"n", "N":"N"}
	for i in rev:
		if d.has_key(i):
			revcomp=revcomp+d[i]
		else:
			revcomp=revcomp+i
	
	return revcomp


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-b", "--bam", action="store", dest="bam", help="input file in sam or bam format (must have suffix .sam or .bam)", default="", metavar="FILE")
	parser.add_option("-o", "--output", action="store", dest="output", help="output prefix for file name(s)", default="", metavar="FILE")
	parser.add_option("-f", "--format", type="choice", dest="fileformat", choices=["fasta","fastq","pairedfastq","sam","bam"], action="store", help="output file format (sam, bam, fastq, pairedfastq, or fasta) [default=%default]", default="pairedfastq", metavar="FORMAT")
	parser.add_option("-t", "--type", action="store", type="choice", dest="outputtype", choices=["unmapped","mapped", "paired", "all","contigs", "minusdup", "onemapped"], help="reads to output (all, mapped, paired, unmapped, minusdup, onemapped or contigs). Note 1: for pairedfastq output, mapped includes all reads where one of the pair is mapped. paired includes only those reads where the reads are in a proper pair. For other output formats the mapped option will only print the individual mapped reads and paired requires the reads to be mapped in proper pairs. Note 2: the contigs option requires the -c option to be set (see below). [default=%default]", default="unmapped", metavar="FILE")
	parser.add_option("-c", "--contigs", action="store", dest="contigs", help="comma separated list of contigs for which mapped reads should be output. Note there should be no whitespace within the list", default="")
	
	
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
	if options.output=="":
		DoError('No output file name (-o) selected')
	if options.outputtype=="contigs" and len(options.contigs)==0:
		DoError('Contigs options (-c) must be used when contigs is selected as read output type')
	if len(options.contigs)>0:
		options.contiglist=options.contigs.split(',')
	
	if options.outputtype=="all" and options.fileformat in ["sam", "bam"]:
		DoError('You do not need to use this program to save all reads in bam/sam format! Use samtools view instead')

	
	return



#############################
# Print read to output file #
#############################

def print_read_to_file(out, samread, format, sammate=False, mateout=False):
	
	
	if samread.is_reverse:
		samreadseq=revcomp(samread.seq)
		samreadqual=samread.qual[::-1]
	else:
		samreadseq=samread.seq
		samreadqual=samread.qual
	
	if sammate and sammate.is_reverse:
		sammateseq=revcomp(sammate.seq)
		sammatequal=sammate.qual[::-1]
	elif sammate!=False:
		sammateseq=sammate.seq
		sammatequal=sammate.qual
		
	
	if format in ["sam", "bam"]:
		out.write(read)
	elif format=="fasta":
		print >>out, ">"+samread.qname
		print >>out, samreadseq
	elif format=="fastq":
		print >>out, "@"+samread.qname
		print >>out, samreadseq
		print >> out, "+"
		print >> out, samreadqual
	elif format=="pairedfastq" and sammate and mateout:
		if samread.is_read1 and sammate.is_read2:
			print >>out, "@"+samread.qname
			print >>out, samreadseq
			print >> out, "+"
			print >> out, samreadqual
			print >>mateout, "@"+sammate.qname
			print >>mateout, sammateseq
			print >> mateout, "+"
			print >> mateout, sammatequal
			
		elif samread.is_read2 and sammate.is_read1:
			print >>out, "@"+sammate.qname
			print >>out, sammateseq
			print >> out, "+"
			print >> out, sammatequal
			print >>mateout, "@"+samread.qname
			print >>mateout, samreadseq
			print >> mateout, "+"
			print >> mateout, samreadqual
			
	
	
	
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
	
	#get reference names and lengths from the sam file header
	refs=samfile.references
	lengths=samfile.lengths
	
	
	#if the contigs option has been chosen, check that the contigs specified are in the sam header
	if options.outputtype=="contigs":
		for contig in options.contiglist:
			if not contig in refs:
				options.contiglist.remove(contig)
		if len(options.contiglist)==0:
			DoError('None of the contigs you specified are in the bam/sam header')

	
	#Create headers for sam/bam output
	if options.fileformat in ["sam", "bam"]:
		print "Adding sam headers"
		newrefs=[]
		newlengths=[]
		if options.outputtype=="contigs":
			
			for x, ref in enumerate(refs):
				if ref in options.contiglist:
					newrefs.append(ref)
					newlengths.append(lengths[x])
			
		else:
			newrefs=refs
			newlengths=lengths

			

	
	#open the output file
	if options.fileformat=="fasta":
		output=open(options.output+".fasta", "w")
	elif options.fileformat=="fastq":
		output=open(options.output+".fastq", "w")
	elif options.fileformat=="pairedfastq":
		output=open(options.output+"_1.fastq", "w")
		routput=open(options.output+"_2.fastq", "w")
	elif options.fileformat=="sam":
		output=pysam.Samfile(options.output+".sam", mode='wh', referencenames=newrefs, referencelengths=newlengths)
	elif options.fileformat=="bam":
		output=pysam.Samfile(options.output+".bam", mode='wb', referencenames=newrefs, referencelengths=newlengths)
	else:
		DoError('Somehow gained a file format')

		
	print "Converting file"
	if options.fileformat=="pairedfastq":
		firstofpair={}
	if options.outputtype=="minusdup":
		prevdup=False
		prevtoprint=False
		remdups=set()
	#Iterate through input file and print appropriate lines
	for read in samfile:
		
		if options.outputtype=="minusdup" and not prevdup:
			prevdup=read
		
		if options.fileformat=="pairedfastq" or options.outputtype=="minusdup":
			if len(read.qname.split("/"))>1:
				readname="/".join(read.qname.split("/")[:-1])
			else:
				readname=read.qname
			
			if options.outputtype=="minusdup":
				if len(prevdup.qname.split("/"))>1:
					prevname="/".join(prevdup.qname.split("/")[:-1])
				else:
					prevname=prevdup.qname
		
		printread=False
		printprev=False
		
		if options.outputtype=="all":
			if options.fileformat=="pairedfastq":
				if readname in firstofpair:
					printread=True
				else:
					firstofpair[readname]=read
			else:
				printread=True
			
		elif read.is_unmapped and options.outputtype=="unmapped":
			if options.fileformat=="pairedfastq":
				if readname in firstofpair:
					printread=True
				else:
					firstofpair[readname]=read
			else:
				printread=True
		
		elif (not read.is_unmapped or not read.mate_is_unmapped) and options.outputtype=="onemapped":
			if options.fileformat=="pairedfastq":
				if readname in firstofpair:
					printread=True
				else:
					firstofpair[readname]=read
			else:
				printread=True
			
		elif not read.is_unmapped and (options.outputtype=="mapped" or (options.outputtype=="contigs" and refs[read.rname] in options.contiglist)):
			if options.fileformat=="pairedfastq":
				if not read.mate_is_unmapped and read.mrnm==read.rname:
					if readname in firstofpair:
						printread=True
					else:
						firstofpair[readname]=read
			else:
				printread=True
		
		elif not read.is_proper_pair and options.outputtype=="paired":
			if options.fileformat=="pairedfastq":
				if not read.mate_is_unmapped and read.mrnm==read.rname:
					if readname in firstofpair:
						printread=True
					else:
						firstofpair[readname]=read
			else:
				printread=True
		
		elif not read.is_unmapped and options.outputtype=="minusdup":
			
			if  readname in remdups:
				remdups.remove(readname)
			
			elif not read.is_reverse and prevdup.pos==read.pos:
				if prevdup.mate_is_unmapped and not read.mate_is_unmapped:
					prevdup=read
					remdups.add(prevname)
				elif prevdup.mate_is_unmapped and not read.mate_is_unmapped:
					remdups.add(readname)
					continue
				elif read.is_proper_pair and not prevdup.is_proper_pair:
					if prevname in firstofpair:
						del firstofpair[prevname]
					remdups.add(prevname)
					prevdup=read
				elif prevdup.is_proper_pair and not read.is_proper_pair:
					remdups.add(readname)
					continue
				elif read.mapq>=prevdup.mapq:
					if prevname in firstofpair:
						del firstofpair[prevname]
					remdups.add(prevname)
					prevdup=read
				else:
					remdups.add(readname)
					continue
					
			elif prevdup.pos!=read.pos or read.is_reverse:
				if options.fileformat=="pairedfastq":
					if not prevdup.mate_is_unmapped and prevdup.mrnm==read.rname:
						if prevname in firstofpair:
							printprev=True
							prevtoprint=prevdup
						else:
							firstofpair[prevname]=prevdup
				else:
					printprev=True
				prevdup=read
				
		
		if printread:
			if options.fileformat=="pairedfastq":
				
				print_read_to_file(output, read, options.fileformat,sammate=firstofpair[readname],mateout=routput)
				del firstofpair[readname]
				
			else:
				print_read_to_file(output, read, options.fileformat)
				
		if printprev:
			if options.fileformat=="pairedfastq":
				print_read_to_file(output, prevtoprint, options.fileformat,sammate=firstofpair[prevname],mateout=routput)
				del firstofpair[prevname]
			else:
				print_read_to_file(output, prevdup, options.fileformat)
	   
	#clean up
	if options.fileformat=="pairedfastq":
		output.close()
		routput.close()
	elif options.fileformat in ["fasta", "fastq"]:
		output.close()
	else:
		output.close
	samfile.close
	    	
			