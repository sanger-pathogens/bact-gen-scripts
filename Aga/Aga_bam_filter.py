#!/usr/bin/env python
#/usr/bin/env python
import string
import os, sys
import pysam
from numpy import std, mean, median, max, min
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
	parser.add_option("-t", "--type", action="store", type="choice", dest="outputtype", choices=["aga", "bothunmapped", "oneunmapped", "atleastoneunmapped", "mapped", "paired", "all","contigs", "minusdup", "notpaired", "unmappedmateofmapped"], help="reads to output (aga, all, mapped, paired, notpaired, oneunmapped, atleastoneunmapped, bothunmapped, minusdup or contigs). Note 1: for pairedfastq output, mapped includes all reads where one of the pair is mapped. paired includes only those reads where the reads are in a proper pair. For other output formats the mapped option will only print the individual mapped reads and paired requires the reads to be mapped in proper pairs. Note 2: the contigs option requires the -c option to be set (see below). [default=%default]", default="oneunmapped", metavar="FILE")
	parser.add_option("-c", "--contigs", action="store", dest="contigs", help="comma separated list of contigs for which mapped reads should be output. Note there should be no whitespace within the list", default="")
	parser.add_option("-T", "--tab", action="store", dest="tab", help="tab file containing regions to treat as unmapped when using the aga option", default="")
	
	
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
	if options.outputtype=="unmappedmateofmapped" and options.fileformat=="pairedfastq":
		options.fileformat="fastq"
	
	if options.tab!="" and not os.path.isfile(options.tab):
		DoError('Cannot find file '+options.tab)
		

	
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
	elif sammate:
		sammateseq=sammate.seq
		sammatequal=sammate.qual
		
	
	if format in ["sam", "bam"]:
		out.write(read)
	elif format=="fasta":
		print >>out, ">"+samread.qname
		print >>out, samreadseq
	elif format=="fastq":
		if not samread.is_paired:
			samname="@"+samread.qname
		if samread.is_read1:
			samname="@"+samread.qname+"/1"
		elif samread.is_read2:
			samname="@"+samread.qname+"/2"
		print >>out, samname
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
	
	
	
	regions=[]
	if options.tab!="":
		try:
			tablines=open(options.tab,"rU").readlines()
		except StandardError:
			print "Failed to read tab file"
			sys.exit()
		ftloc=-1
		note=""
		for line in tablines:
#			         111111111122222222223333333333
#			123456789012345678901234567890123456789
#			FT   misc_feature    2825393..2826315
			line=line.strip()
			if len(line)>22 and line[:2]=="FT" and len(line[2:21].strip())>0:
				try: ftloc=map(int,(line[21:].replace("complement(","").replace(")","").split("..")))
				except StandardError:
					print "Failed to read tab file"
					print line[21:].split("..")
					sys.exit()
				note=""
			elif len(line)>22 and line[:2]=="FT" and len(line[2:21].strip())==0 and ftloc!=-1:
				if line[21:].split("=")[0]=="/note":
					note='='.join(line[22:].split("=")[1:])
					if len(ftloc)==1:
						regions.append([ftloc[0], ftloc[0], note])
					else:
						regions.append([ftloc[0], ftloc[1], note])
						
					ftloc=-1
	
	
	
	
	
	
	
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
	
	contig_start_pos=[]
	curlen=0
	for x, ref in enumerate(refs):
		contig_start_pos.append(curlen)
		curlen+=lengths[x]
	
	
	if options.outputtype=="contigs":
		use_contigs=True
	else:
		use_contigs=False
	
	#if the contigs option has been chosen, check that the contigs specified are in the sam header
	if use_contigs:
		for contig in options.contiglist:
			if not contig in refs:
				options.contiglist.remove(contig)
		if len(options.contiglist)==0:
			DoError('None of the contigs you specified are in the bam/sam header')
	
	
	
	#Go through sam file and identify which contig each tab region is in
	
	contig_regions={}

	contigstart=0
	added=0
	for x, contig in enumerate(refs):
		contigend=contigstart+lengths[x]
		contig_regions[contig]=[]
		for region in regions:
			if region[0]>=contigstart and region[0]<contigend and region[1]>=contigstart and region[1]<contigend:
				contig_regions[contig].append(region)
				added+=1
		
		if added!=len(regions):
			DoError("Regions span contig breaks")
			
		
		contigstart=contigend
		
	reads_to_remove=set([])
	for contig in contig_regions:
		for region in xrange(len(contig_regions[contig])):
			for read in samfile.fetch(contig, contig_regions[contig][region][0], contig_regions[contig][region][1]):
				reads_to_remove.add(read.qname)
	
	print len(reads_to_remove)
	
	
	#Create headers for sam/bam output
	if options.fileformat in ["sam", "bam"]:
		print "Adding sam headers"
		newrefs=[]
		newlengths=[]
		if use_contigs:
			
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
	count=0
	
	sicount=0
	allcount=0
	samfile.reset()
	#Iterate through input file and print appropriate lines
	for read in samfile:
		allcount+=1
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
				
		
		elif options.outputtype=="unmappedmateofmapped" and read.is_unmapped and not read.mate_is_unmapped:
			printread=True
			
		elif (read.is_unmapped and read.mate_is_unmapped) and (options.outputtype=="bothunmapped" or options.outputtype=="aga" or options.outputtype=="atleastoneunmapped"):
			if options.fileformat=="pairedfastq":
				if readname in firstofpair:
					printread=True
				else:
					firstofpair[readname]=read
			else:
				printread=True
		
		elif (read.is_unmapped or read.mate_is_unmapped) and (options.outputtype=="oneunmapped" or options.outputtype=="aga" or options.outputtype=="atleastoneunmapped"):
			if options.fileformat=="pairedfastq":
				if readname in firstofpair:
					printread=True
				else:
					firstofpair[readname]=read
			else:
				printread=True
		
		elif options.outputtype=="aga" and read.qname in reads_to_remove:
			if options.fileformat=="pairedfastq":
				if readname in firstofpair:
					printread=True
				else:
					firstofpair[readname]=read
			else:
				printread=True
			sicount+=1
		
		elif not read.is_unmapped and (options.outputtype=="mapped" or (use_contigs and refs[read.rname] in options.contiglist) ):
			if options.fileformat=="pairedfastq":
				if not read.mate_is_unmapped and read.mrnm==read.rname:
					if readname in firstofpair:
						printread=True
					else:
						firstofpair[readname]=read
			else:
				printread=True

		elif not read.is_proper_pair and options.outputtype=="notpaired":
			if options.fileformat=="pairedfastq":
				if readname in firstofpair:
					printread=True
				else:
					firstofpair[readname]=read
			else:
				printread=True
				
		elif read.is_proper_pair and (options.outputtype=="paired" or (use_contigs and refs[read.rname] in options.contiglist)):
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
				count+=2
				del firstofpair[readname]
				
			else:
				print_read_to_file(output, read, options.fileformat)
				count+=1
				
		if printprev:
			if options.fileformat=="pairedfastq":
				print_read_to_file(output, prevtoprint, options.fileformat,sammate=firstofpair[prevname],mateout=routput)
				count+=2
				del firstofpair[prevname]
			else:
				print_read_to_file(output, prevdup, options.fileformat)
				count+=1
	   
	#clean up
	if options.fileformat=="pairedfastq":
		output.close()
		routput.close()
	elif options.fileformat in ["fasta", "fastq"]:
		output.close()
	else:
		output.close()
	samfile.close()
	print sicount, count, allcount
			
