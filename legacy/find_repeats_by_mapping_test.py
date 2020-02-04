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
SMALT_DIR=""

##########################
# Error message function #
##########################

def DoError(errorstring):
	print "\nError:", errorstring
	print "\nFor help use -h or --help\n"
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


##############################
# Get command line arguments #
##############################

def get_user_options():
	usage = "usage: %prog [options]"
	version="%prog 1.0. Written by Simon Harris, Wellcome Trust Sanger Institute, 2012"
	parser = OptionParser(usage=usage, version=version)
	
	#do not allow arguments to be interspersed. e.g. -a -b arg1 agr2. MUST be -a arg1 -b arg2.
	parser.disable_interspersed_args()
	
	parser.add_option("-q", "--query", action="store", dest="query", help="Query sequence to search for in your reference", default="", metavar="FILE")
	parser.add_option("-r", "--reference", action="store", dest="ref", help="Reference sequence in which to find positions of query", default="", metavar="FILE")
	parser.add_option("-f", "--forwardfastq", action="store", dest="ffastq", help="Forward fastq file (can be gzipped, but must end in .gz)", default="")
	parser.add_option("-R", "--reversefastq", action="store", dest="rfastq", help="Reverse fastq file (can be gzipped, but must end in .gz)", default="")
	parser.add_option("-o", "--output", action="store", dest="output", help="Output file prefix", default="")
	parser.add_option("-i", "--qid", action="store", dest="queryid", help="Identity required for mapping to the query (between 0 and 1) [default= %default]", default=0.99, type='float')
	parser.add_option("-I", "--rid", action="store", dest="refid", help="Identity required for mapping reads back to the reference (between 0 and 1) [default= %default]", default=0.90, type='float')
	parser.add_option("-d", "--distance", action="store", dest="distance", help="Remap reads to reference if they are within this distance (bp) from the end of the query sequence. 0=include all", default=250, type='float')	
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):

	if options.output=='':
		options.output=options.ref.split("/")[-1].split(".")[0]+"_"+options.query.split("/")[-1].split(".")[0]
			

	if options.ref=='':
		DoError('No reference dna file (-r) selected!')
	if not os.path.isfile(options.ref):
		DoError('Cannot find file '+options.ref)
	if options.query=='':
		DoError('No query dna file (-q) selected!')
	if not os.path.isfile(options.query):
		DoError('Cannot find file '+options.query)
	if options.ffastq=='':
		DoError('No forward fastq file (-f) selected!')
	if not os.path.isfile(options.ffastq):
		DoError('Cannot find file '+options.ffastq)
	if options.rfastq=='':
		DoError('No reverse fastq file (-R) selected!')
	if not os.path.isfile(options.rfastq):
		DoError('Cannot find file '+options.rfastq)
	if options.distance<0:
		DoError('Distance options (-d) must be >=0')
	if options.queryid<0 or options.queryid>1:
		DoError('Query id must be between 0 and 1')
	if options.refid<0 or options.refid>1:
		DoError('Reference id must be between 0 and 1')
		
	
	return


def map_reads(freads="", rreads="", ref="", outputname="", maprepeats=False, percentid=0.9):
	
	if freads=="":
		print "No reads given"
		sys.exit()
	if ref=="":
		print "No reference given"
		sys.exit()
	if outputname=="":
		print "No output name given"
	
	#index the reference
	if not os.path.isfile(ref+".index.smi"):
		os.system(SMALT_DIR+"smalt index -k 13 -s 1 "+ref+".index "+ref)
	if not os.path.isfile(ref+".fai"):
		os.system(SAMTOOLS_DIR+"samtools faidx "+ref)
	
	#map the reads to the reference
	if maprepeats:
		os.system(SMALT_DIR+"smalt map -y "+str(percentid)+" -r 0 -f samsoft -o "+outputname+".sam "+ref+".index "+freads+" "+rreads)
	else:
		os.system(SMALT_DIR+"smalt map -y "+str(percentid)+" -f samsoft -o "+outputname+".sam "+ref+".index "+freads+" "+rreads)
	os.system(SAMTOOLS_DIR+"samtools view -F 4 -b -S -T "+ref+" -o "+outputname+".1.bam "+outputname+".sam")
	os.system(SAMTOOLS_DIR+"samtools sort "+outputname+".1.bam "+outputname)
	os.system(SAMTOOLS_DIR+"samtools index "+outputname+".bam")
	os.system("rm -f "+outputname+".1.bam "+outputname+".sam")


def rename_reads(inbam, outbam, contigs=[]):

	try: insamfile = pysam.Samfile( inbam, "rb" )
	except StandardError:
		print bamfile+" not a bam file"
		sys.exit()
	
	refs=insamfile.references
	lengths=insamfile.lengths
	
	myheader={}
	myheader["SQ"]=[]
	for x in xrange(len(refs)):
		myheader["SQ"].append({"SN": refs[x], "LN": lengths[x]})
	
	myheader["RG"]=[]
	for contig in contigs:
		myheader["RG"].append({"ID": contig, "SM": contig})
	
	try: outsamfile = pysam.Samfile( outbam, "wb", header=myheader )
	except StandardError:
		print "Could not create", outbam
		sys.exit()
		
	for read in insamfile:
		read.tags = read.tags + [("RG",read.qname.split("_")[-2])]
		if read.qname.split("_")[-1]=="F":
			read.is_read1=True
			read.is_read2=False
		else:
			read.is_read1=False
			read.is_read2=True
		read.qname="_".join(read.qname.split("_")[:-1])
		outsamfile.write(read)
	
	insamfile.close()
	outsamfile.close()
	os.system(SAMTOOLS_DIR+"samtools index "+outbam)

#############################
# Print read to output file #
#############################

def print_read_to_file(out, samread, refname, clippedseq="", clippedqual="", end="none"):
	
	
	if samread.is_reverse:
		if clippedseq!="":
			samreadseq=revcomp(clippedseq)
			samreadqual=clippedqual[::-1]
		else:
			samreadseq=revcomp(samread.seq)
			samreadqual=samread.qual[::-1]
	else:
		if clippedseq!="":
			samreadseq=clippedseq
			samreadqual=clippedqual
		else:
			samreadseq=samread.seq
			samreadqual=samread.qual
	
#	if samread.is_reverse:
#		samname=samname+"_R"
#	else:
#		samname=samname+"_F"
	samname="@"+samread.qname+"_"+refname	

	if end=="right":
		samname=samname+"_R"
		if not samread.is_reverse:
			samreadseq=revcomp(samreadseq)
			samreadqual=samreadqual[::-1]
	elif end=="left":
		samname=samname+"_F"
		if samread.is_reverse:
			samreadseq=revcomp(samreadseq)
			samreadqual=samreadqual[::-1]
	
	if samread.is_read1:
		samname=samname+"_1"
	elif samread.is_read2:
		samname=samname+"_2"
		
	print end, samname
	print >>out, samname
	
	print >>out, samreadseq
	
	print >> out, "+"
	print >>out, samreadqual
			
	
	return



def create_fastq_from_bam(bamfile, fastqfile):
	
	try: samfile = pysam.Samfile( bamfile, "rb" )
	except StandardError:
		print bamfile+" not a bam file"
		sys.exit()
	count=0
	fastqout=open(fastqfile+".fastq", "w")
	
	rlens=samfile.lengths
	refs=samfile.references
	reftolen={}
	for x in xrange(0,len(refs)):
		reftolen[refs[x]]=rlens[x]
	
	for read in samfile:
			
		#print read.qname, read.is_read1, read.is_read2
#		if read.is_unmapped and not read.mate_is_unmapped:
#			if options.distance==0 or (read.mate_is_reverse and read.pnext<options.distance) or (not read.mate_is_reverse and (read.pnext+read.rlen)>(reftolen[samfile.getrname(read.rnext)]-options.distance)):
#				count+=1
#				print_read_to_file(fastqout, read, samfile.getrname(read.rnext))
		
		if not read.is_unmapped:
			if read.pos==0 and read.cigar[0][0]==4 and read.cigar[0][1]>=10:
				count+=1
				if read.is_reverse:
					print_read_to_file(fastqout, read, samfile.getrname(read.tid), clippedseq=read.seq[:read.cigar[0][1]], clippedqual=read.qual[:read.cigar[0][1]], end="left")
				else:
					print_read_to_file(fastqout, read, samfile.getrname(read.tid), clippedseq=read.seq[:read.cigar[0][1]], clippedqual=read.qual[:read.cigar[0][1]], end="left")
			elif read.pos>(reftolen[samfile.getrname(read.tid)]-read.rlen) and read.cigar[-1][0]==4 and read.cigar[-1][1]>=10:
				count+=1 
				#need to check that read maps right up to end of element
				if read.is_reverse:
					print_read_to_file(fastqout, read, samfile.getrname(read.tid), clippedseq=read.seq[-1*read.cigar[-1][1]:], clippedqual=read.qual[-1*read.cigar[-1][1]:], end="right")
				else:
					print_read_to_file(fastqout, read, samfile.getrname(read.tid), clippedseq=read.seq[-1*read.cigar[-1][1]:], clippedqual=read.qual[-1*read.cigar[-1][1]:], end="right")
#			if options.distance==0 or (read.mate_is_reverse and read.pnext<options.distance) or (not read.mate_is_reverse and (read.pnext+read.rlen)>(reftolen[samfile.getrname(read.rnext)]-options.distance)):
#				count+=1
#				print_read_to_file(fastqout, read, samfile.getrname(read.rnext))
			
	
	fastqout.close()
	samfile.close()
	return count


def get_references(bamfile):
	try: samfile = pysam.Samfile( bamfile, "rb" )
	except StandardError:
		print bamfile+" not a bam file"
		sys.exit()
	refs=samfile.references
	samfile.close()
	return refs
	

########
# Main #
########


if __name__ == "__main__":



	(options, args)=get_user_options()
	check_input_validity(options, args)
	
	#make random name for files
	chars = string.ascii_letters + string.digits
	tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))
#	tmpname="tmpOPLFanZY"
	if options.ffastq.split(".")[-1]=="gz":
		print "Unzipping forward fastq file"
		os.system("zcat "+options.ffastq+" > "+tmpname+"_1.fastq")
		options.ffastq=tmpname+"_1.fastq"
	if options.rfastq.split(".")[-1]=="gz":
		print "Unzipping reverse fastq file"
		os.system("zcat "+options.rfastq+" > "+tmpname+"_2.fastq")
		options.rfastq=tmpname+"_2.fastq"
	
	
	map_reads(freads=options.ffastq, rreads=options.rfastq, ref=options.query, maprepeats=True, outputname=tmpname, percentid=options.queryid)
	
	query_contigs=get_references(tmpname+".bam")
	
	readcount=create_fastq_from_bam(tmpname+".bam", tmpname)
	
	if readcount>0:
		print readcount, "reads mapped to", options.query
		map_reads(freads=tmpname+".fastq", rreads="", ref=options.ref, outputname=tmpname+"_1", maprepeats=True, percentid=options.refid)
		rename_reads(tmpname+"_1.bam", options.output+".bam", contigs=query_contigs) 
	else:
		print "No reads mapped to", options.query
	
	#os.system("rm -f "+tmpname+"*")
	