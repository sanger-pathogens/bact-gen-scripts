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


##############################################
## Function to reverse complement a sequence #
##############################################

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
	
	parser.add_option("-q", "--query", action="store", dest="query", help="Query sequence to search for in your reference", default="", metavar="FILE")
	parser.add_option("-r", "--reference", action="store", dest="ref", help="Reference sequence in which to find positions of query", default="", metavar="FILE")
	parser.add_option("-f", "--forwardfastq", action="store", dest="ffastq", help="Forward fastq file (can be gzipped, but must end in .gz)", default="")
	parser.add_option("-R", "--reversefastq", action="store", dest="rfastq", help="Reverse fastq file (can be gzipped, but must end in .gz)", default="")
	parser.add_option("-o", "--output", action="store", dest="output", help="Output file prefix", default="")
	parser.add_option("-i", "--qid", action="store", dest="queryid", help="Identity required for mapping to the query (between 0 and 1) [default= %default]", default=0.99, type='float')
	parser.add_option("-I", "--rid", action="store", dest="refid", help="Identity required for mapping reads back to the reference (between 0 and 1) [default= %default]", default=0.90, type='float')
	parser.add_option("-d", "--distance", action="store", dest="distance", help="Remap reads to reference if they are within this distance (bp) from the end of the query sequence. 0=include all", default=250, type='float')	
	parser.add_option("-e", "--rep_query", action="store_false", dest="repquery", help="Do not map repeats randomly when mapping to query sequence", default=True)
	parser.add_option("-E", "--rep_ref", action="store_false", dest="repref", help="Do not map repeats randomly when mapping back to reference sequence", default=True)		
	
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


def map_reads(freads="", rreads="", ref="", outputname="", maprepeats=False, percentid=0.9, onlyunmappedpairs=True):
	
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
		os.system(SMALT_DIR+"smalt map -y "+str(percentid)+" -r 0 -f bam -o "+outputname+".bam "+ref+".index "+freads+" "+rreads)
	else:
		os.system(SMALT_DIR+"smalt map -y "+str(percentid)+" -r -1 -f bam -o "+outputname+".bam "+ref+".index "+freads+" "+rreads)
	if onlyunmappedpairs:
		os.system(SAMTOOLS_DIR+"samtools view -b -f 4 -F 8 -o "+outputname+".1.bam "+outputname+".bam")
	else:
		os.system(SAMTOOLS_DIR+"samtools view -b -o "+outputname+".1.bam "+outputname+".bam")
	os.system(SAMTOOLS_DIR+"samtools sort "+outputname+".1.bam "+outputname)
	os.system(SAMTOOLS_DIR+"samtools index "+outputname+".bam")
	os.system("rm -f "+outputname+".1.bam "+outputname+".sam")


def rename_reads(inbam, outbam, idtoref, contigs=[]):

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
		read.tags = read.tags + [("RG",idtoref[int(read.qname.split("_")[-2])])]
		if read.qname.split("_")[-1]=="F":
			read.is_read1=True
			read.is_read2=False
		else:
			read.is_read1=False
			read.is_read2=True
		read.qname="_".join(read.qname.split("_")[:-2]+[idtoref[int(read.qname.split("_")[-2])]])
		outsamfile.write(read)
	
	insamfile.close()
	outsamfile.close()
	os.system(SAMTOOLS_DIR+"samtools index "+outbam)

#############################
# Print read to output file #
#############################

def print_read_to_file(out, samread, refname, reftoid, startpos=0, endpos=-1):
	
	if endpos==-1:
		endpos=samread.rlen
	
	samreadseq=str(samread.seq)[startpos:endpos]
	samreadqual=str(samread.qual)[startpos:endpos]
	
	refconverter={}
	if samread.is_reverse:
		samreadseq=revcomp(samreadseq)
		samreadqual=samreadqual[::-1]
		
	if samread.is_read1:
		samname="@"+samread.qname+"_1"
	elif samread.is_read2:
		samname="@"+samread.qname+"_2"
		
	#This is to avoid problems with names containing _ or other odd characters
	samname=samname+"_"+str(reftoid[refname])
	
	if samread.mate_is_reverse:
		samname=samname+"_F"
	else:
		samname=samname+"_R"
	print >>out, samname
	print >>out, samreadseq
	print >> out, "+"
	print >> out, samreadqual
			
	
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
	idtoref={}
	reftoid={}
	for x in xrange(0,len(refs)):
		reftolen[refs[x]]=rlens[x]
		idtoref[x]=refs[x]
		reftoid[refs[x]]=x
	
	toprint={}
	
	for read in samfile:
			
		#print read.qname, read.is_read1, read.is_read2
		
		if not read.is_unmapped:
		
			if read.mate_is_unmapped and (options.distance==0 or (read.is_reverse and read.pos<options.distance) or (not read.is_reverse and (read.pos+read.rlen)>(reftolen[samfile.getrname(read.tid)]-options.distance))):
				if not read.qname in toprint:
					toprint[read.qname]=set([])
				if read.is_read1:
					toprint[read.qname].add(2)
				else:
					toprint[read.qname].add(1)
			elif read.is_proper_pair:
				print read.pnext, read.rlen, reftolen[samfile.getrname(read.tid)], read.mate_is_reverse
				if read.pnext==0 or (read.pnext+read.rlen)>=(reftolen[samfile.getrname(read.tid)]):
					if not read.qname in toprint:
						toprint[read.qname]=set([])
					if read.is_read1:
						toprint[read.qname].add(2)
					else:
						toprint[read.qname].add(1)
			
			if read.pos==0 or read.is_reverse and (read.pos+read.rlen)>=(reftolen[samfile.getrname(read.tid)]):
				if not read.qname in toprint:
					toprint[read.qname]=set([])
				if read.is_read1:
					toprint[read.qname].add(1)
				else:
					toprint[read.qname].add(2)
		
	
			
	
	print toprint		
	samfile.reset()
	
	for read in samfile:
		if read.is_read1:
			readnum=1
		else:
			readnum=2
		
		
		if read.is_unmapped and read.qname in toprint and readnum in toprint[read.qname]:
			count+=1
			print_read_to_file(fastqout, read, samfile.getrname(read.rnext), reftoid)
		
		elif read.qname in toprint and readnum in toprint[read.qname]:
			if options.distance==0 or read.pos<options.distance or (read.pos+read.rlen)>(reftolen[samfile.getrname(read.tid)]-options.distance):
				mapped=0.0
				matches=0.
				insertions=0.0
				deletions=0.0
				readpos=0
				total=read.rlen
				refpos=read.pos
				readseq=read.seq.upper()
				startpos=0
				endpos=read.rlen
				for cig in read.cigar:
					if cig[0]==0:
						for x in range(0,cig[1]):
							if readseq[readpos]==refseqs[samfile.getrname(read.tid)][refpos] and refseqs[samfile.getrname(read.tid)][refpos] in ["A", "C", "G", "T"]:
								matches+=1
							readpos+=1
							refpos+=1
					elif cig[0]==1:
						readpos+=cig[1]
					elif cig[0]==2:
						refpos+=cig[1]
					elif cig[0]==4:
						readpos+=cig[1]
					else:
						print cig
				
				
				if read.pos==0 and read.cigar[0][0]==4:
					total-=read.cigar[0][1]
					startpos=read.cigar[0][1]
				if read.cigar[-1][0]==4 and refpos==reftolen[samfile.getrname(read.tid)]:
					total-=read.cigar[-1][1]
					endpos=read.rlen-read.cigar[-1][1]
				percentid=matches/total
				
				#print percentid, read.rlen, matches, total, read.cigar, options.queryid, startpos, endpos, read.mate_is_unmapped, read.pos, read.cigar[0][0], refpos, read.cigar[-1][0], reftolen[samfile.getrname(read.tid)]
				
				if percentid>=options.queryid:
					#count+=1
					#print_read_to_file(fastqout, read, samfile.getrname(read.rnext), reftoid)
					if (read.pos==0 and read.cigar[0][0]==4) or (read.cigar[-1][0]==4 and refpos==reftolen[samfile.getrname(read.tid)]):
						count+=1
						print_read_to_file(fastqout, read, samfile.getrname(read.tid), reftoid, startpos=startpos, endpos=endpos)
						print percentid, read.rlen, matches, total, read.cigar, options.queryid, startpos, endpos, read.mate_is_unmapped, read.pos, read.cigar[0][0], refpos, read.cigar[-1][0], reftolen[samfile.getrname(read.tid)]
				
					
					
	
	fastqout.close()
	samfile.close()
	return count, idtoref


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
	
	if options.ffastq.split(".")[-1]=="gz":
		print "Unzipping forward fastq file"
		os.system("zcat "+options.ffastq+" > "+tmpname+"_1.fastq")
		options.ffastq=tmpname+"_1.fastq"
	if options.rfastq.split(".")[-1]=="gz":
		print "Unzipping reverse fastq file"
		os.system("zcat "+options.rfastq+" > "+tmpname+"_2.fastq")
		options.rfastq=tmpname+"_2.fastq"
	
	#tmpname="tmpOPLFanZY"
	map_reads(freads=options.ffastq, rreads=options.rfastq, ref=options.query, maprepeats=options.repquery, outputname=tmpname, percentid=0, onlyunmappedpairs=False)
	
	query_contigs=get_references(tmpname+".bam")
	reflines=open(options.query).read().split(">")[1:]
	
	refseqs={}
	for ref in reflines:
		seqlines=ref.split("\n")
		refname=seqlines[0].split()[0]
		refseqs[refname]=''.join(seqlines[1:]).upper()
		
	readcount,idtoref=create_fastq_from_bam(tmpname+".bam", tmpname)
	
	if readcount>0:
		print readcount, "reads mapped to", options.query
		map_reads(freads=tmpname+".fastq", rreads="", ref=options.ref, outputname=tmpname, maprepeats=options.repref, percentid=options.refid, onlyunmappedpairs=False)
		rename_reads(tmpname+".bam", options.output+".bam", idtoref, contigs=query_contigs) 
	else:
		print "No reads mapped to", options.query
	
	os.system("rm -f "+tmpname+"*")
	