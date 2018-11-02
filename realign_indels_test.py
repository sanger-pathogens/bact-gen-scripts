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
MY_SCRIPTS_DIR="/nfs/pathogen/sh16_scripts/"
BREAKDANCER_DIR="/nfs/users/nfs_s/sh16/breakdancer-1.1_2011_02_21/"
#PINDEL_DIR="/nfs/users/nfs_s/sh16/pindel/trunk/"
PINDEL_DIR="/nfs/users/nfs_s/sh16/pindel024t/"
DINDEL_DIR="/nfs/users/nfs_s/sh16/dindel/binaries/"
DINDEL_SCRIPTS_DIR="/nfs/users/nfs_s/sh16/dindel/dindel-1.01-python/"


##########################
# Error message function #
##########################

def DoError(errorstring):
	print "\nError:", errorstring
	print "\nFor help use -h or --help\n"
	sys.exit()


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
	usage = "usage: %prog [options] <list of bam files>"
	version="%prog 1.0. Written by Simon Harris, Wellcome Trust Sanger Institute, 2011"
	parser = OptionParser(usage=usage, version=version)
	
	#do not allow arguments to be interspersed. e.g. -a -b arg1 agr2. MUST be -a arg1 -b arg2.
	parser.disable_interspersed_args()
	
	parser.add_option("-r", "--reference", action="store", dest="ref", help="Reference DNA sequence (in fasta or multi-fasta format)", default="", metavar="FILE")
	parser.add_option("-b", "--bam", action="store", dest="bam", help="Input bam file", default="")
	parser.add_option("-B", "--bcf", action="store", dest="bcf", help="Input bcf/vcf file", default="")
	parser.add_option("-i", "--isbam", action="store", dest="iselements", help="Input bam IS element file", default="")
	parser.add_option("-I", "--issequence", action="store", dest="ISfasta", help="Input fasta file of IS element file", default="")
	parser.add_option("-v", "--vcf", action="store_true", dest="vcf", help="variant file is in vcf format", default=False)
	parser.add_option("-o", "--output", action="store", dest="output", help="Output file prefix", default="")
	parser.add_option("-t", "--tempname", action="store", dest="tmpname", help="Prefix for temporary files", default="")
	parser.add_option("-p", "--proportion", action="store", dest="proportion", help="Proportion of reads required to support alternate for it to be accepted", default=0.75, type=float)
	parser.add_option("-d", "--depth", action="store", dest="depth", help="Minimum number of reads for an indel to be accepted", default=4, type=int)
	parser.add_option("-P", "--pindel", action="store_false", dest="runpindel", help="Do not run pindel", default=True)
	
	parser.add_option("-D", "--directory", action="store", dest="wd", help="Working directory for analysis", default="")
	
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):

	if options.output=='':
		options.output=options.ref.split("/")[-1].split(".")[0]+"_realigned"
			

	if options.ref=='':
		DoError('No reference dna file (-r) selected!')
	if options.bam=='':
		DoError('No bam file (-b) selected!')
	if not os.path.isfile(options.bam):
		DoError('Cannot find file '+options.bam)
	if options.bcf!='' and not os.path.isfile(options.bcf):
		DoError('Cannot find file '+options.bcf)
	if options.proportion>1 or options.proportion<0:
		DoError('Proportion parameter must be between 0 and 1')
	if options.wd!="" and not os.path.isdir(options.wd):
		DoError('Working directory '+options.wd+" does not exist")
	
	try: samfile = pysam.Samfile( options.bam, "rb" )
	except StandardError:
		print options.bam+" not a bam file"
		sys.exit() 
	
	refs=samfile.references
	lengths=samfile.lengths
	samfile.close()
	
	if not os.path.isfile(options.ref):
		DoError('Cannot find file '+options.ref)
		
	
	return refs, lengths





def samfile_data(region_start, region_end, contig_seqs={}, readlist=[]):
	
	
	try: samfile = pysam.Samfile( tmpname+".bam", "rb" )
	except StandardError:
		print tmpname+".bam not a bam file"
		sys.exit() 
	
	refs=samfile.references
	lengths=samfile.lengths
	
#	tinsertions=0
#	tinslength=0
#	tdeletions=0
#	tdellength=0
#	tSNPs=0
#	tclipped=0
#	tcliplength=0
#	tmapped=0
	readinfo={}
	overlapping_reads=set([])
	read_name_dict={}
	count=0
	count2=0
	for read in samfile:
		
		readname=read.qname
		if read.is_read1:
			readname=readname+"/1"
		elif read.is_read2:
			readname=readname+"/2"
		
		count+=1
		if read.is_unmapped:
			readinfo[readname]=1000
			read_name_dict[readname]=read.qname
			continue
		start=read.pos
		readpos=0
		refpos=start
		insertions=0
#		inslength=0
		deletions=0
#		dellength=0
		SNPs=0
		clipped=0
#		cliplength=0
#		mapped=0
		for cig in read.cigar:
			if cig[0]==0:
				for x in range(0,cig[1]):
					try:
						if read.seq[readpos].upper()!=contig_seqs[refs[read.tid]][refpos].upper():
							SNPs+=1
					except IndexError:
						print read.pos, readpos, refpos, read.seq, read.tid, refs[read.tid], len(contig_seqs[refs[read.tid]]), read.seq[readpos], contig_seqs[refs[read.tid]][refpos]
						sys.exit()
						#print read.seq[readpos].upper(), contigs[refs[read.tid]][refpos].upper(),
#					else:
#						mapped+=1
					readpos+=1
					refpos+=1
			elif cig[0]==1:
				insertions+=1
#				inslength+=cig[1]
				readpos+=cig[1]
			elif cig[0]==2:
				deletions+=1
#				dellength+=cig[1]
				refpos+=cig[1]
			elif cig[0]==4:
				clipped+=1
#				cliplength+=cig[1]
				readpos+=cig[1]
			else:
				print cig
		end=refpos
		count2+=1
#		print readname, start, end, region_start, region_end
		if (start<region_start and end>region_start) or (start<region_end and end>region_end):
			
			overlapping_reads.add(readname)
					
#		tinsertions+=insertions
#		tinslength+=inslength
#		tdeletions+=deletions
#		tdellength+=dellength
#		tSNPs+=SNPs
#		tclipped+=clipped
#		tcliplength+=cliplength
#		tmapped+=mapped
		
		readinfo[readname]=insertions+deletions+SNPs+clipped
		read_name_dict[readname]=read.qname
#		print
#		print readname, start, end, insertions, deletions, SNPs, clipped
	
	samfile.close()
	#print len(overlapping_reads), len(readinfo), count, count2
	#print read_name_dict
	return readinfo, overlapping_reads, read_name_dict
		
		


def map_reads(freads="", rreads="", ref=""):
	
	if freads=="":
		print "No reads given"
		sys.exit()
	if ref=="":
		print "No reference given"
		sys.exit()
	
	#index the reference
	os.system(SMALT_DIR+"smalt index -k 13 -s 1 "+ref+".index "+ref+" > /dev/null 2>&1")
	#map the reads to the reference
	os.system(SMALT_DIR+"smalt map -y 0.9 -f samsoft -o "+tmpname+".sam "+ref+".index "+freads+" "+rreads+" > /dev/null 2>&1")
	os.system("samtools view -b -S -T "+ref+" -o "+tmpname+".1.bam "+tmpname+".sam > /dev/null 2>&1")
	os.system("samtools sort "+tmpname+".1.bam "+tmpname+" > /dev/null 2>&1")
	os.system("samtools index "+tmpname+".bam > /dev/null 2>&1")





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
		if samread.is_read1:
			samname="@"+samread.qname+"_1"
		elif samread.is_read2:
			samname="@"+samread.qname+"_2"
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



def create_fastq_from_bam(bamfile, fastqfile):
	
	try: samfile = pysam.Samfile( bamfile, "rb" )
	except StandardError:
		print bamfile+" not a bam file"
		sys.exit()
	
	fastqout1=open(fastqfile+"_1.fastq", "w")
	fastqout2=open(fastqfile+"_2.fastq", "w")
	foundreads={}
	for read in samfile:
			
		#print read.qname, read.is_read1, read.is_read2
		if read.qname in foundreads:
			print_read_to_file(fastqout1, read, "pairedfastq", sammate=foundreads[read.qname], mateout=fastqout2)
		else:
			foundreads[read.qname]=read
			
	
	fastqout1.close()
	fastqout2.close()
	samfile.close()


def separate_reads_to_map(readstomap, fastq_file, bam_file, outputunmapped=False, outputothers=False):
	
	
	try: samfile = pysam.Samfile( options.bam, "rb" )
	except StandardError:
		print tmpname+".bam not a bam file"
		sys.exit() 
	
	refs=samfile.references
	lengths=samfile.lengths
	
	overlapping_reads=set([])

	foundreads={}
	
	out=open(fastq_file+"_1.fastq", "w")
	mate=open(fastq_file+"_2.fastq", "w")
	if outputothers:
		sam=pysam.Samfile(bam_file, mode='wb', referencenames=refs, referencelengths=lengths)
	
	for read in samfile:
		if read.qname in readstomap or (outputunmapped and ( read.is_unmapped or read.mate_is_unmapped)):
			
			#print read.qname, read.is_read1, read.is_read2
			if read.qname in foundreads:
				print_read_to_file(out, read, "pairedfastq", sammate=foundreads[read.qname], mateout=mate)
			else:
				foundreads[read.qname]=read
			
		elif outputothers:
			#continue
			sam.write(read)

	out.close()
	mate.close()
	samfile.close()
	if outputothers:
		sam.close()


########
# Main #
########


if __name__ == "__main__":



	(options, args)=get_user_options()
	#print options, args
	reference_contig_names, reference_contig_lengths=check_input_validity(options, args)
	
	#make random name for files
	chars = string.ascii_letters + string.digits
	if options.tmpname=="":
		tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))
	else:
		tmpname=options.tmpname
	
	if options.wd!="":
		tmpname=options.wd+"/"+tmpname
	
	if not os.path.isfile(options.bam+".bai"):
			print >> bashfile, SAMTOOLS_DIR+"samtools index", options.bam
	
	
	print "Reading reference contigs"
	
	sequences={}
	altsequences={}
	try:
		reffile=open(options.ref, "rU")
	except StandardError:
		print "Could not open reference sequence file"
		sys.exit()
	lines=reffile.read().split('>')[1:]
	for line in lines:
		words=line.strip().split('\n')
		sequences[words[0].split()[0]]=''.join(words[1:]).upper()
	reffile.close()
	
	
	
	
	
	ISsequences={}
	if options.ISfasta!="" and os.path.isfile(options.ISfasta):
		print "Reading IS element sequences"
		try:
			ISfile=open(options.ISfasta, "rU")
		except StandardError:
			print "Could not open IS sequence file"
			sys.exit()
		lines=ISfile.read().split('>')[1:]
		for line in lines:
			words=line.strip().split('\n')
			ISsequences[words[0].split()[0]]=''.join(words[1:]).upper()
		ISfile.close()
	
	#get the maximum read length
	
	print "Extracting read mapping information"
	sys.stdout.flush()
	
	try: samfile = pysam.Samfile( options.bam, "rb" )
	except StandardError:
		print tmpname+".bam not a bam file"
		sys.exit() 
		
	maxrlen=0
	insertsizes=[]
		
	for read in samfile:
		if read.rlen>maxrlen:
			maxrlen=read.rlen
		if read.tlen>0 and read.tlen<1000:
			insertsizes.append(read.tlen)
		elif read.tlen<0 and read.tlen>-1000:
			insertsizes.append(read.tlen*-1)
	
	meaninsert=mean(insertsizes)
	medianinsert=median(insertsizes)
	stdinsert=std(insertsizes)
	print "Max read length =", maxrlen
	print "Mean insert size =", meaninsert
	print "Median insert size =", medianinsert
	print "Insert size standard deviation =", stdinsert
	
	depths=[]
	
	for pileupcolumn in samfile.pileup():
		depths.append(pileupcolumn.n)
		
	meandepth=mean(depths)
	mediandepth=median(depths)
	stddepth=std(depths)
	
	print "Mean read depth =", meandepth
	print "Median read depth =", mediandepth
	print "Read depth standard deviation =", stddepth	
		
	samfile.close()
	
	
	sys.stdout.flush()
	
	
	indels={}
	indelnum=0
	
	
	
	
	
	
	if options.iselements!="" and len(ISsequences)>0:
		try: samfile = pysam.Samfile( options.iselements, "rb" )
		except StandardError:
			print bamfile+" not a bam file"
			sys.exit()
		coverage={}
		for read in samfile:
			
			ISname=read.qname.split("_")[-2]
			direction=read.qname.split("_")[-1]
			if not read.tid in coverage:
				coverage[read.tid]={}
			
			toadd=[]
			
			pos=read.pos
			for cig in read.cigar:
				if cig[0]==0:
					for x in xrange(0,cig[1]):
						toadd.append(pos)
						pos+=1
				elif cig[0]==2:
					for x in xrange(0,cig[1]):
						pos+=1
			#print read.pos, read.cigar, toadd	
			
			if read.is_reverse and direction=="R":
				status="forward"
				for pos in toadd:
					if not pos in coverage[read.tid]:
						coverage[read.tid][pos]={}
					if not ISname in coverage[read.tid][pos]:
						coverage[read.tid][pos][ISname]={}
					if not "elementforward" in coverage[read.tid][pos][ISname]:
						coverage[read.tid][pos][ISname]["elementforward"]={}
					if not "readreverse" in coverage[read.tid][pos][ISname]["elementforward"]:
						coverage[read.tid][pos][ISname]["elementforward"]["readreverse"]=0
					coverage[read.tid][pos][ISname]["elementforward"]["readreverse"]+=1
			elif read.is_reverse and direction=="F":
				status="reverse"
				for pos in toadd:
					if not pos in coverage[read.tid]:
						coverage[read.tid][pos]={}
					if not ISname in coverage[read.tid][pos]:
						coverage[read.tid][pos][ISname]={}
					if not "elementreverse" in coverage[read.tid][pos][ISname]:
						coverage[read.tid][pos][ISname]["elementreverse"]={}
					if not "readreverse" in coverage[read.tid][pos][ISname]["elementreverse"]:
						coverage[read.tid][pos][ISname]["elementreverse"]["readreverse"]=0
					coverage[read.tid][pos][ISname]["elementreverse"]["readreverse"]+=1
			elif not read.is_reverse and direction=="F":
				status="forward"
				for pos in toadd:
					if not pos in coverage[read.tid]:
						coverage[read.tid][pos]={}
					if not ISname in coverage[read.tid][pos]:
						coverage[read.tid][pos][ISname]={}
					if not "elementforward" in coverage[read.tid][pos][ISname]:
						coverage[read.tid][pos][ISname]["elementforward"]={}
					if not "readforward" in coverage[read.tid][pos][ISname]["elementforward"]:
						coverage[read.tid][pos][ISname]["elementforward"]["readforward"]=0
					coverage[read.tid][pos][ISname]["elementforward"]["readforward"]+=1
					coverage[read.tid][pos][ISname]["elementforward"]["readforward"]+=1
			elif not read.is_reverse and direction=="R":
				status="reverse"
				for pos in toadd:
					if not pos in coverage[read.tid]:
						coverage[read.tid][pos]={}
					if not ISname in coverage[read.tid][pos]:
						coverage[read.tid][pos][ISname]={}
					if not "elementreverse" in coverage[read.tid][pos][ISname]:
						coverage[read.tid][pos][ISname]["elementreverse"]={}
					if not "readforward" in coverage[read.tid][pos][ISname]["elementreverse"]:
						coverage[read.tid][pos][ISname]["elementreverse"]["readforward"]=0
					coverage[read.tid][pos][ISname]["elementreverse"]["readforward"]+=1
					coverage[read.tid][pos][ISname]["elementreverse"]["readforward"]+=1
			#print ISname, direction, read.is_reverse, status
	
		tidlist=coverage.keys()
		for tid in tidlist:
			poslist=coverage[tid].keys()
			for pos in poslist:
				elementlist=coverage[tid][pos].keys()
				for element in elementlist:
					eldrnlist=coverage[tid][pos][element].keys()
					for eldrn in eldrnlist:
						#print coverage[tid][pos][element][eldrn].keys()
						if not "readforward" in  coverage[tid][pos][element][eldrn] or not "readreverse" in coverage[tid][pos][element][eldrn]:
							del coverage[tid][pos][element][eldrn]
					if len(coverage[tid][pos][element].keys())==0:
						del coverage[tid][pos][element]
				if len(coverage[tid][pos].keys())==0:
					del coverage[tid][pos]
			if len(coverage[tid].keys())==0:
				del coverage[tid]
		
		overlaps=[]
		for tid in coverage:
			for pos in coverage[tid]:
				for element in coverage[tid][pos]:
					for eldrn in coverage[tid][pos][element]:
						if not "readforward" in  coverage[tid][pos][element][eldrn] or not "readreverse" in coverage[tid][pos][element][eldrn]:
							continue
						else:
							count=coverage[tid][pos][element][eldrn]["readforward"]+coverage[tid][pos][element][eldrn]["readreverse"]
							overlaps.append([element, eldrn, samfile.getrname(tid), pos, count])
						
		overlaps.sort()
		possibleinserts=[]
		currinsert=[]
		prevoverlap=[]
		for overlap in overlaps:
			if prevoverlap==[]:
				prevoverlap=overlap
				continue
			if not (overlap[0]==prevoverlap[0] and overlap[3]==prevoverlap[3]+1 and overlap[2]==prevoverlap[2] and overlap[1]==prevoverlap[1]):
				possibleinserts.append(currinsert)
				currinsert=[]
				
			currinsert.append(overlap)
			prevoverlap=overlap
		possibleinserts.append(currinsert)
		
		for possibleinsert in  possibleinserts:
			print possibleinsert
			indelnum+=1
			indels[indelnum]={}
			indels[indelnum]["info"]={}
			indels[indelnum]["chrom"]=possibleinsert[0][2]
			indels[indelnum]['start']=possibleinsert[-1][3]+1
			indels[indelnum]['info']['END']=possibleinsert[-1][3]+2
			indels[indelnum]["info"]["SVTYPE"]="INS"
			indels[indelnum]["info"]["info"]="from ISscan: "+possibleinsert[0][0]
			indels[indelnum]['ref']=sequences[indels[indelnum]["chrom"]][indels[indelnum]['start']-1]
			if possibleinsert[0][1]=="elementforward":
				indels[indelnum]['alt']=sequences[indels[indelnum]["chrom"]][possibleinsert[-1][3]-1]+ISsequences[possibleinsert[0][0]]+sequences[indels[indelnum]["chrom"]][possibleinsert[0][3]-1:possibleinsert[-1][3]+1]
			else:
				indels[indelnum]['alt']=sequences[indels[indelnum]["chrom"]][possibleinsert[-1][3]]+revcomp(ISsequences[possibleinsert[0][0]])+sequences[indels[indelnum]["chrom"]][possibleinsert[0][3]-1:possibleinsert[-1][3]+1]
			indels[indelnum]["info"]["SVLEN"]=len(indels[indelnum]['alt'])-len(indels[indelnum]['ref'])
				
			
			
			
			print indels[indelnum]
		
#		indels[indelnum]['alt']=BASEINFO["ALT"].split(",")[0][0].upper()+altseq
#		indels[indelnum]['ref']=BASEINFO["REF"].split(",")[0][0].upper()+refseq
#		indels[indelnum]['start']=int(BASEINFO["POS"])
#		indels[indelnum]['info']['END']=int(BASEINFO["POS"]+len(indels[indelnum]['ref']))
#		indels[indelnum]["info"]["SVLEN"]=len(indels[indelnum]['alt'])-len(indels[indelnum]['ref'])
#		if indels[indelnum]["info"]["SVLEN"]>0:
#			indels[indelnum]["info"]["SVTYPE"]="INS"
#		else:
#			indels[indelnum]["info"]["SVTYPE"]="DEL"
#		indels[indelnum]["chrom"]=BASEINFO["CHROM"]
#		indels[indelnum]["info"]["info"]="from bcf"
	
	#sys.exit()
	
	
	
	
	
	
	if options.bcf!="":
		if options.vcf:
			try:
				bcffile=open(options.bcf, "rU")
			except StandardError:
				DoError("Cannot open vcf file")
		else:
			try:
				bcffile=os.popen(BCFTOOLS_DIR+"bcftools view "+options.bcf)
			except StandardError:
				DoError("Cannot open bcf file")
	
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
		
			if 'INDEL' in BASEINFO['INFO'] and BASEINFO["ALT"].split(",")[0]!=".":
				indelnum+=1
				indels[indelnum]={}
				indels[indelnum]["info"]={}
				
				altseq=BASEINFO["ALT"].split(",")[0][1:].upper()
				refseq=BASEINFO["REF"].split(",")[0][1:].upper()
				
				if len(altseq)>len(refseq) and len(refseq)>0:
					if altseq[len(refseq)*-1:]==refseq:
						altseq=altseq[:len(refseq)*-1]
						refseq=refseq[:len(refseq)*-1]
				elif len(refseq)>len(altseq) and len(altseq)>0:
					if refseq[len(altseq)*-1:]==altseq:
						refseq=refseq[:len(altseq)*-1]
						altseq=altseq[:len(altseq)*-1]
				
				indels[indelnum]['alt']=BASEINFO["ALT"].split(",")[0][0].upper()+altseq
				indels[indelnum]['ref']=BASEINFO["REF"].split(",")[0][0].upper()+refseq
				indels[indelnum]['start']=int(BASEINFO["POS"])
				indels[indelnum]['info']['END']=int(BASEINFO["POS"]+len(indels[indelnum]['ref']))
				indels[indelnum]["info"]["SVLEN"]=len(indels[indelnum]['alt'])-len(indels[indelnum]['ref'])
				if indels[indelnum]["info"]["SVLEN"]>0:
					indels[indelnum]["info"]["SVTYPE"]="INS"
				else:
					indels[indelnum]["info"]["SVTYPE"]="DEL"
				indels[indelnum]["chrom"]=BASEINFO["CHROM"]
				indels[indelnum]["info"]["info"]="from bcf"
				
				#print indels[indelnum]['ref'], indels[indelnum]['alt'], indels[indelnum]['start'], indels[indelnum]['info']['END']
				# the rest of the loop moves indels back to the start of repeats
				
				x=indels[indelnum]['start']
				if indels[indelnum]["info"]["SVLEN"]>0:
					indellength=indels[indelnum]["info"]["SVLEN"]
				else:
					indellength=indels[indelnum]["info"]["SVLEN"]*-1
				if indels[indelnum]["info"]["SVTYPE"]=="INS":
					alt=indels[indelnum]["alt"][indellength*-1:]
				elif indels[indelnum]["info"]["SVTYPE"]=="DEL":
					alt=indels[indelnum]["ref"][indellength*-1:]
				
				
				
				x-=indellength
				while sequences[indels[indelnum]["chrom"]][x:x+indellength]==alt:
					x-=indellength
				indels[indelnum]['start']=x+indellength
				
				if indels[indelnum]["info"]["SVTYPE"]=="INS":
					indels[indelnum]['alt']=sequences[indels[indelnum]["chrom"]][x+indellength-1:x+indellength]+alt
					indels[indelnum]['ref']=sequences[indels[indelnum]["chrom"]][x+indellength-1]
				elif indels[indelnum]["info"]["SVTYPE"]=="DEL":
					indels[indelnum]['ref']=sequences[indels[indelnum]["chrom"]][x+indellength-1:x+indellength]+alt
					indels[indelnum]['alt']=sequences[indels[indelnum]["chrom"]][x+indellength-1]
					
				indels[indelnum]['info']['END']=indels[indelnum]['start']+len(indels[indelnum]['ref'])
				#print indels[indelnum]['ref'], indels[indelnum]['alt'], indels[indelnum]['start'], indels[indelnum]['info']['END']
				
		bcffile.close()
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	#print len(indels)
	
	#run pindel
	
	if options.runpindel:
		print "Running pindel"
		sys.stdout.flush()
		os.system(' '.join(map(str, [SAMTOOLS_DIR+"samtools view ",  options.bam, "|", PINDEL_DIR+"sam2pindel -", tmpname+".pindelin "+str(meaninsert)+" "+tmpname+" 0 Illumina-PairEnd"])))
		os.system(' '.join(map(str, [PINDEL_DIR+"pindel -u 0.2 -e 0.05 -m 2 -x 6 -f "+options.ref+" -p", tmpname+".pindelin -c ALL -o", tmpname+".pindel.out"])))
		os.system(' '.join(map(str, ["cat", tmpname+".pindel.out_D", tmpname+".pindel.out_LI", tmpname+".pindel.out_SI >", tmpname+".pindel.out"])))
		os.system(' '.join(map(str, [PINDEL_DIR+"pindel2vcf -p", tmpname+".pindel.out -r", options.ref, "-R", options.ref.split("/")[-1].split(".")[0], "-d", ''.join(map(str,[1979, 1, 30]))])))
		os.system('rm -f '+tmpname+'.pindelin')
	
	
	#read in the pindel output file
	print "Reading pindel output"
	sys.stdout.flush()
	
	
	pindelout=open(tmpname+".pindel.out.vcf", "rU")
	for line in pindelout:
		if len(line)>0 and line[0]=="#":
			continue
		words=line.strip().split()
		indelnum+=1
		indels[indelnum]={}
				
		indels[indelnum]["chrom"]=words[0]
		indels[indelnum]["start"]=int(words[1])
		indels[indelnum]["ref"]=words[3]
		indels[indelnum]["alt"]=words[4]
		indels[indelnum]["info"]={}
		indels[indelnum]["info"]["info"]="from pindel"
		for i in words[7].split(";"):
			bits=i.split("=")
			if bits[0]=="END":
				indels[indelnum]["info"][bits[0]]=int(bits[1])+1
			elif bits[0]=="SVLEN":
				indels[indelnum]["info"][bits[0]]=int(bits[1])
			else:
				indels[indelnum]["info"][bits[0]]=bits[1]
		if words[4]=="<INS>":
			indels[indelnum]["info"]["SVTYPE"]="LI"
	pindelout.close()
	
	#print len(indels)
	
	indelpositions=[]
	
	for indel in indels:
		indelpositions.append([int(indels[indel]["info"]["END"]), indel])
	
	
	
	
	
	
	removedindels=[]
	
	#for now I'll only look at insertions, deletions and make an attempt at replacements
	toremove=set([])
	toadd=[]
	
	print "Creating proposed indels from replacements"
	sys.stdout.flush()
	
	for indel in indels:
		if indels[indel]["info"]["SVTYPE"] not in ["INS","DEL", "RPL"]:
			removedindels.append(indels[indel])
			toremove.add(indel)
		elif indels[indel]["info"]["SVTYPE"] in ["RPL"]:
			toremove.add(indel)
			seqs={}
#			print len(indels[indel]["ref"]), len(indels[indel]["alt"])
			if len(indels[indel]["ref"])>500 and len(indels[indel]["ref"])>len(indels[indel]["alt"])*2:
				refbitlen=int(len(indels[indel]["ref"]))/2
				if refbitlen>len(indels[indel]["alt"])*4:
					refbitlen=len(indels[indel]["alt"])*4
				refstart=indels[indel]["ref"][:refbitlen]
				refend=indels[indel]["ref"][-1*refbitlen:]
				#print "\n",indels[indel]["alt"]
				#print refstart
				#print refend
				scores=[]
				alns=[]
				for refseq in [refstart, refend]:
					fastout=open(tmpname+".fasta", "w")
					print >> fastout, ">alt"
					print >> fastout, indels[indel]["alt"]
					print >> fastout, ">ref"
					print >> fastout, refseq
					fastout.close()
					
					os.system("muscle -in "+tmpname+".fasta -out "+tmpname+".aln > /dev/null 2>&1")
					
					
					muscleout=open(tmpname+".aln", "rU")
				
					curseq=""
					for line in muscleout:
						words=line.strip().split()
						if len(words)>0 and len(words[0])>0:
							if words[0][0]==">":
								curseq=words[0][1:]
								seqs[curseq]=""
							elif curseq in seqs:
								seqs[curseq]+="".join(words)
					alns.append(seqs)
					score=0
					if "alt" in seqs and "ref" in seqs:
#						print seqs["ref"]
#						print seqs["alt"]
						for x in range(0, len(seqs["ref"])):
							if seqs["ref"][x]==seqs["alt"][x]:
								score+=1
					else:
						print seqs
						sys.exit()
					#print score, len(indels[indel]["alt"])-score
					if len(indels[indel]["alt"])>len(indels[indel]["ref"]):
						#print score, len(indels[indel]["ref"]), float(score)/len(indels[indel]["ref"])
						scores.append(float(score)/len(indels[indel]["ref"]))
					else:
						#print score, len(indels[indel]["alt"]), float(score)/len(indels[indel]["alt"])
						scores.append(float(score)/len(indels[indel]["alt"]))
				#print scores
				if scores[0]>scores[1]:
					seqs={}
					seqs["alt"]=alns[0]["alt"]+(len(indels[indel]["ref"])-refbitlen)*"-"
					seqs["ref"]=alns[0]["ref"]+indels[indel]["ref"][:len(indels[indel]["ref"])-refbitlen]
#					print scores[0]
				else:
					seqs={}
					seqs["alt"]=(len(indels[indel]["ref"])-refbitlen)*"-"+alns[1]["alt"]
					seqs["ref"]=indels[indel]["ref"][:len(indels[indel]["ref"])-refbitlen]+alns[1]["ref"]
					
			elif len(indels[indel]["alt"])==0:
				seqs["ref"]=indels[indel]["ref"]
				seqs["alt"]=len(indels[indel]["ref"])*"-"
			elif len(indels[indel]["ref"])==0:
				seqs["alt"]=indels[indel]["alt"]
				seqs["ref"]=len(indels[indel]["alt"])*"-"
			else:
				fastout=open(tmpname+".fasta", "w")
				print >> fastout, ">alt"
				print >> fastout, indels[indel]["alt"]
				print >> fastout, ">ref"
				print >> fastout, indels[indel]["ref"]
				fastout.close()
				
				os.system("muscle -in "+tmpname+".fasta -out "+tmpname+".aln > /dev/null 2>&1")
				
				muscleout=open(tmpname+".aln", "rU")
				
				curseq=""
				for line in muscleout:
					words=line.strip().split()
					if len(words)>0 and len(words[0])>0:
						if words[0][0]==">":
							curseq=words[0][1:]
							seqs[curseq]=""
						elif curseq in seqs:
							seqs[curseq]+="".join(words)
				score=0
				if "alt" in seqs and "ref" in seqs:
					for x in range(0, len(seqs["ref"])):
						if seqs["ref"][x]==seqs["alt"][x]:
							score+=1
				else:
					print seqs
					sys.exit()
#			print
#			if len(seqs["ref"])<550:
#				print seqs["ref"]
#				print seqs["alt"]
#				print indels[indel]["start"]
			ininsertion=False
			indeletion=False
			insertionstart=0
			deletionstart=0
			for x in range(0, len(seqs["ref"])):
				if seqs["ref"][x]=="-" and not ininsertion:
					ininsertion=True
					insertionstart=x
				elif seqs["ref"][x]!="-" and ininsertion:
					ininsertion=False
#					print "INS", insertionstart, x
					
					newindel={}
					newindel["chrom"]=indels[indel]["chrom"]
					newindel["start"]=insertionstart-1+indels[indel]["start"]
					newindel["info"]={}
					newindel["info"]["SVTYPE"]="INS"
					if 'info' in indels[indel]["info"]:
						newindel["info"]["info"]=indels[indel]["info"]["info"]+" replacement"
					else:
						newindel["info"]["info"]="from replacement"
					newindel["info"]["END"]=newindel["start"]+1
					refindelpos=(indels[indel]["start"]+len(seqs["ref"][:insertionstart].replace("-","")))-2

					newindel["ref"]=sequences[newindel["chrom"]][refindelpos]+seqs["ref"][insertionstart:x].replace("-","")
					newindel["alt"]=sequences[newindel["chrom"]][refindelpos]+seqs["alt"][insertionstart:x].replace("-","")
#					if len(newindel["ref"])<550 and len(newindel["alt"])<550:
#						print newindel["ref"]
#						print newindel["alt"]
#						print newindel["start"], len(newindel["alt"])-len(newindel["ref"]), x, newindel["info"]["END"]
#						print indels[indel]["start"]
					newindel["info"]["SVLEN"]=len(newindel["alt"])-len(newindel["ref"])
					toadd.append(newindel)
							
				if seqs["alt"][x]=="-" and not indeletion:
					indeletion=True
					deletionstart=x
				elif seqs["alt"][x]!="-" and indeletion:
					indeletion=False
#					print "DEL", deletionstart, x
					newindel={}
					newindel["chrom"]=indels[indel]["chrom"]
					newindel["start"]=deletionstart+indels[indel]["start"]-1
					newindel["info"]={}
					newindel["info"]["SVTYPE"]="DEL"
					if 'info' in indels[indel]["info"]:
						newindel["info"]["info"]=indels[indel]["info"]["info"]+" replacement"
					else:
						newindel["info"]["info"]="from replacement"
					newindel["info"]["END"]=x+indels[indel]["start"]
					refindelpos=(indels[indel]["start"]+len(seqs["ref"][:deletionstart].replace("-","")))-2

					newindel["ref"]=sequences[newindel["chrom"]][refindelpos]+seqs["ref"][deletionstart:x].replace("-","")
					newindel["alt"]=sequences[newindel["chrom"]][refindelpos]+seqs["alt"][deletionstart:x].replace("-","")
#					if len(newindel["ref"])<550 and len(newindel["alt"])<550:
#						print newindel["ref"]
#						print newindel["alt"]
#						print newindel["start"], len(newindel["alt"])-len(newindel["ref"]), x, newindel["info"]["END"]
#						print indels[indel]["start"]
					newindel["info"]["SVLEN"]=len(newindel["alt"])-len(newindel["ref"])
					toadd.append(newindel)
			
			if indeletion:
				indeletion=False
				x+=1
#				print "DEL", deletionstart, x
				newindel={}
				newindel["chrom"]=indels[indel]["chrom"]
				newindel["start"]=deletionstart+indels[indel]["start"]-1
				newindel["info"]={}
				newindel["info"]["SVTYPE"]="DEL"
				if 'info' in indels[indel]["info"]:
					newindel["info"]["info"]=indels[indel]["info"]["info"]+" last of replacement"
				else:
					newindel["info"]["info"]="from replacement"
				newindel["info"]["END"]=x+indels[indel]["start"]
				refindelpos=(indels[indel]["start"]+len(seqs["ref"][:deletionstart].replace("-","")))-2

				newindel["ref"]=sequences[newindel["chrom"]][refindelpos]+seqs["ref"][deletionstart:x+1].replace("-","")
				newindel["alt"]=sequences[newindel["chrom"]][refindelpos]+seqs["alt"][deletionstart:x+1].replace("-","")
#				if len(newindel["ref"])<550 and len(newindel["alt"])<550:
#					print newindel["ref"]
#					print newindel["alt"]
#					print newindel["start"], len(newindel["alt"])-len(newindel["ref"]), x, newindel["info"]["END"]
#					print indels[indel]["start"]
				newindel["info"]["SVLEN"]=len(newindel["alt"])-len(newindel["ref"])
				toadd.append(newindel)
			if ininsertion:
				ininsertion=False
				x+=1
#				print "INS", insertionstart, x
				newindel={}
				newindel["chrom"]=indels[indel]["chrom"]
				newindel["start"]=insertionstart-1+indels[indel]["start"]
				newindel["info"]={}
				newindel["info"]["SVTYPE"]="INS"
				if 'info' in indels[indel]["info"]:
					newindel["info"]["info"]=indels[indel]["info"]["info"]+" last of replacement"
				else:
					newindel["info"]["info"]="from replacement"
				newindel["info"]["END"]=newindel["start"]+1
				refindelpos=(indels[indel]["start"]+len(seqs["ref"][:insertionstart].replace("-","")))-2
				newindel["ref"]=sequences[newindel["chrom"]][refindelpos]+seqs["ref"][insertionstart:x+1].replace("-","")
				newindel["alt"]=sequences[newindel["chrom"]][refindelpos]+seqs["alt"][insertionstart:x+1].replace("-","")
#				if len(newindel["ref"])<550 and len(newindel["alt"])<550:
#					print newindel["ref"]
#					print newindel["alt"]
#					print newindel["start"], len(newindel["alt"])-len(newindel["ref"]), x, newindel["info"]["END"]
#					print indels[indel]["start"]
				newindel["info"]["SVLEN"]=len(newindel["alt"])-len(newindel["ref"])
				toadd.append(newindel)
			
	
	
	for indel in toremove:
		del indels[indel]

	#indels={}
	for x, indel in enumerate(toadd):
		#if indel["start"]==1254048:
		#if 'info' in indel["info"] and indel["info"]["info"]==" from pindel replacement":
		indelnum+=1
		indels[indelnum]=indel


	indelpositions=[]
	
	for indel in indels:
		indelpositions.append([int(indels[indel]["info"]["END"]), indel])
	
	#print "Before removing identicals", len(indels)
	
	
	indelpositions.sort()
	indelpositions.reverse()
	#sort the pindel indels
	
	
	
	print "Removing identical indels"
	sys.stdout.flush()
	toremove=set([])
	for i, indelpos in enumerate(indelpositions):#[:5]:
		
		if indelpos[1] in toremove:
			continue
		
		j=i+1
		
		while j<len(indelpositions) and indelpos[0]==indelpositions[j][0]:
			
			if indels[indelpositions[i][1]]["chrom"]==indels[indelpositions[j][1]]["chrom"] and indels[indelpositions[i][1]]["start"]==indels[indelpositions[j][1]]["start"] and indels[indelpositions[i][1]]["alt"].upper()==indels[indelpositions[j][1]]["alt"].upper():
				toremove.add(indelpositions[j][1])
				
			j+=1
	

	
	for indel in toremove:
		del indels[indel]


	indelpositions=[]
	
	for indel in indels:
		indelpositions.append([int(indels[indel]["info"]["END"]), indel])
	
#	for indel in indelpositions:
#		#if (indels[indel[1]]["info"]["SVLEN"]<50 and indels[indel[1]]["info"]["SVLEN"]>0) or (indels[indel[1]]["info"]["SVLEN"]>-50 and indels[indel[1]]["info"]["SVLEN"]<0):
#		if not 'info' in indels[indel[1]]["info"] or (indels[indel[1]]["info"]["info"]!="from pindel replacement" and indels[indel[1]]["info"]["info"]!="from pindel last of replacement"):
#			del indels[indel[1]]
#			continue
##		if indel[0]<1250000 or indel[0]>1260000:
##			del indels[indel[1]]
#	
#	indelpositions=[]
#	
#	for indel in indels:
#		indelpositions.append([int(indels[indel]["info"]["END"]), indel])
#	
	print len(indels), "indels proposed"
	sys.stdout.flush()
	
	indelpositions.sort()
	#indelpositions=indelpositions[3:5]
	indelpositions.reverse()
	#indelpositions=indelpositions[78:80]

	#read in the reference genome sequence
	
	altsequences={}	
	



	# print a tab file of indel locations
	toremove=set([])
	readstomap=set([])
	
	
	
	ref_contigs={}
	#Write the reference to file
	output=open(tmpname+".ref.dna", "w")
	for sequence in sequences:
		print >> output, ">"+sequence
		print >> output, sequences[sequence]
		ref_contigs[sequence]=sequences[sequence]
	output.close()
	
	
	#for each indel
	print "Assessing proposed indels"
	sys.stdout.flush()
	
	for indelpos in indelpositions:
		indel=indels[indelpos[1]]
#		if not indel["info"]["info"]=="from ISscan: is6110":
#			toremove.add(indelpos[1])
#			continue
		#first extract reads around the indel
		if (int(indel["info"]["SVLEN"])>0 and int(indel["info"]["SVLEN"])<maxrlen) or int(indel["info"]["SVLEN"])*-1<maxrlen:
			os.system(SAMTOOLS_DIR+"samtools view -b -o "+tmpname+".bam "+options.bam+" '"+indel["chrom"]+":"+str(int(indel["start"])-100)+"-"+str(int(indel["info"]["END"])+100)+"'")
		else:
			os.system(SAMTOOLS_DIR+"samtools view -b -o "+tmpname+".bam "+options.bam+" '"+indel["chrom"]+":"+str(int(indel["start"])-100)+"-"+str(int(indel["start"])+100)+"' '"+indel["chrom"]+":"+str(int(indel["info"]["END"])-100)+"-"+str(int(indel["info"]["END"])+100)+"'")
		
		
		regionreads=set([])
		
		try: samfile = pysam.Samfile( tmpname+".bam", "rb" )
		except StandardError:
			print tmpname+".bam not a bam file"
			sys.exit() 
		
		for read in samfile:
			regionreads.add(read.qname)
		samfile.close()
		
#		print len(regionreads)
		
		if len(regionreads)<options.depth:
			toremove.add(indelpos[1])
			continue
		
		separate_reads_to_map(regionreads, tmpname, "")
		
		
	
		#map the reads
		map_reads(freads=tmpname+"_1.fastq", rreads=tmpname+"_2.fastq", ref=tmpname+".ref.dna")
		
		results=samfile_data(int(indel["start"]), int(indel["info"]["END"]), contig_seqs=ref_contigs)
		
		refreadinfo=results[0]
		refset=results[1]
		refdict=results[2]
		#print len(results[1])
		
		#write the alternate sequences
		contigs={}
		output=open(tmpname+".alt.dna", "w")
		for sequence in reference_contig_names:
			if sequence==indel["chrom"]:
				print >> output, ">"+sequence
				print >> output, sequences[sequence][:indel["start"]-1]+indel["alt"]+sequences[sequence][int(indel["info"]["END"])-1:]
				contigs[sequence]=sequences[sequence][:indel["start"]-1]+indel["alt"]+sequences[sequence][int(indel["info"]["END"])-1:]
			else:
				print >> output, ">"+sequence
				print >> output, sequences[sequence]
				contigs[sequence]=sequences[sequence]
		output.close()
		
		
		
		
		
		
		
#		sys.exit()
		#map the reads
		map_reads(freads=tmpname+"_1.fastq", rreads=tmpname+"_2.fastq", ref=tmpname+".alt.dna")
		
		results=samfile_data(int(indel["start"]), int(indel["info"]["END"])+int(indel["info"]["SVLEN"]), contig_seqs=contigs)
		
		altreadinfo=results[0]
		#print len(results[1])
		altset=results[1]
		altdict=results[2]
		
		alloverlaps=altset.union(refset)
		
		awins=0
		rwins=0
		draws=0
		awinset=set([])
		drawset=set([])
		
		for read in alloverlaps:
#			print refreadinfo[read], altreadinfo[read]
			if refreadinfo[read]<altreadinfo[read]:
				rwins+=1
			elif refreadinfo[read]>altreadinfo[read]:
				awins+=1
				awinset.add(altdict[read])
			else:
				draws+=1
				drawset.add(altdict[read])
		
		indel["awins"]=awins
		indel["rwins"]=rwins
		indel["draws"]=draws
		
		#print indel["info"]["SVTYPE"], indel["start"], indel["info"]["END"], indel["ref_score"], indel["alt_score"], indel["rwins"], indel["awins"], indel["draws"]
		if "info" in indel["info"]:
			print indel["chrom"], indel["info"]["SVLEN"], indel["info"]["SVTYPE"], indel["start"], indel["info"]["END"], indel["rwins"], indel["awins"], indel["draws"], indel["info"]["info"],
		else:
			print indel["chrom"], indel["info"]["SVLEN"], indel["info"]["SVTYPE"], indel["start"], indel["info"]["END"], indel["rwins"], indel["awins"], indel["draws"], #, len(alloverlaps)
		
		
		if indel.has_key("awins"):
			#print indel["awins"]>indel["rwins"], float(indel["awins"]), (indel["awins"]+indel["rwins"]), options.proportion, indel["awins"]+indel["rwins"]+indel["draws"], options.depth,
			if indel["awins"]>indel["rwins"] and (options.proportion==0 or (float(indel["awins"])/(indel["awins"]+indel["rwins"]))>=options.proportion) and (indel["awins"]+indel["rwins"]+indel["draws"]>options.depth):
#				if indel["info"]["SVTYPE"]=="DEL" and indel["info"]["END"]-indel["start"]>10:
#					#print "del", 
#					
#					#os.system(SAMTOOLS_DIR+"samtools view -b -o "+tmpname+".bam "+options.bam+" '"+indel["chrom"]+":"+str(int(indel["start"]))+"-"+str(int(indel["info"]["END"]))+"'")
#		
#					
#					try: samfile = pysam.Samfile( options.bam, "rb" )
#					except StandardError:
#						print tmpname+".bam not a bam file"
#						sys.exit()
#					lastcolumn=indel["start"]-1
#					cov_list=[]
#					for pileupcolumn in samfile.pileup(indel["chrom"], start=indel["start"], end=indel["info"]["END"], truncate=True):
#						while pileupcolumn.pos!=lastcolumn+1:
#							cov_list.append(0)
#							lastcolumn+=1
#						cov_list.append(pileupcolumn.n)
#						lastcolumn=pileupcolumn.pos
#					#print max(cov_list), mean(cov_list), median(cov_list),
#					
#					
##					for read in samfile:
##						regionreads.add(read.qname)
#					samfile.close()
					#sys.exit()
					
				print "...Accepted"
				readstomap.update(awinset)
				readstomap.update(drawset)
			else:
				toremove.add(indelpos[1])
				print "...Rejected"
		else:
			toremove.add(indelpos[1])
			print "...Rejected"
		sys.stdout.flush()
		
		indels[indelpos[1]]=indel
	

	#remove any indels that have been rejected
				
	for indel in toremove:
		removedindels.append(indels[indel])
		del indels[indel]
	
	#print indels
	indelpositions=[]
	
	for indel in indels:
		indelpositions.append([int(indels[indel]["info"]["END"]), indel])
	
	indelpositions.sort()
	indelpositions.reverse()
	
	
	#	tabout=open(options.output+"_b.tab", "w")
	
	toremove=set([])
	
	doreassessment=True
	
	if doreassessment:
	
		#check that overlapping indels aren't actually one indel
		print "Reassessing indels within a readlength"
		sys.stdout.flush()
		
		for i, indelpos in enumerate(indelpositions):#[:5]:
			
			if indelpos[1] in toremove:
				continue
			
			j=i+1
			
			while j<len(indelpositions) and (indelpositions[i][0]-indelpositions[j][0]<=maxrlen or (indels[indelpos[1]]["start"])-indelpositions[j][0]<=maxrlen):
			
				indel2pos=indelpositions[j]
				j+=1
				
				indel=indels[indelpos[1]]
				indel2=indels[indel2pos[1]]
				
				if indel["chrom"]!=indel2["chrom"]:
					continue
				
				minstart=int(indel2["start"])
				if int(indel["start"])<minstart:
					#toremove.add(indel2pos[1])
					continue
				
				#first extract reads around the two indels
				os.system(SAMTOOLS_DIR+"samtools view -b -o "+tmpname+".bam "+options.bam+" '"+indel["chrom"]+":"+str(minstart-100)+"-"+str(int(indel["info"]["END"])+100)+"'")
				
				regionreads=set([])
				
				try: samfile = pysam.Samfile( tmpname+".bam", "rb" )
				except StandardError:
					print tmpname+".bam not a bam file"
					sys.exit() 
				
				for read in samfile:
					regionreads.add(read.qname)
				samfile.close()
				
				separate_reads_to_map(regionreads, tmpname, "")
				
				
				contigs={}
				#Write the reference to file
				output=open(tmpname+".ref.dna", "w")
				for sequence in reference_contig_names:
					if sequence==indel["chrom"]:
						print >> output, ">"+sequence
						print >> output, sequences[sequence][:indel["start"]-1]+indel["alt"]+sequences[sequence][int(indel["info"]["END"])-1:]
						contigs[sequence]=sequences[sequence][:indel["start"]-1]+indel["alt"]+sequences[sequence][int(indel["info"]["END"])-1:]
					else:
						print >> output, ">"+sequence
						print >> output, sequences[sequence]
						contigs[sequence]=sequences[sequence]
				output.close()
			
				#map the reads
				map_reads(freads=tmpname+"_1.fastq", rreads=tmpname+"_2.fastq", ref=tmpname+".ref.dna")
				
				results=samfile_data(minstart, int(indel["info"]["END"])+int(indel["info"]["SVLEN"]), contig_seqs=contigs)
				
				indelreadinfo=results[0]
				indelset=results[1]
				indeldict=results[2]
				#print len(results[1])
				
				#write the alternate sequences
				contigs={}
				output=open(tmpname+".ref.dna", "w")
				for sequence in reference_contig_names:
					if sequence==indel["chrom"]:
						print >> output, ">"+sequence
						print >> output, sequences[sequence][:indel2["start"]-1]+indel2["alt"]+sequences[sequence][int(indel2["info"]["END"])-1:]
						contigs[sequence]=sequences[sequence][:indel2["start"]-1]+indel2["alt"]+sequences[sequence][int(indel2["info"]["END"])-1:]
					else:
						print >> output, ">"+sequence
						print >> output, sequences[sequence]
						contigs[sequence]=sequences[sequence]
				output.close()
				
				#map the reads
				map_reads(freads=tmpname+"_1.fastq", rreads=tmpname+"_2.fastq", ref=tmpname+".ref.dna")
				
				results=samfile_data(minstart, int(indel["info"]["END"])+int(indel["info"]["SVLEN"]), contig_seqs=contigs)
				
				indel2readinfo=results[0]
				#print len(results[1])
				indel2set=results[1]
				indel2dict=results[2]
				
				
				contigs={}
				for contig in reference_contig_names:
					contigs[contig]=sequences[contig]
				
				
				if int(indel2["info"]["END"])<indel["start"]:
					
					
					contigs[indel["chrom"]]=contigs[indel["chrom"]][:indel["start"]-1]+indel["alt"]+contigs[indel["chrom"]][int(indel["info"]["END"])-1:]
					contigs[indel2["chrom"]]=contigs[indel2["chrom"]][:indel2["start"]-1]+indel2["alt"]+contigs[indel2["chrom"]][int(indel2["info"]["END"])-1:]
					
							
					
					output=open(tmpname+".ref.dna", "w")
					for contig in reference_contig_names:
						print >> output, ">"+contig
						print >> output, contigs[contig]
					output.close()
					
					
					
					#map the reads
					map_reads(freads=tmpname+"_1.fastq", rreads=tmpname+"_2.fastq", ref=tmpname+".ref.dna")
					
					results=samfile_data(minstart, int(indel["info"]["END"])+int(indel["info"]["SVLEN"]), contig_seqs=contigs)
					
					bothreadinfo=results[0]
					#print len(results[1])
					bothset=results[1]
					bothdict=results[2]
					
					
					
					#compare both indels with each of the individual indels
					
					
					alloverlaps=bothset.union(indelset)
					
					biwins=0
					ibwins=0
				
					
					
					for read in alloverlaps:
						if bothreadinfo[read]<indelreadinfo[read]:
							biwins+=1
						elif bothreadinfo[read]>indelreadinfo[read]:
							ibwins+=1
					
					alloverlaps=bothset.union(indel2set)
					
					bi2wins=0
					i2bwins=0
					
					for read in alloverlaps:
						if bothreadinfo[read]<indel2readinfo[read]:
							bi2wins+=1
						elif bothreadinfo[read]>indel2readinfo[read]:
							i2bwins+=1
					
					
					
				
				else:
					ibwins=1000
					biwins=0
					i2bwins=1000
					bi2wins=0
				
				alloverlaps=indelset.union(indel2set)
				
				ii2wins=0
				i2iwins=0
				
				for read in alloverlaps:
					if indelreadinfo[read]<indel2readinfo[read]:
						ii2wins+=1
					elif indelreadinfo[read]>indel2readinfo[read]:
						i2iwins+=1
				
				print indel["chrom"], minstart, int(indel["info"]["END"]), biwins, ibwins, bi2wins, i2bwins, ii2wins, i2iwins,
				
				if ibwins>biwins and ii2wins>i2iwins:
					toremove.add(indel2pos[1])
					print "...indel 1 kept"
					sys.stdout.flush()
				elif i2bwins>bi2wins and i2iwins>=ii2wins:
					toremove.add(indelpos[1])
					print "...indel 2 kept"
					sys.stdout.flush()
					continue
				else:
					print "...both indels kept"
					sys.stdout.flush()
				

	if doreassessment:
	
		#check that overlapping indels aren't actually one indel
		print "Removing indels that are within another indel"
		sys.stdout.flush()
		
		for i, indelpos in enumerate(indelpositions):#[:5]:
			
			if indelpos[1] in toremove:
				continue
			
			j=i+1
			
			while j<len(indelpositions) and (indelpositions[i][0]-indelpositions[j][0]<=maxrlen or (indels[indelpos[1]]["start"])-indelpositions[j][0]<=maxrlen):
			
				indel2pos=indelpositions[j]
				j+=1
				
				indel=indels[indelpos[1]]
				indel2=indels[indel2pos[1]]
				
				if indel["chrom"]!=indel2["chrom"]:
					continue
				
				minstart=int(indel2["start"])
				if int(indel["start"])<minstart:
					toremove.add(indel2pos[1])
					continue
	
	
	
	#remove any indels that have been rejected
				
	for indel in toremove:
		removedindels.append(indels[indel])
		del indels[indel]
	
	
	print "Writing tab and txt output files"
	sys.stdout.flush()
	
	#print a tab file of indels
	
	tabout=open(options.output+"_indels.tab", "w")
	for indel in removedindels:
		if indel["info"]["SVTYPE"]=="DEL":
			print >> tabout, "FT   DELETION        "+str(int(indel["start"])+1)+".."+str(int(indel["info"]["END"])-1)
			print >> tabout, "FT                   /colour=16"
		elif indel["info"]["SVTYPE"]=="INS":
			print >> tabout, "FT   INSERTION       "+str(int(indel["start"])+1)+".."+str(int(indel["info"]["END"])-1)
			print >> tabout, "FT                   /colour=9"
		elif indel["info"]["SVTYPE"]=="RPL":
			print >> tabout, "FT   REPLACEMENT     "+str(int(indel["start"])+1)+".."+str(int(indel["info"]["END"])-1)
		elif indel["info"]["SVTYPE"]=="LI":
			print >> tabout, "FT   LARGE_INSERTION "+str(int(indel["start"])+1)+".."+str(int(indel["start"])+1)
			print >> tabout, "FT                   /colour=3"
		
		for i in indel:
			if i=="info":
				for j in indel[i]:
					print >> tabout, "FT                   /"+str(j).lower()+"="+str(indel[i][j])
			else:
				print >> tabout, "FT                   /"+str(i)+"="+str(indel[i])
		print >> tabout, "FT                   /note=Rejected"

	
	indelpositions=[]
	indelout=open(options.output+"_indels.txt", "w")
	
	indellist=indels.keys()
	indellist.sort()
	
	for k in indellist:
		#print indel
		indel=indels[k]
		indeloutlist=[]
		indeloutlist.append(indel["chrom"])
		indeloutlist.append(indel["start"])
		indeloutlist.append(indel["awins"])
		indeloutlist.append(indel["rwins"])
		indeloutlist.append(indel["draws"])
		indeloutlist.append(indel["info"]["SVLEN"])
		indeloutlist.append(indel["ref"])
		indeloutlist.append(indel["alt"])
		
		if indel["info"]["SVTYPE"]=="DEL":
			print >> tabout, "FT   DELETION        "+str(int(indel["start"])+1)+".."+str(int(indel["info"]["END"])-1)
			print >> tabout, "FT                   /colour=2"
		else:
			print >> tabout, "FT   INSERTION       "+str(int(indel["start"])+1)+".."+str(int(indel["info"]["END"])-1)
			print >> tabout, "FT                   /colour=4"
		
		for i in indel:
			if i=="info":
				for j in indel[i]:
					print >> tabout, "FT                   /"+str(j).lower()+"="+str(indel[i][j])
			else:
				print >> tabout, "FT                   /"+str(i)+"="+str(indel[i])
			
		
		
		print >> indelout, '\t'.join(map(str,indeloutlist))
		
		indelpositions.append([int(indel["info"]["END"]), k])
	
	indelout.close()
	tabout.close()
	indelpositions.sort()
	indelpositions.reverse()
	
	
	
	
	
	
	
	#Change all accepted indels in reference and remap data.
	print "Remapping data to altered reference"
	sys.stdout.flush()
	
	new_contigs={}
	for contig in sequences:
		new_contigs[contig]=sequences[contig]
	
	for indelpos in indelpositions:#[:4]:
			indel=indels[indelpos[1]]
			new_contigs[indel["chrom"]]=new_contigs[indel["chrom"]][:indel["start"]-1]+indel["alt"]+new_contigs[indel["chrom"]][int(indel["info"]["END"])-1:]
				
	
	output=open(tmpname+".ref.dna", "w")
	for contig in reference_contig_names:
		print >> output, ">"+contig
		print >> output, new_contigs[contig]
	output.close()
			
	
	
	separate_reads_to_map(readstomap, tmpname, tmpname+"_others.bam", outputunmapped=True, outputothers=True)
	os.system("samtools index "+tmpname+"_others.bam")
	
	#map the reads
	map_reads(freads=tmpname+"_1.fastq", rreads=tmpname+"_2.fastq", ref=tmpname+".ref.dna")
	
	
	
	print "Realiging reads around indels"
	sys.stdout.flush()
	
	try: samfile = pysam.Samfile( tmpname+".bam", "rb" )
	except StandardError:
		print tmpname+".bam not a bam file"
		sys.exit() 
	
	
	
	refs=samfile.references
	lengths=samfile.lengths
	
	startposns={}
	for ref in refs:
		startposns[ref]={}
		
	for read in samfile:
		if not read.is_unmapped and read.pos not in startposns[refs[read.tid]]:
			startposns[refs[read.tid]][read.pos]=-1
	samfile.reset()
	sortedstarts={}
	for ref in startposns:
		sortedstarts[ref]= startposns[ref].keys()
		sortedstarts[ref].sort()
	#indelpositions=indelpositions[:4]
	indelpositions.reverse()
	x=0
	toadd=0
	for ref in sortedstarts:
		for pos in sortedstarts[ref]:
		
			while x<len(indelpositions) and ((indels[indelpositions[x][1]]["info"]["SVTYPE"]=="DEL" and (pos+toadd)>=indels[indelpositions[x][1]]["start"]) or (indels[indelpositions[x][1]]["info"]["SVTYPE"]=="INS" and (pos+toadd)>=(indels[indelpositions[x][1]]["start"]+int(indels[indelpositions[x][1]]["info"]["SVLEN"]))) ):
				toadd+=int(indels[indelpositions[x][1]]["info"]["SVLEN"])*-1
				x+=1
				
			diff=0
			insdiff=0
			readlenwithins=maxrlen
			
			#print x, indelpositions[x], indelpositions[x][1]
			
			if x<len(indelpositions) and indels[indelpositions[x][1]]["info"]["SVTYPE"]=="INS" and (pos+toadd)>=indels[indelpositions[x][1]]["start"] and (pos+toadd)<(indels[indelpositions[x][1]]["start"]+int(indels[indelpositions[x][1]]["info"]["SVLEN"])):
				diff=indels[indelpositions[x][1]]["start"]+int(indels[indelpositions[x][1]]["info"]["SVLEN"])-(pos+toadd)
				insdiff=int(indels[indelpositions[x][1]]["info"]["SVLEN"])-diff
				startposns[ref][pos]=[pos+toadd-insdiff]
			else:
				startposns[ref][pos]=[pos+toadd]
			y=x
			readindels=[]
			
			while y<len(indelpositions) and ((startposns[ref][pos][0])+readlenwithins)>=indels[indelpositions[y][1]]["start"]:
				if indels[indelpositions[y][1]]["info"]["SVTYPE"]=="DEL":
					readindels.append([indels[indelpositions[y][1]]["start"]-(startposns[ref][pos][0]),int(indels[indelpositions[y][1]]["info"]["SVLEN"])*-1,indels[indelpositions[y][1]]["info"]["SVTYPE"] ])
					readlenwithins+=int(indels[indelpositions[y][1]]["info"]["SVLEN"])*-1
				else:
					if y==x:
						readindels.append([indels[indelpositions[y][1]]["start"]-(startposns[ref][pos][0]),int(indels[indelpositions[y][1]]["info"]["SVLEN"])-insdiff,indels[indelpositions[y][1]]["info"]["SVTYPE"] ])
					else:
						#print y, indels[indelpositions[y][1]]["start"]-(startposns[ref][pos][0]),int(indels[indelpositions[y][1]]["info"]["SVLEN"])-insdiff,indels[indelpositions[y][1]]["info"]["SVTYPE"]
						readindels.append([indels[indelpositions[y][1]]["start"]-(startposns[ref][pos][0]),int(indels[indelpositions[y][1]]["info"]["SVLEN"]),indels[indelpositions[y][1]]["info"]["SVTYPE"] ])
				y+=1
			startposns[ref][pos].append(readindels)
	
	
	reference_contig_tids={}
	for x, ref in enumerate(reference_contig_names):
		reference_contig_tids[ref]=x
	#pysam cigar states: 0 = Mapped, 1 = Insertion, 2 = Deletion, 3 = N, 4 = Soft clipped, 5 = Hard clipped, 6 = Padded
	#i[0] = indel start position after read start position. i[1] = indel length
	towrite={}
	sam=pysam.Samfile(tmpname+"_Si.bam", mode='wb', referencenames=reference_contig_names, referencelengths=reference_contig_lengths)
	for read in samfile:
		if not read.is_unmapped:
			pos=read.pos
			oricigar=read.cigar[:]
			cigcount=0
			for cig in read.cigar:
				if cig[0] in [0,1,4]:
					cigcount+=cig[1]
			
			read.tid=reference_contig_tids[refs[read.tid]]
			
			read.pos=startposns[reference_contig_names[read.tid]][pos][0]
			
			for inum in xrange (0, len(startposns[reference_contig_names[read.tid]][pos][1])):
				
				i=startposns[reference_contig_names[read.tid]][pos][1][inum][:]
				newcigar=[]
				cigar=read.cigar[:]
				refpos=0
				readpos=0
				
				for cig in cigar:
					
					
					if cig[0]==0 and cig[1]+refpos>i[0] and i[0]>=refpos:
						if (i[0]-refpos)>0:
							newcigar.append((cig[0],i[0]-refpos))
						if i[2]=="DEL":
							if i[1]>0:
								newcigar.append((2,i[1]))
							if (cig[1])-(i[0]-refpos)>0:
								
								newcigar.append((cig[0],(cig[1])-(i[0]-refpos)))
						elif i[2]=="INS":
							if ((i[0]-refpos)+i[1])+readpos>read.rlen:
								#print i[1]
								i[1]=read.rlen-((i[0]-refpos)+readpos)
								#print readpos, refpos, i[0], i[0]-refpos, ((i[0]-refpos)+readpos), read.rlen, i[1]
							if i[1]>0:
								newcigar.append((1,i[1]))
							if (cig[1])-((i[0]-refpos)+i[1])>0:
								newcigar.append((cig[0],(cig[1])-((i[0]-refpos)+i[1])))
							else:
								i[0]=float("Inf")
					
					elif cig[0]==2 and cig[1]+refpos>i[0] and i[0]>=refpos:
						if i[2]=="DEL":
							newcigar.append((cig[0],i[1]+cig[1]))
						elif i[2]=="INS":
							if i[1]>cig[1] and i[1]-cig[1]>0:
								newcigar.append((1,i[1]-cig[1]))
								i[0]=float("Inf")
							elif cig[1]-i[1]>0:
								newcigar.append((cig[0],cig[1]-i[1]))
							else:
								i[0]=float("Inf")
								
					
					elif cig[0]==1 and i[0]==refpos:
						if i[2]=="INS":
							newcigar.append((cig[0],i[1]+cig[1]))
							i[0]=float("Inf")
						elif i[2]=="DEL":
							print i, cig
							if i[1]>cig[1] and i[1]-cig[1]>0:
								newcigar.append((2,i[1]-cig[1]))
							elif cig[1]-i[1]>0:
								newcigar.append((cig[0],cig[1]-i[1]))
							else:
								i[0]=float("Inf")
							
					elif cig[0] in [0,1,4]  and cig[1]+readpos<read.rlen:
						newcigar.append(cig)
					elif cig[0] in [0,1,4]  and readpos<read.rlen:
						newcigar.append((cig[0],read.rlen-readpos))
					elif readpos<read.rlen:
						newcigar.append(cig)
					
					refpos=0
					readpos=0
					for ncig in newcigar:
						if ncig[0]==0:
							refpos+=ncig[1]
							readpos+=ncig[1]
						elif ncig[0]==2:
							refpos+=ncig[1]
						elif ncig[0] in [1,4]:
							readpos+=ncig[1]
										
					
				
				try:
					read.cigar=newcigar[:]
				except StandardError:
					print i
					print cigar, read.cigar
				
				newcigcount=0
				for cig in read.cigar:
					if cig[0] in [0,1,4]:
						newcigcount+=cig[1]
				if cigcount!=newcigcount:
				#if len(startposns[reference_contig_names[read.tid]][pos][1])>1:
					print read.qname, read.rlen, cigcount, newcigcount, oricigar, read.cigar,  startposns[reference_contig_names[read.tid]][pos]
		
		if read.qname in towrite:
			if not read.is_unmapped and not towrite[read.qname].is_unmapped:
				read.pnext=towrite[read.qname].pos
				towrite[read.qname].pnext=read.pos
				isize=(towrite[read.qname].pnext-towrite[read.qname].pos)+read.qlen
				read.tlen=isize*-1
				towrite[read.qname].tlen=isize
			sam.write(towrite[read.qname])
			sam.write(read)
		else:
			towrite[read.qname]=read
		
	
	sam.close()
	
	print "Creating final bam"
	sys.stdout.flush()
	indels=""
	towrite=""
	os.system("rm -f "+tmpname+"*.sam "+tmpname+"*.fastq")
	os.system("samtools sort "+tmpname+"_Si.bam "+tmpname+"_sort")
	os.system("rm -f "+tmpname+"_Si.bam")
	os.system("samtools index "+tmpname+"_sort.bam")
	os.system("samtools merge -f "+options.output+".bam "+tmpname+"_sort.bam "+tmpname+"_others.bam")
	os.system("samtools index "+options.output+".bam")
	
	print "Cleaning up"
	sys.stdout.flush()
	if options.tmpname=="":
		os.system("rm -f "+tmpname+"*")
	print "Done"
	sys.stdout.flush()

