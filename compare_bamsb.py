#!/usr/bin/env python
#/usr/bin/env python
import string, re
import os, sys, getopt, random, math
import time
import pysam


#CIGAR operators: 0=match, 1=insertion, 2=deletion, 3=skipped region from reference, 4=soft clip on read (i.e. sequence still in SEQ, 5=hard clip on read (sequence not in SEQ), 6=Padding (silent deletion from padded reference - for multiple alignment)

def calculate_cigar_length(cigar_sequence):
	length=0
	
	for f in cigar_sequence:
		if f[0] in [0,1,3,4,5]:
			length+=f[1]
	
	return length


#open the two SAM/BAM format input files (must be sorted by name)

if sys.argv[1].split(".")[-1]=="bam":
	samfile1 = pysam.Samfile( sys.argv[1], "rb" )
	#pysam.sort( "-n", sys.argv[1], "sort1" )
elif sys.argv[1].split(".")[-1]=="sam":
	samfile1 = pysam.Samfile( sys.argv[1], "r" )

if sys.argv[2].split(".")[-1]=="bam":
	samfile2 = pysam.Samfile( sys.argv[2], "rb" )
	#pysam.sort( "-n", sys.argv[2], "sort2" )
elif sys.argv[2].split(".")[-1]=="sam":
	samfile2 = pysam.Samfile( sys.argv[2], "r" )

samout = pysam.Samfile("test.bam", "wb", template=samfile1)



contigmatch={}

for x in range(len(samfile1.references)):
	contigmatch[x]={}

#reads_iter1=samfile1
reads_iter1=samfile1

reads_iter2=samfile2


samforpileup1 = pysam.Samfile( "tmp.bam", "rb" )
samforpileup2 = pysam.Samfile( "tmp2.bam", "rb" )




potential_joins={}

for count, read1 in enumerate(reads_iter1):
	
	read2 = reads_iter2.next()
	read1mate = reads_iter1.next()
	read2mate = reads_iter2.next()
	
	
	if read1.is_read2:
		readtmp=read1mate
		read1mate=read1
		read1=readtmp
	if read2.is_read2:
		readtmp=read2mate
		read2mate=read2
		read2=readtmp
#	
	
#
#	
	###########################################################
	# Deal with reads that map perfectly and the same in both #
	###########################################################


	if  not read1.is_unmapped and read1.rname==0 and not read2.is_unmapped and read1.is_proper_pair and read2.is_proper_pair:
		
		if not contigmatch[read1.rname].has_key(read1.pos):
			contigmatch[read1.rname][read1.pos]={}
		if not contigmatch[read1.rname][read1.pos].has_key(read2.rname):
			contigmatch[read1.rname][read1.pos][read2.rname]=[]
		contigmatch[read1.rname][read1.pos][read2.rname].append(read2.pos)
		
		
		pileup_iter1=samforpileup1.fetch(samfile1.references[0], read1.pos, read1.pos)
		
		
		for f in pileup_iter1:
			print str(f)
		pileup_iter2=samforpileup2.fetch(samfile2.references[read2.rname], read2.pos, read2.pos)
		for f in pileup_iter2:
			print str(f)
			
		sys.exit()
	
	if  not read1mate.is_unmapped and read1mate.rname==0 and not read2mate.mate_is_unmapped and read1mate.is_proper_pair and read2mate.is_proper_pair:
	
		if not contigmatch[read1.rname].has_key(read1.mpos):
			contigmatch[read1.rname][read1.mpos]={}
		if not contigmatch[read1.rname][read1.mpos].has_key(read2.rname):
			contigmatch[read1.rname][read1.mpos][read2.rname]=[]
		contigmatch[read1.rname][read1.mpos][read2.rname].append(read2.mpos)
	

		pileup_iter1=samforpileup1.fetch(samfile1.references[0], read1.mpos, read1.mpos)
		for f in pileup_iter1:
			print str(f)
		pileup_iter2=samforpileup2.fetch(samfile2.references[read2.rname], read2.mpos, read2.mpos)
		for f in pileup_iter2:
			print str(f)

		sys.exit()
#	for x, f in enumerate(read1.cigar):
#		if f[1]!=2:
#			if not contigmatch[read1.rname].has_key(x+read1.pos):
#				contigmatch[read1.rname][x+read1.pos]={}
#			if not contigmatch[read1.rname][x+read1.pos].has_key(read2.rname):
#				contigmatch[read1.rname][x+read1.pos][read2.rname]=0
#			contigmatch[read1.rname][x+read1.pos][read2.rname] += 1

	
#for n in range(0, samfile1.lengths[0]):
#	if contigmatch[0].has_key(n):
#		print n, contigmatch[0][n]
		
		
		
		
		
		