#!/usr/bin/env python
#/usr/bin/env python
import string
import os, sys
import pysam
from numpy import mean, std, median

#
##############################################
## Function to reverse complement a sequence #
##############################################
#
#def revcomp(sequence):
#	rev=sequence[::-1]
#	revcomp=''
#	d={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'g':'c', 'c':'g', "n":"n", "N":"N"}
#	for i in rev:
#		if d.has_key(i):
#			revcomp=revcomp+d[i]
#		else:
#			revcomp=revcomp+i
#	
#	return revcomp
refseqs={}
if sys.argv[1]=="-e":
	doerrors=True
	reflines=open(sys.argv[2]).read().split(">")[1:]
	for ref in reflines:
		seqlines=ref.split("\n")
		refname=seqlines[0].split()[0]
		refseqs[refname]=''.join(seqlines[1:]).upper()
	#if len(refseq.split(">"))>1:
	#	print "Ref must be one contig"
	#	sys.exit()

	filestart=3
else:
	doerrors=False
	filestart=1
if doerrors:
	print "File\tContig\tContig length\tTotal Reads\tMean length\tMapped\tProper pairs\tUnmapped\tChimeras\tErrors\tInsertions\tDeletions\tErrors per mapped base\tInsertions per mapped base\tDeletions per mapped base\tMean depth\tStd depth\tMedian depth\tCount of unmapped bases in contig"
else:
	print "File\tContig\tContig length\tTotal\tMean length\tMapped\tProper pairs\tUnmapped\tChimeras"

for filename in sys.argv[filestart:]:
	if filename.split(".")[-1]=="bam":
		samfile = pysam.Samfile( filename, "rb" )
	elif filename.split(".")[-1]=="sam":
		samfile = pysam.Samfile( filename, "r" )
	else:
		print "Not a bam file"
		sys.exit()
		
	refs=samfile.references
	lengths=samfile.lengths
	
	for x, ref in enumerate(refs):
		if not ref in refseqs:
			print "Error! Reference sequences in fasta and bam do not match"
			sys.exit()
		elif len(refseqs[ref])!=lengths[x]:
			print "Error! Reference sequences in fasta and bam are not the same length"
			sys.exit()
			
	
	refstats={}
	
	unmapped=0.0
	chimera=0.0
	proper_pair=0.0
	mapped=0.0
	total=0.0
	errors=0.0
	insertions=0.0
	deletions=0.0
	mappedlen=0.0
	totallength=0.0
	depth={}
	
	totreflen=0
	for x, ref in enumerate(refs):
		refstats[ref]={"proper_pair":0.0, "mapped":0.0, "total":0.0, "errors":0.0, "mappedlen":0.0, "insertions":0.0, "deletions":0.0}
		depth[ref]=[0]*lengths[x]
		totreflen+=lengths[x]
	
	
	
	for read in samfile:
		total+=1
		totallength+=len(read.seq)
		if read.is_proper_pair and read.rname==read.mrnm:
			refstats[samfile.getrname(read.rname)]["proper_pair"]+=0.5
			refstats[samfile.getrname(read.rname)]["mapped"]+=1
			mapped+=1
			proper_pair+=0.5
		elif read.is_unmapped:
			unmapped+=1
		elif not read.is_unmapped and read.rname!=read.mrnm:
			refstats[samfile.getrname(read.rname)]["mapped"]+=1
			mapped+=1
			if not read.mate_is_unmapped:
				chimera+=0.5
		
		else:
			refstats[samfile.getrname(read.rname)]["mapped"]+=1
			mapped+=1
		
		if read.rname!=-1:
			refstats[samfile.getrname(read.rname)]["total"]+=1

		if doerrors and not read.is_unmapped:# and not read.is_reverse:
			start=read.pos
			readpos=0
			refpos=start
#			if read.is_reverse:
#				readseq=revcomp(read.seq.upper())
#			else:
#				readseq=read.seq.upper()
			readseq=read.seq.upper()
#			print read.is_reverse, readseq
#			print refseq[refpos:refpos+100]
			
			for cig in read.cigar:
				if cig[0]==0:
					for x in range(0,cig[1]):
						if readseq[readpos]!=refseqs[samfile.getrname(read.rname)][refpos] and refseqs[samfile.getrname(read.rname)][refpos] in ["A", "C", "G", "T"]:
							errors+=1
							refstats[samfile.getrname(read.rname)]["errors"]+=1
						readpos+=1
						
						depth[samfile.getrname(read.rname)][refpos]+=1
						refpos+=1
					mappedlen+=cig[1]
					refstats[samfile.getrname(read.rname)]["mappedlen"]+=cig[1]
				elif cig[0]==1:
					insertions+=1
					refstats[samfile.getrname(read.rname)]["insertions"]+=1
					readpos+=cig[1]
				elif cig[0]==2:
					deletions+=1
					refstats[samfile.getrname(read.rname)]["deletions"]+=1
					refpos+=cig[1]
				elif cig[0]==4:
					readpos+=cig[1]
				else:
					print cig
		elif not read.is_unmapped:
			for cig in read.cigar:
				if cig[0]==0:
					mappedlen+=cig[1]
					refstats[samfile.getrname(read.rname)]["mappedlen"]+=cig[1]

	samfile.close()	
			#print errors, errors/mappedlen, insertions, deletions
	if doerrors:
		totdepth=[]
		depthstats={}
		
		for ref in depth:
			totdepth+=depth[ref]
			depthstats[ref]=[mean(depth[ref]), std(depth[ref]), median(depth[ref]), depth[ref].count(0)]
			
		totdepthstats=[mean(totdepth), std(totdepth), median(totdepth), totdepth.count(0)]
	
	
	
	if not doerrors:
		if len(refs)>1:
			if total==0:
				totaltoreport="-"
			else:
				totaltoreport=str(totallength/total)
			print "\t".join([filename, "All contigs", str(totreflen), str(total), totaltoreport, str(mapped), str(proper_pair), str(unmapped), str(chimera)])
		for x, ref in enumerate(refs):
			print "\t".join([filename, ref, str(lengths[x]), "-", "-", str(refstats[ref]["mapped"]), str(refstats[ref]["proper_pair"]), "-", "-"])
	else:
		if len(refs)>1:
			
			if total==0:
				totaltoreport="-"
			else:
				totaltoreport=str(totallength/total)
			if mappedlen==0:
				errorproportion="-"
				insertproportion="-"
				deletionproportion="-"
			else:
				errorproportion=str(errors/mappedlen)
				insertproportion=str(insertions/mappedlen)
				deletionproportion=str(deletions/mappedlen)
				
			print "\t".join([filename, "All contigs", str(totreflen), str(total), totaltoreport, str(mapped), str(proper_pair), str(unmapped), str(chimera), str(errors), str(insertions), str(deletions), errorproportion, insertproportion, deletionproportion, str(totdepthstats[0]), str(totdepthstats[1]), str(totdepthstats[2]), str(totdepthstats[3])])
		for x, ref in enumerate(refs):
			
			if refstats[ref]["mappedlen"]==0:
				errorproportion="-"
				insertproportion="-"
				deletionproportion="-"
			else:
				errorproportion=str(refstats[ref]["errors"]/refstats[ref]["mappedlen"])
				insertproportion=str(refstats[ref]["insertions"]/refstats[ref]["mappedlen"])
				deletionproportion=str(refstats[ref]["deletions"]/refstats[ref]["mappedlen"])
		
		
			print "\t".join([filename, ref, str(lengths[x]), "-", "-", str(refstats[ref]["mapped"]), str(refstats[ref]["proper_pair"]), "-", "-", str(refstats[ref]["errors"]), str(refstats[ref]["insertions"]), str(refstats[ref]["deletions"]), errorproportion, insertproportion, deletionproportion, str(depthstats[ref][0]), str(depthstats[ref][1]), str(depthstats[ref][2]), str(depthstats[ref][3])])
#	if doerrors:
#		print "Errors\tInsertions\tDeletions\tMapped bases"
#		print str(errors), str(insertions), str(deletions)
#		print "Errors per mapped base\tInsertions per mapped base\tDeletions per mapped base"
#		print str(errors/mappedlen), str(insertions/mappedlen), str(deletions/mappedlen)
	
	    	
			
