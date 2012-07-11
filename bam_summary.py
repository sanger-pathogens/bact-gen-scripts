#!/usr/bin/env python
#/usr/bin/env python
import string
import os, sys
import pysam


if sys.argv[1]=="-e":
	doerrors=True
	reflines=open(sys.argv[2]).read().split("\n")
	refname=reflines[0].split()[0][1:]
	refseq=''.join(reflines[1:])
	#if len(refseq.split(">"))>1:
	#	print "Ref must be one contig"
	#	sys.exit()

	filestart=3
else:
	doerrors=False
	filestart=1
if doerrors:
	print "File\tContig\tTotal\tMean length\tMapped\tProper_pairs\tUnmapped\tChimeras\tErrors\tInsertions\tDeletions"
else:
	print "File\tContig\tTotal\tMean length\tMapped\tProper_pairs\tUnmapped\tChimeras"

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
	
	for ref in refs:
		refstats[ref]={"proper_pair":0.0, "mapped":0.0, "total":0.0}
	
	
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
		elif read.rname!=read.mrnm:
			refstats[samfile.getrname(read.rname)]["mapped"]+=1
			mapped+=1
			chimera+=0.5
		
		else:
			refstats[samfile.getrname(read.rname)]["mapped"]+=1
			mapped+=1
		
		if read.rname!=-1:
			refstats[samfile.getrname(read.rname)]["total"]+=1

		if doerrors and not read.is_unmapped and not read.is_reverse:
			start=read.pos
			readpos=0
			refpos=start
			
			for cig in read.cigar:
				if cig[0]==0:
					for x in range(0,cig[1]):
						if read.seq[readpos].upper()!=refseq[refpos].upper():
							errors+=1
						readpos+=1
						refpos+=1
					mappedlen+=cig[1]
				elif cig[0]==1:
					insertions+=1
					readpos+=cig[1]
				elif cig[0]==2:
					deletions+=1
					refpos+=cig[1]
				elif cig[0]==4:
					readpos+=cig[1]
				else:
					print cig
		elif not read.is_unmapped:
			for cig in read.cigar:
				if cig[0]==0:
					mappedlen+=cig[1]

			
			#print errors, errors/mappedlen, insertions, deletions
	
	if doerrors:
		if len(refs)>1:
			print "\t".join([filename, "All", str(total), str(totallength/total), str(mapped), str(proper_pair), str(unmapped), str(chimera)])
		for ref in refs:
			print "\t".join([filename, ref, str(refstats[ref]["total"]), str(refstats[ref]["mapped"]), str(refstats[ref]["proper_pair"]), "-", "-"])
	else:
		if len(refs)>1:
			print "\t".join([filename, "All", str(total), str(totallength/total), str(mapped), str(proper_pair), str(unmapped), str(chimera), str(errors/mappedlen), str(insertions/mappedlen), str(deletions/mappedlen)])
		for ref in refs:
			print "\t".join([filename, ref, str(refstats[ref]["total"]), str(refstats[ref]["mapped"]), str(refstats[ref]["proper_pair"]), "-", "-"])
	if doerrors:
		print "Errors\tInsertions\tDeletions\tMapped bases"
		print str(errors), str(insertions), str(deletions)
		print "Errors per mapped base\tInsertions per mapped base\tDeletions per mapped base"
		print str(errors/mappedlen), str(insertions/mappedlen), str(deletions/mappedlen)
	samfile.close()
	    	
			
