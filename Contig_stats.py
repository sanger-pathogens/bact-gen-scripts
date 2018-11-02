#!/usr/bin/env python
import os, sys, string, numpy
from Bio.Seq import Seq
sys.path.extend(map(os.path.abspath, ['/nfs/pathogen/sh16_scripts/modules/']))
from Si_SeqIO import *




for filename in sys.argv[1:]:
	try:
		contigs=SeqIO.parse(open(filename), "fasta")
	except StandardError:
		print "Could not open file", filename
		continue
	
	print filename
	heading=['contig','length','length without Ns','GC%']
	
	lengths=[]
	GCs=[]
	
	for x in ["A","C","G","T"]:
		for y in ["A","C","G","T"]:
			heading.append(x+y)
	print '\t'.join(heading)
	for contig in contigs:
	
		dinucleotides={}
		for x in ["A","C","G","T"]:
			dinucleotides[x]={}
			for y in ["A","C","G","T"]:
				dinucleotides[x][y]=0
		
		seq=str(contig.seq).upper()
		length=len(seq)
		
		lengths.append(length)
		
		lengthnons=len(seq.upper().replace("N",""))
		GC=(float(len(seq.upper().replace("N","").replace("A","").replace("T","")))/lengthnons)*100
		
		GCs.append(GC)
		
		for x in range(0, length-1):
			if seq[x]!="N" and seq[x+1]!="N":
				dinucleotides[seq[x]][seq[x+1]]+=1
		outline=[contig.id,str(length),str(lengthnons),str(GC)]
		
		for x in ["A","C","G","T"]:
			for y in ["A","C","G","T"]:
				outline.append(str(dinucleotides[x][y]))
		print '\t'.join(outline)
		
	print "Total length =", numpy.sum(lengths)
	print "Number of contigs =", len(lengths)
	print "Mean length =", numpy.mean(lengths)
	print "Standard deviation of lengths =", numpy.std(lengths)
	print "Maximum length =", numpy.max(lengths)
	print "Minimum length =", numpy.min(lengths)
	lengths.sort()
	lengths.reverse()
	fifty=float(numpy.sum(lengths))/2
	count=0
	sum=0
	
	while sum<fifty:
		N50=lengths[count]
		sum+=lengths[count]
		count+=1
	
	print "N50 =", N50
	print "N50n =", count
	print "Mean GC% =", numpy.mean(GCs)
	print "GC% standard deviation =", numpy.std(GCs)
	print "Median GC% =", numpy.median(GCs)
	print "Maximum GC% =", numpy.max(GCs)
	print "Minimum GC% =", numpy.min(GCs)
	
	
	
