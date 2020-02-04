#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math


def revcomp(sequence):
	rev=sequence[::-1]
	revcomp=''
	d={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'g':'c', 'c':'g'}
	for i in rev:
		if d.has_key(i):
			revcomp=revcomp+d[i]
		else:
			revcomp=revcomp+i
	
	return revcomp



if (len (sys.argv)!=3) or '-h' in sys.argv[1:]:
	print "clusters_from_blast.py <blast file> <contigs file>"
	sys.exit()
	
lines=open(sys.argv[2], "rU").read().split('>')[1:]

contigseqs={}

for line in lines:
	words=line.strip().split('\n')
	contigseqs[words[0].split()[0]]=''.join(words[1:])


blastlines=open(sys.argv[1],'rU').readlines()

hits={}

conversiondict={}

for line in blastlines:
	words=line.split()

	if not hits.has_key(words[0].split('_')[1]):
		hits[words[0].split('_')[1]]=[]
	if words[0] not in hits[words[0].split('_')[1]]:
		hits[words[0].split('_')[1]].append(words[0])


for hit in hits.keys():
	print hit
	output=open(hit+'.fasta','w')
	
	for contig in hits[hit]:
		print >> output, '>'+contig
		print >> output, contigseqs[contig]
	
	output.close()
