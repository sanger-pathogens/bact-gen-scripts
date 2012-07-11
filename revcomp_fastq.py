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
	
fastqfile=sys.argv[1]

lines=open(fastqfile,'rU').readlines()

output=open(fastqfile.replace('.fastq', '_revcomp.fastq'), 'w')

for x in range(0,len(lines),4):

	print >> output, lines[x].strip()
	
	revseq=revcomp(lines[x+1].strip())
	
	print >> output, revseq
	
	print >> output, lines[x+2].strip()
	
	print >> output, lines[x+3].strip()[::-1]

output.close()
		
