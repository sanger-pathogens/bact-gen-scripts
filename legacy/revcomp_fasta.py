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


data=open(sys.argv[1], 'rU').read().split(">")[1:]

sequences={}
for datum in data:
	name=datum.split("\n")[0].split()[0]
	sequences[name]=''.join(datum.split("\n")[1:])

output=open(sys.argv[2], 'w')

for sequence in sequences:

	print >> output, ">"+sequence
	
	revseq=revcomp(sequences[sequence])
	
	print >> output, revseq
	

output.close()
		
