#!/usr/bin/env python
import string, re, gzip
import os, sys


if len(sys.argv)!=4 or '-h' in sys.argv[1:]:
	print "cigar2plot.py <reference fasta file> <cigar file (can be zipped)> <output file>"
	sys.exit()

if sys.argv[1].split('.')[-1]=='gz':
	reflen=len(''.join(gzip.open(sys.argv[1], 'r').read().split('>')[1].split('\n')[1:]))	
else:
	reflen=len(''.join(open(sys.argv[1], 'rU').read().split('>')[1].split('\n')[1:]))

depth=[0]*reflen


if sys.argv[2].split('.')[-1]=='gz':
	lines=gzip.open(sys.argv[2], 'r').readlines()
else:
	lines=open(sys.argv[2], 'r').readlines()

for line in lines:
	words=line.split()
	for x in range(int(words[6])-1,int(words[7])):
		depth[x]=depth[x]+1

output=open(sys.argv[3], "w")

for x in range(0,len(depth)):
	print >> output, depth[x]

output.close()