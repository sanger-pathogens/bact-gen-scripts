#!/usr/bin/env python
import string, re, gzip
import os, sys


if len(sys.argv)!=4 or '-h' in sys.argv[1:]:
	print "concatenatefastqs.py <forward fastq file> <reverse fastq file> <output file>"
	sys.exit()
	
if sys.argv[1].split('.')[-1]=='gz':
	forward=gzip.open(sys.argv[1], 'r')
else:
	forward=open(sys.argv[1], 'r')
	
if sys.argv[2].split('.')[-1]=='gz':
	reverse=gzip.open(sys.argv[2], 'r')
else:
	reverse=open(sys.argv[2], 'r')

output=open(sys.argv[3],'w')


for line in forward:
	flines=[line.strip()]
	for y in range(0,3):
		flines.append(forward.next().strip())
	
	rlines=[]
	for y in range(0,4):
		rlines.append(reverse.next().strip())

	print >> output, flines[0]
	print >> output, flines[1]+rlines[1]
	print >> output, flines[2]
	print >> output, flines[3]+rlines[1]

output.close()