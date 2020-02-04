#!/usr/bin/env python

import os, sys

suffix=""

if len(sys.argv[1:])==1:
	inputfile=sys.argv[1]
	outputfile=sys.argv[1]+".oneseq.fasta"
elif len(sys.argv[1:])==2:
	inputfile=sys.argv[1]
	outputfile=sys.argv[2]
elif len(sys.argv[1:])==3:
	inputfile=sys.argv[1]
	outputfile=sys.argv[2]
	suffix=sys.argv[3]
else:
	print "mfa_2_oneseq.py multifasta [output name] [suffix to add to sequence name]"
	sys.exit()

try:
	output=open(outputfile, "w")
except StandardError:
	print "Error: Could not open", outputfile
	sys.exit()

print >> output, ">"+os.path.basename(inputfile)+suffix

for line in open(inputfile, 'rU'):
	line=line.rstrip()
	if len(line)>0 and line[0]!=">":
		print >> output, line
	