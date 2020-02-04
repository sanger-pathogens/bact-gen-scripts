#!/usr/bin/env python
import string, re
import os, sys

if len (sys.argv)!=3 or '-h' in sys.argv[1:]:
	print "Alignment_to_coveragecount.py <alignment file> <outfile name>"
	sys.exit()


alignment=sys.argv[1]

output=open(sys.argv[2],"w")

lines=open(alignment, "rU").read().split('>')[1:]

sequences={}
otherseq=''

for line in lines:
	words=line.strip().split('\n')
	name=words[0]
	sequences[name]=''.join(words[1:])


y=0
for x in range(len(sequences[name])):
	count=0
	for sequence in sequences:
		if sequences[sequence][x] in ["A", "C", "G", "T", "a", "c", "g", "t"]:
			count+=1
	print >> output, count

