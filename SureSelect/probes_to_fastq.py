#!/usr/bin/env python
import os, sys

if len(sys.argv)<3:
	print "probes_to_fastq.py <probes> <output>"
	sys.exit()

output=open(sys.argv[2],"w")

for line in open(sys.argv[1]):
	if len(line)==0 or line[0]=="#":
		continue
	words=line.strip().replace(" ", "_").split("\t")
	if len(words)<8:
		continue
	print >> output, ">"+'_'.join([words[1], words[0], words[6], words[7]])
	print >> output, words[2]
output.close()
