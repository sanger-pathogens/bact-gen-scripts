#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math


if (len (sys.argv)!=2) or '-h' in sys.argv[1:]:
	print "GC_moving_average.py <fasta file>"
	sys.exit()

print sys.argv[1]

seq=''.join(open(sys.argv[1], "rU").read().split('\n')[1:])

output=open("GCmovingwindow.plot", "w")

for x in range(0,len(seq)):

	if x-27<0:
		tempseq=seq[x-27:]+seq[0:x+27]
	elif x+27>len(seq):
		tempseq=seq[x-27:]+seq[0:(x+27)-len(seq)]
	else:
		tempseq=seq[x-27:x+27]


	tempseq=tempseq.upper()

	
	GCcount=0
	GCpercent=0.0
	for y in tempseq:
		if y in ['G','C']:
			GCcount=GCcount+1
	GCpercent=int(float(GCcount)/55*100)
	print >> output, GCpercent
	