#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math


if ((len (sys.argv)!=2) and (len (sys.argv)!=3)) or '-h' in sys.argv[1:]:
	print "GC_moving_average.py <fasta file> <window size>"
	sys.exit()

if len (sys.argv)==3:
	try:
		window=float(int(sys.argv[2]))
	except StandardError:
		print "Window size must be a positive integer"
		sys.exit()
else:
	window=100

window+=1
halfwindow=int(window)/2

#print sys.argv[1]

seq=''.join(open(sys.argv[1], "rU").read().split('\n')[1:])

#output=open("GCmovingwindow.plot", "w")

for x in range(0,len(seq)):

	if x-49<0:
		tempseq=seq[x-halfwindow:]+seq[0:x+halfwindow]
	elif x+49>len(seq):
		tempseq=seq[x-halfwindow:]+seq[0:(x+halfwindow)-len(seq)]
	else:
		tempseq=seq[x-halfwindow:x+halfwindow]


	tempseq=tempseq.upper()

	
	GCcount=0
	GCpercent=0.0
	for y in tempseq:
		if y in ['G','C']:
			GCcount=GCcount+1
	GCpercent=int(float(GCcount)/window*100)
	print GCpercent
	
