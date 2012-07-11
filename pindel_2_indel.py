#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math
import time

if len(sys.argv)!=3:
	print "pindel_2_tab.py <breakdancer output file> <output indel file name>"

try:
	infile=open(sys.argv[1], "rU")
except StandardError:
	print "Cannot open file", sys.argv[1]
	sys.exit()

output=open(sys.argv[2], "w")


for line in infile:
	if len(line)==0 or line[0]=="#":
		continue
	
	
	words=line.strip().split()
	
	if len(words)<2 or words[1] not in ["I", "D"]:
		continue
	
	
	if words[1]=="I":
		indeltype="+"
	else:
		indeltype="-"
	
	print >> output, words[7], words[9], indeltype, words[5], (int(words[10])-(int(words[9])))-1
	
		


output.close()
	

	