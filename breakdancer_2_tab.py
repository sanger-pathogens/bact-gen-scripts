#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math
import time

if len(sys.argv)!=3:
	print "breakdancer_2_tab.py <breakdancer output file> <output tab file name>"

try:
	infile=open(sys.argv[1], "rU")
except StandardError:
	print "Cannot open file", sys.argv[1]
	sys.exit()

output=open(sys.argv[2], "w")

print >> output, "ID   SV"

for line in infile:
	if len(line)==0 or line[0]=="#":
		continue
	
	words=line.strip().split()
	if len(words)<7:
		print words
	if words[6] in ["DEL", "INS"]:
		location=words[1]+'..'+words[4]
	else:
		location=words[1]
	
	print >> output, "FT   "+words[6]+" "*(16-len(words[6]))+location
	
	if words[6] in ["DEL", "INS"]:
		print >> output, "FT                   /size="+words[7]
	
	print >> output, "FT                   /confidence="+words[8]
	print >> output, "FT                   /supporting_reads="+words[9]
	
	if words[6]=="DEL":
		if int(words[8])==99:
			colour="6"
		else:
			colour="12"
	elif words[6]=="INS":
		if int(words[8])==99:
			colour="4"
		else:
			colour="5"
	elif words[6]=="INV":
		if int(words[8])==99:
			colour="3"
		else:
			colour="8"
	elif words[6]=="ITX":
		if int(words[8])==99:
			colour="10"
		else:
			colour="7"
	elif words[6]=="ITX":
		if int(words[8])==99:
			colour="1"
		else:
			colour="13"
	else:
		colour="11"
		
	print >> output, "FT                   /colour="+colour

output.close()
	

	