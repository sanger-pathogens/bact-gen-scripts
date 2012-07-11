#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math
import time

if len(sys.argv)!=3:
	print "pindel_2_tab.py <pindel output file> <output tab file name>"

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
	
	if len(words)<2 or words[1] not in ["I", "D", "LI", "INV", "TD"]:
		continue
	
	
	if words[1] in ["D", "I", "INV", "TD"]:
		location=str(int(words[9])+1)+'..'+str(int(words[10])-1)
	elif words[1]=="LI":
		location=str(int(words[4])+1)+'..'+str(int(words[6])-1)
	
	print >> output, "FT   "+words[1]+" "*(16-len(words[1]))+location
	
	if words[1] in ["D", "I", "INV", "TD"]:
		if words[4]!="0" and words[1]!="INV":
			print >> output, "FT                   /ins_length="+words[4]
			print >> output, "FT                   /ins_sequence="+words[5]
		print >> output, "FT                   /supports="+words[15]
		print >> output, "FT                   /unique_supports="+words[16]
		print >> output, "FT                   /upstream_supports="+words[18]
		print >> output, "FT                   /upstream_unique_supports="+words[19]
		print >> output, "FT                   /downstream_supports="+words[21]
		print >> output, "FT                   /downstream_unique_supports="+words[22]
		print >> output, "FT                   /S1="+words[24]
		print >> output, "FT                   /SUM_mapping_score="+words[26]
	
	if words[1]=="D":
		print >> output, "FT                   /colour=6"
	elif words[1]=="I":
		print >> output, "FT                   /colour=4"
	elif words[1]=="LI":
		print >> output, "FT                   /colour=5"
	elif words[1]=="INV":
		print >> output, "FT                   /colour=3"
	elif words[1]=="TD":
		print >> output, "FT                   /colour=9"
	
	


output.close()
	

	