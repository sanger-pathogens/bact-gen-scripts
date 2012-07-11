#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math


if (len (sys.argv)!=2) or '-h' in sys.argv[1:]:
	print "fastq_GC.py <fastq file>"
	sys.exit()

print sys.argv[1]

lines=open(sys.argv[1], "rU").read().split('\n')

basehash={"A":0,"C":0,"G":0,"T":0,"N":0}
GChisto=[0]*55

for x in range(1,len(lines),4):
#	if len(lines[x].split("#"))<2:
#		continue
	GCcount=0
	for base in lines[x]:
		basehash[base]=basehash[base]+1
		if base in ["G", "C"]:
			GCcount=GCcount+1
	GChisto[GCcount]=GChisto[GCcount]+1

output=open("GCplots.csv", "a")
print >> output, sys.argv[1]+","+",".join(map(str, GChisto))
output.close()


total=0.0
GC=0.0
for base in basehash.keys():
	print base, basehash[base]
	if base in ["G", "C"]:
		GC=GC+basehash[base]
		total=total+basehash[base]
	elif base in ["A", "T"]:
		total=total+basehash[base]
		

print "GC content = "+str(((GC/total)*100))+"%"