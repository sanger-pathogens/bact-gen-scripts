#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math


if  len (sys.argv)!=4 or '-h' in sys.argv[1:]:
	print "get_USA300_snpsites.py <USA300.aln fastafile> <SNP file> <output file>"
	sys.exit()

alnfile=sys.argv[1]

lines=open(alnfile, "rU").read().split('>')[1:]

sequences={}


for line in lines:
	words=line.strip().split('\n')
	
	if not sequences.has_key(words[0]):
		sequences[words[0]]=''.join(words[1:])
	else:
		sequences[words[0]]=sequences[words[0]]+''.join(words[1:])

USA300=[]
TW20=[]
for x in range(0, len(sequences['My_core'])):
	if sequences['My_core'][x]!='-':
		USA300.append(sequences['USA300_TCH1516'][x])
		TW20.append(sequences['TW20'][x])


snpfile=sys.argv[2]

lines=open(snpfile, 'rU').readlines()

output=open(sys.argv[3],'w')

outline=''
for line in lines[1:]:
	words=line.split()
	posn=int(words[0])-1
	outline=outline+USA300[posn]

print >> output, ">USA300"
print >> output, outline
