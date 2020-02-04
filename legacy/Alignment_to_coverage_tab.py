#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math

if len (sys.argv)!=3 or '-h' in sys.argv[1:]:
	print "Alignment_to_coverage_tab.py <reference name> <alignment file> > <outfile name>"
	sys.exit()

ref=sys.argv[1]
alignment=sys.argv[2]

lines=open(alignment, "rU").read().split('>')[1:]


print "ID   Mapping"

sequences={}
otherseq=''

for line in lines:
	words=line.strip().split('\n')
	name=words[0]
	sequences[name]=''.join(words[1:])
	if name!=ref:
		otherseq=name

inblock='n'
y=0
for x in range(len(sequences[ref])):
	if sequences[ref][x]=='-':
		continue
	elif sequences[otherseq][x]!='-':
		if inblock=='n':
			inblock='y'
			outstring= "FT   contig          "+str(y+1)+'..'
	else:
		if inblock=='y':
			inblock='n'
			outstring=outstring+str(y+1)
			print outstring
			print 'FT                   /colour="4"'
	y=y+1