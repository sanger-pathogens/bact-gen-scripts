#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math


if (len (sys.argv)!=2) or '-h' in sys.argv[1:]:
	print "strain_snps.py <snp file>"
	sys.exit()

lines=open(sys.argv[1], "rU").readlines()

strains=lines[0].split()[7:-1]

crunchstuff=[]

snps={}

for strain in strains:
	snps[strain]=[]

for line in lines:
	words=line.split()
		
	for x,snp in enumerate(words[7:-1]):
		snps[strains[x]].append(snp]
		
print snps
		