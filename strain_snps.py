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

for line in lines[1:]:
	words=line.split()
		
	for x,snp in enumerate(words[7:-1]):
		if snp!='-' and snp!=words[4]:
			snps[strains[x]].append(words[0])


for strain in snps.keys():
	output=open(strain+'_snps.tab','w')
	print >> output, 'ID   SNP'
	
	for snp in snps[strain]:
		print >> output, 'FT   snp             '+snp
	
	
	
	output.close()
		