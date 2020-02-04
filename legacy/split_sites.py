#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math


if (len (sys.argv)!=3) or '-h' in sys.argv[1:]:
	print "split_sites.py <phylip alignment file> <partitions file>"
	sys.exit()


align=sys.argv[1]
partitions=sys.argv[2]

lines=open(partitions, "rU").readlines()

partitions={}
for line in lines:
	name=line.split("=")[0].strip().split()[1]
	partitions[name]=line.split("=")[1].strip().split()


strains=[]
sequences={}

lines=open(align, "rU").readlines()

numstrains=int(lines[0].split()[0])

straincount=0
secondset='n'
extras=0

for line in lines[1:]:
	if len(line.split())==0:
		continue
	
	
	if straincount>=numstrains:
		secondset='y'
	
	if secondset=='n':
		strain=line.split()[0]
	else:
		if straincount==numstrains:
			straincount=0
		strain=strains[straincount]
	straincount=straincount+1
	print strain
	
	if secondset=='n':
		strains.append(strain)
		sequences[strain]=''.join(line.strip().split()[1:])
	else:
		sequences[strain]=sequences[strain]+''.join(line.strip().split())
	
for partition in partitions.keys():
	output=open(align.split('/')[-1].split('.')[0]+"_"+partition+".aln", "w")
	print >> output, len(strains), len(partitions[partition])
	
	for strain in strains:
		seqstring=strain+' '
		for site in partitions[partition]:
			seqstring=seqstring+sequences[strain][int(site.replace(',',''))-1]
	
		print >> output, seqstring
	output.close()
			







