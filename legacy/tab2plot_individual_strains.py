#!/usr/bin/env python

import os, sys
sequences={}
if len(sys.argv)!=3:
	print "Usage: tab2plot.py <tab file> <output suffix>"
	sys.exit()


count={}

for line in open(sys.argv[1],"r"):
	words=line.split()
	
	if len(words)==3 and words[1].lower() in ["snp", "misc_feature", "core", "cds", "gene"]:
			
		position=words[2].replace("complement(","").replace(")","")
		if len(position.split(".."))>1:
			start=int(position.split("..")[0])
			end=int(position.split("..")[1])
		else:
			start=int(position.split("..")[0])
			end=int(position.split("..")[0])+1
		
		
		
		
	
	elif len(words)>1 and words[1].split("=")[0]=="/strains":
		strains=line.strip().split("=")[1].split(", ")

		for strain in strains:
			strain=strain.replace('"','')
			if not count.has_key(strain):
				count[strain]=[]
			
			while end>len(count[strain]):
				count[strain].append(0)
			
			for f in range(start,end):
				count[strain][f]+=1
	

for strain in count.keys():

	output=open(strain+"_"+sys.argv[2],"w")
	for f in count[strain]:
		print >> output, f
	
	output.close()
				