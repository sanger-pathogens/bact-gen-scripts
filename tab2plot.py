#!/usr/bin/env python

import os, sys
sequences={}
if len(sys.argv)!=3 and len(sys.argv)!=4 and len(sys.argv)!=5:
	print "Usage: tab2plot.py <tab file> <output file name> <sum/present> <optional genome length>"
	sys.exit()


if len(sys.argv)>3 and not sys.argv[3] in ["sum", "present"]:
	print "arg 3 must be sum or present"
	sys.exit()

count=[]
length=False
if len(sys.argv)==5:
	try:
		length=int(sys.argv[4])
	except StandardError:
		print "arg 4 must be int"
		sys.exit()

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
		
		while end>len(count):
			count.append(0)
		
		for f in range(start,end):
			if len(sys.argv)>3 and sys.argv[3]=="present":
				count[f]=1
			else:
				count[f]+=1

output=open(sys.argv[2],"w")
for f in count:
	print >> output, f
if length:
	for x in range(length-len(count)):
		print >> output, 0

output.close()
				
