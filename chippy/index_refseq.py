#!/usr/bin/env python

import os, sys
import cPickle

index = {}
curspecies=""
curid=""
offset=0
line_number=0

with open(sys.argv[1], "rU") as f:
	line=f.readline()
	while line:
		words=line.strip().split()
		if len(words[0])>0 and words[0][0]==">":
			species = ' '.join(words[1:3])
			accession=words[0][1:]
			if curspecies!="":
				index[curspecies][curid].append(offset)
			if not species in index:
				index[species]={}
			index[species][accession]=[offset]
				
			curspecies=species
			curid=accession
		offset=f.tell()
		line=f.readline()
	if curspecies!="":
		index[curspecies][curid].append(offset)

f.close()

idx=open(sys.argv[1]+".chippy.index", "w")
cPickle.dump(index, idx)
idx.close()

