#!/usr/bin/env python

import os, sys
import cPickle
try:
	f=open(sys.argv[1], "rU")
except:
	print "Could not open file "+sys.argv[1]

try:
	idx=open(sys.argv[1]+".chippy.index", "rU")
except:
	print "Could not open Chippy index file "+sys.argv[1]+".chippy.index"
	sys.exit()
index=cPickle.load(idx)
idx.close()
for accession in sys.argv[2:]:
	found=False
	for species in index:
		if accession in index[species]:
			print accession, species
			found=True
			break
	if found==False:
		print accession, "Not found"
