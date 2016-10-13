#!/usr/bin/env python

import os, sys
import cPickle
try:
	f=open(sys.argv[1], "rU")
except:
	print "Could not open file "+sys.argv[1]
my_species=sys.argv[2]
try:
	idx=open(sys.argv[1]+".chippy.index", "rU")
except:
	print "Could not open Chippy index file "+sys.argv[1]+".chippy.index"
	sys.exit()
index=cPickle.load(idx)
idx.close()
if not my_species in index:
	print "Could not find", my_species, "in index file"
	sys.exit()
for my_accession in index[my_species]:
	f.seek(index[my_species][my_accession][0])
	while f.tell()<index[my_species][my_accession][1]:
		line = f.readline()
		print line.strip()

f.close()
