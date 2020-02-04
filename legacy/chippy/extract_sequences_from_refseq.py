#!/usr/bin/env python

import os, sys
import cPickle
try:
	f=open(sys.argv[1], "rU")
except:
	print "Could not open file "+sys.argv[1]
my_species=sys.argv[2]
toprint="a"
if len(sys.argv)>3 and sys.argv[3] in ['p','c']:
	toprint=sys.argv[3]

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
	if toprint=="c" and index[my_species][my_accession]['plasmid']==True:
		continue
	elif toprint=="p" and index[my_species][my_accession]['plasmid']==False:
		continue
	f.seek(index[my_species][my_accession]["locations"][0])
	while f.tell()<index[my_species][my_accession]["locations"][1]:
		line = f.readline()
		print line.strip()

f.close()
