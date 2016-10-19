#!/usr/bin/env python

import os, sys
import cPickle

try:
	idx=open(sys.argv[1]+".chippy.index", "rU")
except:
	print "Could not open Chippy index file "+sys.argv[1]+".chippy.index"
	sys.exit()
index=cPickle.load(idx)
idx.close()

print "Index stats:"
print len(index), "species"
seqs=0
lengths=[]
for s in index:
	seqs+=len(index[s])
	#for a in index[s]:
	#	lengths.append([index[s][a]["locations"][1]-index[s][a]["locations"][0], a, s])
print seqs, "sequences"
		
