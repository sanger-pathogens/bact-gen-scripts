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
for s in index:
	seqs+=len(index[s])
print seqs, "sequences"
		
