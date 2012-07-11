#!/usr/bin/env python
import string, re, copy
import os, sys
from Bio import SeqIO
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_general import *
from Si_SeqIO import *

#Read the alignment file
	
try:
	alignment=read_alignment(sys.argv[1])
except StandardError:
	DoError("Cannot open alignment file")

print "Writing fasta format"

output=open(sys.argv[2], "w")

SeqIO.write(alignment,output,"fasta")

output.close()

print "Done"
