#!/usr/bin/env python
import os, sys, string, numpy
from Bio.Seq import Seq
from Bio import SeqIO
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_SeqIO import *


try:
	contigs=SeqIO.parse(open(sys.argv[1]), "fasta")
except StandardError:
	print "Could not open file", sys.argv[1]
	sys.exit()

contignames=[]
newcontigs={}
for contig in contigs:
	contignames.append([contig.id, contig])
	

contignames.sort()

printcontigs=[]

for contig in contignames:
	printcontigs.append(contig[1])
	
output_handle = open(sys.argv[2], "w")
SeqIO.write(printcontigs, output_handle, "fasta")
output_handle.close()

