#!/usr/bin/env python

import os, sys
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_SeqIO import *
from Si_general import *
from Bio import SeqIO


#sequences=SeqIO.parse(open(sys.argv[0], "rU"), "fasta")

if (len(sys.argv) != 3 and len(sys.argv)!=4 ) or "-h" in sys.argv:
	print "~sh16/scripts/fix_embl.py <embl file name> <output file name> <optional alignment>"


if len(sys.argv)==4:
	sequences=read_alignment(sys.argv[3])
	refseq=sequences[0].seq
	try:
		emblrecord=open_annotation(sys.argv[1], refseq)
	except StandardError:
		DoError("Cannot open annotation file "+sys.argv[1]+" please check the format")
elif len(sys.argv)==3:
	try:
		emblrecord=open_annotation(sys.argv[1])
	except StandardError:
		DoError("Cannot open annotation file "+sys.argv[1]+" please check the format")
	
handle =open(sys.argv[2],"w")
SeqIO.write([emblrecord],handle, "gb")
handle.close()