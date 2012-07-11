#!/usr/bin/env python
#/usr/bin/env python
import string
import os, sys
import pysam



if len(sys.argv)!=2 or "-h" in sys.argv:
	print "~sh16/scripts/bam2tradisoutput.py <input bam/sam file>"
	sys.exit()


if sys.argv[1].split(".")[-1]=="bam":
	samfile = pysam.Samfile( sys.argv[1], "rb" )
elif sys.argv[1].split(".")[-1]=="sam":
	samfile = pysam.Samfile( sys.argv[1], "r" )
else:
	print "Not a bam file"
	sys.exit()
	
refs=samfile.references
lengths=samfile.lengths



print '\t'.join(["readname", "refname", "orientation", "position"])

for read in samfile:

	if not read.is_unmapped:
	
		readname=read.qname
		
		refname=refs[read.rname]
		
		if read.is_reverse:
			orientation="-"
			position=read.pos+read.qend+1
		else:
			orientation="+"
			position=read.pos+1
		
	print '\t'.join(map(str,[readname, refname, orientation, position]))




    	
		