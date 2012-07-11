#!/usr/bin/env python
#/usr/bin/env python
import string
import os, sys
import pysam


#function to add extended cigar list to a plot list

def addcigarlisttoplotlist(cigarlist, contig, startingbase, plotlist):

	position=startingbase

	for x in cigarlist:
		if x[0]==0:
			for y in range(0, x[1]):
				if position<len(plotlist[contig]):
					plotlist[contig][position]=plotlist[contig][position]+1
				position=position+1
		elif x[0]==1:
			for y in range(0, x[1]):
				position=position+1
	
	return plotlist


if len(sys.argv)!=3 or "-h" in sys.argv:
	print "~sh16/scripts/bam2mappingplot.py <input bam/sam file> <output plot file name>"
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

depths={}

for x, ref in enumerate(refs):
	depths[ref]=[0]*lengths[x]


for read in samfile:

	if not read.is_unmapped:
#		cigarlist=readExtendedCIGAR(read.cigar)
#		print read.cigar, cigarlist
		depths=addcigarlisttoplotlist(read.cigar, refs[read.rname], read.pos, depths)
		


outfile=open(sys.argv[2], "w")

for ref in refs:
	print >> outfile, '\n'.join(map(str,depths[ref]))

outfile.close()

    	
		