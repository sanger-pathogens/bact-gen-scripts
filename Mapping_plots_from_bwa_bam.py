#!/usr/bin/env python
#/usr/bin/env python
import string
import os, sys
from random import *
import pysam





#function to add extended cigar list to a plot list

def addcigarlisttoplotlist(cigarlist, contig, startingbase, plotlist):

	position=startingbase-1

	for x in cigarlist:
		if x[0]==0:
			for y in range(0, x[1]):
				if position<len(plotlist[contig]):
					plotlist[contig][position]=plotlist[contig][position]+1
				position=position+1
		elif x[0]==2:
			for y in range(0, x[1]):
				position=position+1
	
	return plotlist
	
	

	
#function to print a list of plot lists to file			

def print_multiple_plots(contiglist, plotlistlist, outputfilename):
	output=open(outputfilename, "w")
	for contig in contiglist:
		for x in range(0,len(plotlistlist[0][contig])):
			for plotlist in plotlistlist:
				print >> output, plotlist[contig][x],
			print >> output, "\n",
	output.close()



if sys.argv[1].split(".")[-1]=="bam":
	samfile = pysam.Samfile( sys.argv[1], "rb" )
elif sys.argv[1].split(".")[-1]=="sam":
	samfile = pysam.Samfile( sys.argv[1], "r" )
else:
	print "Not a bam file"
	sys.exit()



refs=samfile.references
lengths=samfile.lengths
UniqueMapping={}
NonUniqueMapping={}
for x, ref in enumerate(refs):
	UniqueMapping[ref]=[0]*lengths[x]
	NonUniqueMapping[ref]=[0]*lengths[x]


for read in samfile:

    if not read.is_unmapped:
#	print read.opt("XT")
        if read.opt("XT")==85:
		UniqueMapping=addcigarlisttoplotlist(read.cigar, samfile.getrname(read.rname), read.pos, UniqueMapping)
#		print read.pos, samfile.getrname(read.rname)
		
	elif read.opt("XT")==82:
		NonUniqueMapping=addcigarlisttoplotlist(read.cigar, samfile.getrname(read.rname), read.pos, NonUniqueMapping)
		#print read.pos, read.cigar, samfile.getrname(read.rname)

	#elif read.opt("XT")==77:
		#print read.pos, samfile.getrname(read.rname)

	
print_multiple_plots(refs, [UniqueMapping, NonUniqueMapping], sys.argv[2])
