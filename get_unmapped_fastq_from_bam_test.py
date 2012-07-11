#!/usr/bin/env python
#/usr/bin/env python

import string
import os, sys
from random import *
import pysam
#import numpy
#from scipy import stats

MAX_INSERT=1000
MIN_INSERT=0

#CIGAR operators: 0=match, 1=insertion, 2=deletion, 3=skipped region from reference, 4=soft clip on read (i.e. sequence still in SEQ, 5=hard clip on read (sequence not in SEQ), 6=Padding (silent deletion from padded reference - for multiple alignment)
def calculate_cigar_length(cigar_sequence):
	length=0
	
	for f in cigar_sequence:
		if f[0] in [0,1,3,4,5]:
			length+=f[1]
	
	return length
	
	
def bamline2fastq(bamline, direction, handle=""):
	if handle=="":
		print "@"+bamline.qname+"/"+direction
		print bamline.seq
		print "+"
		print bamline.qual
	else:
		print >> handle, "@"+bamline.qname+"/"+direction
		print >> handle, bamline.seq
		print >> handle, "+"
		print >> handle, bamline.qual
		

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




#open the two SAM/BAM format input files (must be sorted by name)
chars = string.ascii_letters + string.digits
tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))



if sys.argv[1].split(".")[-1]=="bam":
	#pysam.sort( "-n", sys.argv[1], tmpname )
#	os.system("samtools sort -n "+sys.argv[1]+" "+tmpname )
#	os.system("cp "+sys.argv[1]+".bai "+tmpname+".bam.bai" )
#	samfile = pysam.Samfile( tmpname+".bam", "rb" )
#	os.system("rm "+tmpname+".bam*")


	samfile = pysam.Samfile( sys.argv[1], "rb" )


elif sys.argv[1].split(".")[-1]=="sam":
	samfile = pysam.Samfile( sys.argv[1], "r" )
else:
	print "Not a sam or bam file"
	sys.exit()	



refs=samfile.references
lengths=samfile.lengths
Mapping={}
for x, ref in enumerate(refs):
	Mapping[ref]=[0]*lengths[x]

output=open(sys.argv[1].split("/")[-1].split(".")[0]+"_unmapped.fastq", "w")
	
reads_iter=samfile

addedcount=0

for count, read in enumerate(reads_iter):
	
	mate = reads_iter.next()
	
	if read.is_reverse:
		readtmp=mate
		mate=read
		read=readtmp
	
	#print read, mate
	
	#if  (not read.is_proper_pair or read.opt("XT")!=85) or (not mate.is_proper_pair or  mate.opt("XT")!=85) or read.isize>MAX_INSERT or read.isize<MIN_INSERT:
	
	if not read.is_proper_pair or not mate.is_proper_pair:

		readqual=0
		toprint=True
		for base, basequal in enumerate(read.qual):
			if ord(basequal)-33<15 or ord(mate.qual[base])<15:
				toprint=False
				#print read, mate
				break	
		
		if toprint:
			bamline2fastq(read, "1", output)
			bamline2fastq(mate, "2", output)
			addedcount+=1
	
	elif read.is_proper_pair and mate.is_proper_pair:
		addcigarlisttoplotlist(read.cigar, samfile.getrname(read.rname), read.pos, Mapping)
		addcigarlisttoplotlist(mate.cigar, samfile.getrname(mate.rname), mate.pos, Mapping)


#means={}
#SDs={}
#variances={}
#for x, contig in enumerate(Mapping.keys()):
#	print stats.normaltest(Mapping[contig])
#	print stats.describe(Mapping[contig])
#	means[contig]=numpy.mean(Mapping[contig])
#	SDs[contig]=numpy.std(Mapping[contig])
#	variances[contig]=numpy.var(Mapping[contig])
#print means, SDs, variances
#sys.exit()
#
#bad_regions=[]
#windowsize=1000
#step=1000
#
#for refnum, contig in enumerate(refs):
#	for x in range(0, len(Mapping[contig])-windowsize,step):
#		windowmean=numpy.mean(Mapping[contig][x:x+windowsize])
#		windowSD=numpy.std(Mapping[contig][x:x+windowsize])
#		windowvar=numpy.var(Mapping[contig][x:x+windowsize])
#		
#		diff_of_means=float(means[contig]-windowmean)
#		
#		otherbit=numpy.sqrt((windowvar/windowsize)+(variances[contig]/lengths[refnum]))
#		
#		t= (means[contig]-windowmean)/(windowSD/numpy.sqrt(windowsize))
#		
#		print x, x+windowsize, windowmean, windowSD, diff_of_means, otherbit, diff_of_means/otherbit, t
#
#
#sys.exit()
#
#
##bad_regions=[]
##for contig in Mapping.keys():
##
##	badmapping=False
##	distance_since_last_bad_mapping=0
##	length_of_bad_mapping=0
##	
##	for x, base in enumerate(Mapping[contig]):
##		if base<means[contig]-SDs[contig] and not badmapping:
##			badmapping=True
##			length_of_bad_mapping=1
##			distance_since_last_bad_mapping=0
##			count_greater_than_mean_since_bad_mapping=0
##			start=x
##		elif base<means[contig]-SDs[contig]:
##			length_of_bad_mapping+=1
##			distance_since_last_bad_mapping=0
##		elif badmapping:
##		
##			distance_since_last_bad_mapping+=1
##			if base>means[contig]:
##				count_greater_than_mean_since_bad_mapping+=1
##			if distance_since_last_bad_mapping>1000 and length_of_bad_mapping>=((x-distance_since_last_bad_mapping)-start)/10 and ((x-distance_since_last_bad_mapping)-start)>1000:
##				bad_regions.append([contig,start,x-distance_since_last_bad_mapping, length_of_bad_mapping])
##				badmapping=False
##				length_of_bad_mapping=0
##			elif distance_since_last_bad_mapping>1000:
##				badmapping=False
##				length_of_bad_mapping=0
##		else:
##			distance_since_last_bad_mapping+=1
#		
#print bad_regions
#testout=open("test.tab", "w")
#for region in bad_regions:
#	print>>testout, 'FT   misc_feature    '+str(region[1]+1)+".."+str(region[2]+1)
#
#testout.close()
#print addedcount, "pairs added to unmapped reads fastq file"
#output.close()