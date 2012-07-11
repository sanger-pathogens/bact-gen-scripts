#!/usr/bin/env python
#/usr/bin/env python
import string
import os, sys
from random import *
import pysam

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
#	os.system("samtools view -S -h -b -o "+sys.argv[1].split(".")[-1]+".bam "+sys.argv[1] )
#	os.system("samtools sort -n "+sys.argv[1].split(".")[-1]+".bam "+tmpname )
#	os.system("cp "+sys.argv[1]+".bai "+tmpname+".bam.bai" )
##	os.system("cp "+sys.argv[1]+".bai "+tmpname+".bam.bai" )
#	samfile = pysam.Samfile( tmpname+".bam", "rb" )

	samfile = pysam.Samfile( sys.argv[1], "r" )
else:
	print "Not a sam or bam file"
	sys.exit()	




refs=samfile.references
lengths=samfile.lengths

reads_iter=samfile


readstoadd={}


for count, read in enumerate(reads_iter):
	
	if not read.is_proper_pair:
		
		toprint=True
		
		if read.rname==read.mrnm and not read.is_unmapped and not read.mate_is_unmapped:
			if ( not read.is_reverse and read.pos>(lengths[read.rname]-1000)) and  (read.mate_is_reverse and read.mpos<1000) or (read.is_reverse and read.pos<1000) and (not read.mate_is_reverse and read.mpos>(lengths[read.rname]-1000)):
				toprint=False

		if toprint:
			for base, basequal in enumerate(read.qual):
				if ord(basequal)-33<15:
					toprint=False
					break
		if toprint:
			if not readstoadd.has_key(read.qname):
				readstoadd[read.qname]=[['',''],['','']]
			if read.is_read1:
				readstoadd[read.qname][1]=[read.seq, read.qual]
			else:
				readstoadd[read.qname][0]=[read.seq, read.qual]


		
handle=open(sys.argv[1].split(".")[0]+"_unmapped.fastq", "w")		
addedcount=0
for readname in readstoadd.keys():		

	if readstoadd[readname][0][0]!="" and readstoadd[readname][1][0]!="":
		print >> handle, "@"+readname+"/1"
		print >> handle, readstoadd[readname][0][0]
		print >> handle, "+"
		print >> handle, readstoadd[readname][0][1]
		print >> handle, "@"+readname+"/2"
		print >> handle, readstoadd[readname][1][0]
		print >> handle, "+"
		print >> handle, readstoadd[readname][1][1]
		addedcount+=1
		
		
#		readqual=0
#		toprint=True
#		for base, basequal in enumerate(read.qual):
#			if ord(basequal)-33<15 or ord(mate.qual[base])<15:
#				toprint=False
#				#print read, mate
#				break	
#		
#		if toprint:
#			bamline2fastq(read, "1", output)
#			bamline2fastq(mate, "2", output)
#			
#			addedcount+=1
#	
#
print addedcount, "pairs added to", sys.argv[1].split(".")[0]+"_unmapped.fastq"
handle.close()