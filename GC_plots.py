#!/usr/bin/env python
import string, re, copy
import os, sys
from Bio import SeqIO
from Bio import AlignIO
import pysam
from Bio.SeqUtils import GC
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_SeqIO import *
from Si_general import *
#Read the alignment file
	
try:
	chromosomes=read_seq_file(sys.argv[1])
except StandardError:
	DoError("Cannot open alignment file")


GCs=[]
count=0
append=GCs.append
for chromosome in chromosomes:
	
	for x in xrange(0,len(chromosome)-56,56):
		count+=1
		append(GC(chromosome.seq[x:x+56]))
		if count>100000:
			break
print len(GCs)

GCout=open("GCs.txt", "w")
print >> GCout, ','.join(map(str,GCs))
GCout.close()




if sys.argv[2].split(".")[-1]=="bam":
	samfile = pysam.Samfile( sys.argv[2], "rb" )
elif sys.argv[2].split(".")[-1]=="sam":
	samfile = pysam.Samfile( sys.argv[2], "r" )
else:
	print "Not a bam file"
	sys.exit()
	
refs=samfile.references
lengths=samfile.lengths


fGC=[]
rGC=[]
samGC=[]
count=0
append=samGCs.append
for read in samfile:

	if not read.is_unmapped:
		count+=1
		refseq=read.seq.upper()
		append(100*((float(refseq.count("G")+refseq.count("C")))/56))
	if count>100000:
		break


GCout=open("samGCs.txt", "w")
print >> GCout, ','.join(map(str,samGC))
GCout.close()






rout=open("rfile.prog", "w")

print >> rout, "gc.file<-read.csv(file=\"GCs.txt\",header=FALSE);"
print >> rout, "gc.list<-as.matrix(gc.file);"
print >> rout, "samgc.file<-read.csv(file=\"samGCs.txt\",header=FALSE);"
print >> rout, "samgc.list<-as.matrix(samgc.file);"
#print >> rout, "gc.file<-read.csv(file=\"GCs.txt\",header=FALSE);"
#print >> rout, "gc.v<-as.vector(gc.file);"
#print >> rout, "gc.list<-t(gc.v);"
#print >> rout, "samgc.file<-read.csv(file=\"samGCs.txt\",header=FALSE);"
#print >> rout, "samgc.v<-as.vector(samgc.file);"
#print >> rout, "samgc.list<-t(samgc.v);"
print >> rout, "ks.test(gc.list,\"pnorm\");"
print >> rout, "ks.test(gc.list,samgc.list);"
	
#print >> rout, "dev.off()"
rout.close()

redraw=True
#while redraw:
os.system("R CMD BATCH rfile.prog")
#	os.system("xpdf "+options.prefix+"_MST.pdf")
#	outopt=""
#	outopt=raw_input('\nTo redraw press r. To keep this output and exit press any other key: ')
#	if outopt=='r':
#		redraw=True
#	else:
#		redraw=False







#	tree=branchlengths_to_SNP_count(tree)#, lengthtype="insertion_locations")
#	tree=support_to_node_names(tree)
#	length=get_total_tree_length(tree)
#	print "Total tree length =", length

#sys.exit()

#print time.clock()-starttime
#sys.exit()