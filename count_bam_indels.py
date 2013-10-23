#!/usr/bin/env python
#/usr/bin/env python
import string
import os, sys
import pysam
from numpy import mean, std, median
from optparse import OptionParser


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options] <list of mfa files>"
	parser = OptionParser(usage=usage)

	parser.add_option("-r", "--reference", action="store", dest="ref", help="Reference fasta (or multifasta). Must be the one used for mapping", default="")
	parser.add_option("-b", "--bam", action="store", dest="bam", help="bam file", default="")
	parser.add_option("-o", "--output", action="store", dest="output", help="output file name", default="")
	
	
	
	return parser.parse_args()



################
# Main program #
################		

if __name__ == "__main__":

	#Get command line arguments

	(options, args) = main()
	
	if options.ref=="":
		DoError("Reference file (-r) required")
	
	if not os.path.isfile(options.ref):
		DoError("Cannot find file "+options.ref)
		
	if options.bam=="":
		DoError("Bam file (-b) required")
	
	if not os.path.isfile(options.bam):
		DoError("Cannot find file "+options.bam)
	
	#Read the reference file
	

	refseqs={}
	reforder=[]
	
	reflines=open(options.ref).read().split(">")[1:]
	for ref in reflines:
		seqlines=ref.split("\n")
		refname=seqlines[0].split()[0]
		refseqs[refname]=''.join(seqlines[1:]).upper()
		reforder.append(refname)
	
	filename=options.bam
	
	if filename.split(".")[-1]=="bam":
		samfile = pysam.Samfile( filename, "rb" )
	elif filename.split(".")[-1]=="sam":
		samfile = pysam.Samfile( filename, "r" )
	else:
		print "Not a bam file"
		sys.exit()
		
	refs=samfile.references
	lengths=samfile.lengths
	
	insertions={}
	deletions={}
	SNPs={}
	depth={}
	toremove=[]
	for ref in refseqs:
		if not ref in refs:
			print >> sys.stderr,  "Error! Reference sequences in fasta and bam do not match:",ref, "in fasta but not in bam"
#			sys.exit()
			toremove.append(ref)
	
	for x in toremove:
		del refseqs[x]
		reforder.remove(x)
				
	for x, ref in enumerate(refs):
		if not ref in refseqs:
			print >> sys.stderr, "Error! Reference sequences in fasta and bam do not match:",ref, "in bam but not in fasta"
			sys.exit()
		elif len(refseqs[ref])!=lengths[x]:
			print >> sys.stderr,  "Error! Reference sequences in fasta and bam are not the same length:", ref
			sys.exit()
		depth[ref]=[0]*lengths[x]
		deletions[ref]=[0]*lengths[x]
		insertions[ref]=[0]*lengths[x]
		SNPs[ref]=[0]*lengths[x]
			
	
	
	
	
	
	for read in samfile:

		if not read.is_unmapped:# and not read.is_reverse:
			start=read.pos
			readpos=0
			refpos=start
			readseq=read.seq.upper()
			
			for cig in read.cigar:
				if cig[0]==0:
					for x in range(0,cig[1]):
						if readseq[readpos]!=refseqs[samfile.getrname(read.rname)][refpos] and refseqs[samfile.getrname(read.rname)][refpos] in ["A", "C", "G", "T"]:
							SNPs[samfile.getrname(read.rname)][refpos]+=1
						depth[samfile.getrname(read.rname)][refpos]+=1
						readpos+=1
						refpos+=1
				elif cig[0]==1:
					insertions[samfile.getrname(read.rname)][refpos]+=1
					readpos+=cig[1]
				elif cig[0]==2:
					for x in range(0,cig[1]):
						deletions[samfile.getrname(read.rname)][refpos]+=1
						refpos+=1
				elif cig[0]==4:
					readpos+=cig[1]
				elif cig[0]==5:
					continue
				else:
					print cig
	
	samfile.close()
	
	tot=0
	print '\t'.join(["Total_position", "Contig", "Contig_position", "Base", "Depth", "SNPs", "Insertions", "Deletions"])
	if options.output!="":
		output=open(options.output, "w")
		print >> output, '\t'.join(["Total_position", "Contig", "Contig_position", "Base", "Depth", "SNPs", "Insertions", "Deletions"])
	else:
		print '\t'.join(["Total_position", "Contig", "Contig_position", "Base", "Depth", "SNPs", "Insertions", "Deletions"])
		
	for ref in reforder:
		for i, base in enumerate(refseqs[ref]):
			tot+=1
			if depth[ref][i]>0:
				if options.output!="":
					print >> output, '\t'.join(map(str,[tot, ref, i+1, refseqs[ref][i], depth[ref][i], SNPs[ref][i], insertions[ref][i], deletions[ref][i]]))
				else:
					print '\t'.join(map(str,[tot, ref, i+1, refseqs[ref][i], depth[ref][i], SNPs[ref][i], insertions[ref][i], deletions[ref][i]]))
					
		
		