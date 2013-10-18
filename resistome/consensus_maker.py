#!/usr/bin/env python

##################
# Import modules #
##################

import string, re, os, sys
from Bio import AlignIO

#################################################################################
# Function to return a strict consensus sequence when given an alignment object #
#################################################################################

def consensus_from_alignment(alignmentObject):
	
	consensus_sequence=""
	
	print "Creating consensus..."
	sys.stdout.flush()
	
	for x in range(alignmentObject.get_alignment_length()):

		foundbases={"A":0, "C":0,"G":0,"T":0}
		for record in alignmentObject:
			base=record.seq[x].upper()
			if base in foundbases: 
				foundbases[base]+=1
		maxbase="N"
		maxbasenum=0
		
		for base in foundbases:
			if foundbases[base]>maxbasenum:
				maxbase=base
				maxbasenum=foundbases[base]
		
		if maxbasenum>0:
			consensus_sequence=consensus_sequence+maxbase
	

	#print "100.00% complete"#Found %d SNP locations" % len(SNPlocations),
	
	
	print "100% complete"
	sys.stdout.flush()
	
	return consensus_sequence



alignment = AlignIO.read(sys.argv[1], "fasta")

consensus=consensus_from_alignment(alignment)


output=open(sys.argv[2], "w")
print >> output, ">"+'.'.join(sys.argv[1].split('/')[-1].split('.')[:-1])
print >> output, consensus
output.close()