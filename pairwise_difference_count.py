#!/usr/bin/env python

#Create plot of dissimilarity between pairs of strains in an alignment



#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------

import string, re
import os, sys
from optparse import OptionParser
from Bio.Align import AlignInfo
from modules.Si_general import *
from modules.Si_SeqIO import *

##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options] args"
	parser = OptionParser(usage=usage)

	parser.add_option("-a", "--alignment", action="store", dest="alignment", help="alignment file name", default="", metavar="FILE")
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


#	if options.embl!="" and options.reference=='':
#		DoError('No reference selected! If you give an embl file you must specify which sequence it is linked to (i.e. the reference)')
	if options.alignment=='':
		DoError('No alignment file selected!')
	elif not os.path.isfile(options.alignment):
		DoError('Cannot find file '+options.alignment+'!')


gap_and_missing=set(["-", "N", "?"])

################
# Main program #
################		

if __name__ == "__main__":

	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)
	
	
	#Read the alignment file
	
	try:
		alignment=read_alignment(options.alignment)
	except StandardError:
		DoError("Cannot open alignment file")

	seqnames=[]
	for sequence in alignment:
		if not sequence.name in seqnames:
			seqnames.append(sequence.name)
	
	print '\t'.join(["Taxon1", "Taxon2", "SNPs", "%ID", "Aligned bases"])
	for x, taxon in enumerate(seqnames):
		taxonseq=""
#		print x+1
		for sequence in alignment:
			if sequence.name==taxon:
				taxonseq=str(sequence.seq).upper()
				break
		if taxonseq=="":
			print "Cannot find ", taxon
			continue
#		print x
		for taxonb in seqnames[x+1:]:
						
			if taxon==taxonb:
				continue
			
			taxonbseq=""
			for sequence in alignment:
				if sequence.name==taxonb:
					taxonbseq=str(sequence.seq).upper()
					break
			
			if taxonbseq=="":
				print "Cannot find ", taxonb
				continue
			
			count=0
			ncount=0
			idcount=0			
			
			if len(taxonseq)==len(taxonbseq):
			
				for y, base in enumerate(taxonseq):
					
					if base not in gap_and_missing and taxonbseq[y] not in gap_and_missing:
						idcount+=1
					
					if base!=taxonbseq[y] and base not in gap_and_missing and taxonbseq[y] not in gap_and_missing:
						count+=1
						#print y+1, base, taxonbseq[y]
					elif base=="N":
						ncount+=1
						
			else:
				print taxon, "and", taxonb, "are not the same length"
#			pairwise_distances[taxonb][taxon]=count
			print '\t'.join(map(str, [taxon, taxonb, count, ((float(idcount)-float(count))/idcount)*100, idcount]))
