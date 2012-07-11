#!/usr/bin/env python

#Create plot of dissimilarity between pairs of strains in an alignment



#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------

import string, re
import os, sys
from optparse import OptionParser
from Bio.Align import AlignInfo
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_general import *
from Si_SeqIO import *

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
		
	
	for x, taxon in enumerate(args):
		taxonseq=""
		for sequence in alignment:
			if sequence.name==taxon:
				taxonseq=str(sequence.seq).upper()
				break
		if taxonseq=="":
			print "Cannot find ", taxon
			continue
		
		for taxonb in args[x+1:]:
			print taxon, taxonb
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
			
			if len(taxonseq)==len(taxonbseq):
			
				output=open(taxon+"_"+taxonb+"_pairwise.plot", "w")
			
				for y, base in enumerate(taxonseq):
					if base==taxonbseq[y]:
						print >> output, 0
					elif base not in gap_and_missing and taxonbseq[y] not in gap_and_missing:
						print >> output, 1
					else:
						print >> output, 0
						
				output.close()
			else:
				print taxon, "and", taxonb, "are not the same length"
			
			
	
	