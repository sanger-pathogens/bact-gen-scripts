#!/usr/bin/env python
import string, re, copy
import os, sys
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import AlignInfo
from Bio.Align.Generic import Alignment
from optparse import OptionParser
#sys.path.extend(map(os.path.abspath, ['/usr/lib/python2.4/site-packages/']))
#sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/lib/python2.5/site-packages/']))

sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_general import *
from Si_SeqIO import *



GAPS_AND_MISSING=["N", "-", "?", "X"]


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options] <taxon list>"
	parser = OptionParser(usage=usage)

	parser.add_option("-a", "--alignment", action="store", dest="alignment", help="alignment file name", default="", metavar="FILE")

	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):



	if options.alignment=='':
		DoError('No alignment file selected!')
	elif not os.path.isfile(options.alignment):
		DoError('Cannot find file '+options.alignment+'!')

		
		
	return





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
		alignment=read_alignment(options.alignment, quiet=True)
	except StandardError:
		DoError("Cannot open alignment file")
	
	SNPcount=0
	taxaSNPcount=0
	
	for columnnumber in range(0,alignment.get_alignment_length()):
		bases=[]
		taxabases=[]
		count=0
		taxacount=0
		for record in alignment:
			base=record.seq[columnnumber]
			if base not in bases and base not in GAPS_AND_MISSING:
				bases.append(base)
				count+=1
			if record.id in args and base not in taxabases and base not in GAPS_AND_MISSING:
				taxabases.append(base)
				taxacount+=1
			if taxacount>1:
				break
		
		if count>1:
			SNPcount+=1
			
		if taxacount>1:
			taxaSNPcount+=1
			
	print alignment.get_alignment_length(), SNPcount, taxaSNPcount
				
				
