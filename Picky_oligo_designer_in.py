#!/usr/bin/env python
# !/usr/bin/env python
import string
import os, sys
from Bio import Seq
from Bio import AlignIO
from optparse import OptionParser, OptionGroup
import tre

sys.path.extend(map(os.path.abspath, ['/nfs/pathogen/sh16_scripts/modules/']))
from Si_general import *
from Si_SeqIO import *
from Si_SNPs_temp import *
#from  multiprocessing import cpu_count



import time



####################
# Set some globals #
####################


PRIMER_3_LOCATION="/nfs/users/nfs_s/sh16/primer3-2.2.2-beta/src/"






##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	
	group = OptionGroup(parser, "General options")
	group.add_option("-i", "--input", action="store", dest="alignment", help="Input alignment file name", default="", metavar="FILE")
	group.add_option("-o", "--output", action="store", dest="output", help="Prefix for output files", default="")
	group.add_option("-m", "--max_size", action="store", dest="max_size", help="Maximum oligo size [default= %default]", default=70, type="int", metavar="INT")
	parser.add_option_group(group)
	
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


#	if options.processors>cpu_count():
#		DoError('You do not have '+str(options.processors)+' processors')
#
	if options.alignment=='':
		DoError('No alignment file selected')
	elif not os.path.isfile(options.alignment):
		DoError('Cannot find file '+options.alignment)
	
	if options.max_size<1 or options.max_size>1000:
		DoError('Oligo maximum size must be between 10 and 1000')
	
	if options.output=='':
		DoError('No output file prefix selected')
	elif options.output[-1]!="_":
		options.output=options.output+"_"
	return



################
# Main program #
################		

if __name__ == "__main__":
	
	starttime=time.clock()

	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)
	
	#make random name for temporary files
	chars = string.ascii_letters + string.digits
	tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))
	
	
	#Read the alignment file

	try:
		alignment=read_alignment(options.alignment)
	except StandardError:
		DoError("Cannot open alignment file")

	consensus=consensus_from_alignment(alignment, unknowns=["N"]).replace("-", "N")
	
	consensus_chunks=consensus.split("N")
	consensusout=open(options.output+"consensus.fasta","w")
	print >> consensusout, ">consensus"
	print >> consensusout, consensus
	consensusout.close()
	
	pickyout=open(options.output+"picky_input.seq","w")
	conversion=open(options.output+"picky_backconverter.txt","w")
	print >> conversion, len(consensus)
	count=0
	position=1
	for chunk in consensus_chunks:
		if len(chunk)>options.max_size:
			for x in range(0,len(chunk),options.max_size):
				if x+(options.max_size*2)<len(chunk):
					print >> pickyout, ">Picky"+str(count)
					print >> pickyout, chunk[x:x+(options.max_size*2)]
					print >> conversion, "Picky"+str(count), str(position+x), str(position+x+(options.max_size*2))
					count+=1
				elif len(chunk)>x+options.max_size:
					print >> pickyout, ">Picky"+str(count)
					print >> pickyout, chunk[x:]
					print >> conversion, "Picky"+str(count), str(position+x), str(position+len(chunk))
					count+=1
				else:
					continue
				
				
			count+=1
		position+=len(chunk)+1
	pickyout.close()
	conversion.close()
	print count+1, "sequences created for Picky input"
	print "Open "+options.output+"picky_input.seq in Picky and identify your oligos"
	print "Then run Picky_oligo_designer_out.py with the output from Picky and the conversion file "+options.output+"picky_backconverter.txt"
	sys.exit()
	
	
