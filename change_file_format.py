#!/usr/bin/env python
import string, re, copy
import os, sys
from Bio import SeqIO
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_general import *
from Si_SeqIO import *
from optparse import OptionParser
#Read the alignment file

##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-i", "--inputfile", action="store", dest="infile", help="Name of input file", default="", metavar="FILE")
	parser.add_option("-t", "--type", action="store", dest="outputtype", help="Output file type. (choices= %choices) [default= %default]", default="GTRGAMMA", type="choice", choices=["phylip", "fasta", "clustal", "nexus", "emboss", "stockholm", "fasta-m10", "ig"], default="fasta")
	parser.add_option("-o", "--outputfilename", action="store", dest="outfile", help="Output file name. [default= change suffix to new file format name]", default="", metavar="FILE")
	
	
	return parser.parse_args()



################################
# Check command line arguments #
################################

def check_input_validity(options, args):


#	if options.embl!="" and options.reference=='':
#		DoError('No reference selected! If you give an embl file you must specify which sequence it is linked to (i.e. the reference)')
	if options.infile=='':
		DoError('No alignment file selected!')
	elif not os.path.isfile(options.infile):
		DoError('Cannot find file '+options.infile+'!')
	if options.outfile=="":
		options.outfile='.'.join(options.infile.split(".")[-1])+"."+options.outputtype
	while os.path.isfile(options.outfile) and options.force==False:
		outopt=""
		outopt=raw_input('\nA file with the chosen output file name already exists.\n\nWould you like to overwrite (o), choose a new output file name (n) or quit (Q): ')
		if outopt=='Q':
			sys.exit()
		elif outopt=="o":
			break
		elif outopt=="n":
			options.outfile=raw_input('Enter a new output file name: ')

if __name__ == "__main__":



	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)


	try:
		alignment=read_alignment(options.infile)
	except StandardError:
		DoError("Cannot open alignment file")
	
	print "Writing fasta format"
	
	output=open(options.outfile, "w")
	
	SeqIO.write(alignment,output,options.outputtype)
	
	output.close()
	
	print "Done"