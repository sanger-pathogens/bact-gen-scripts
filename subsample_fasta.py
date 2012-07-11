#!/usr/bin/env python
import os, sys, string, numpy
from Bio.Seq import Seq
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_SeqIO import *
from optparse import OptionParser

##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	
	parser.add_option("-f", "--fasta", action="store", dest="fasta", help="Name of input fasta file", default="")
	parser.add_option("-m", "--min_contig_len", action="store", dest="minlength", help="Minimum contig length to include in stats", default=200, type="int")
	parser.add_option("-o", "--outfile", action="store", dest="outfile", help="Name for output file", default="")
	
	return parser.parse_args()



def check_input_options(options, args):
	
	
	if options.outfile=="":
		DoError("Output file (-o) required")
	
	
	if options.fasta=="":
		DoError("Fasta file (-f) required")
	
	
	if not os.path.isfile(options.fasta):
		DoError("Cannot find file "+options.fasta)
	
	if options.minlength<=0:
		print "Minimum length of contigs to include must be >0"



(options, args) = main()
check_input_options(options, args)


try:
	contigs=SeqIO.parse(open(options.fasta), "fasta")
except StandardError:
	print "Could not open file", options.fasta
	sys.exit()


output=open(options.outfile, "w")


for contig in contigs:

	
	seq=str(contig.seq)
	
	if len(seq)>=options.minlength:
		print >> output, ">"+contig.name
		print >> output, seq
	
	

output.close()

