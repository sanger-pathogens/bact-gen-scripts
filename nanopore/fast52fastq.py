#!/usr/bin/env python

#################################
# Import some necessary modules #
#################################
import h5py
import os, sys
from optparse import OptionParser, OptionGroup



################################
# Get the command line options #
################################


def main():
        usage = "usage: %prog [options]"
        parser = OptionParser(usage=usage)
        
        parser.add_option("-i", "--input", action="store", dest="fast5", help="input fast5 file", type="string", metavar="fast5 file", default="")
        parser.add_option("-o", "--output", action="store", dest="fastq", help="output fastq file name", type="string", metavar="fastq file", default="")
	parser.add_option("-d", "--dataset", action="store", dest="dataset", help="dataset path within fast5 file [default = %default]", metavar="path", default="/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq")
        
        return parser.parse_args()


################
# Main program #
################                

if __name__ == "__main__":

	(options, args) = main()
 	
	if options.fast5=="":
		print "Error: No fast5 file provided"
		sys.exit()
	elif not os.path.isfile(options.fast5):
		print "Error: File not found:", options.fast5
		sys.exit()
	if options.fastq=="":
		print "Error: No name provided for output fastq file"
		sys.exit()
	try:
		hdf = h5py.File(options.fast5, 'r')
	except StandardError:
		print "Error: Unable to read fast5 file:", options.fast5
		sys.exit()
	try:
		fq = hdf[options.dataset][()]
	except StandardError:
		print "Error: Cannot find dataset:", options.dataset, "within fast5 file:", options.fast5
		sys.exit()
	print "Fast5 file", options.fast5, "successfully read"
	output=open(options.fastq,"w")
	print >> output, fq
	print "Output written to", options.fastq
	output.close()
	hdf.close()

