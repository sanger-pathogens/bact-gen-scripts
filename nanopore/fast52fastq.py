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
        usage = "usage: %prog [options] args"
        parser = OptionParser(usage=usage)
        
        parser.add_option("-o", "--output", action="store", dest="fastq", help="output fastq file name", type="string", metavar="fastq file", default="")
	parser.add_option("-d", "--dataset", action="store", dest="dataset", help="dataset path within fast5 file [default = %default]", metavar="path", default="/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq")
        
        return parser.parse_args()


################
# Main program #
################                

if __name__ == "__main__":

	(options, args) = main()
 	
	if options.fastq=="":
		print "Error: No name provided for output fastq file"
		sys.exit()
	
	output=open(options.fastq, "w")
	print "Reading", len(args), "fast5 files"
	nods=0
	fnf=0
	utr=0
	success=0
	for fast5 in args:
		if not os.path.isfile(fast5):
			fnf+=1
			continue
		try:
			hdf = h5py.File(fast5, 'r')
		except StandardError:
			utr+=1
			continue
		try:
			fq = hdf[options.dataset][()]
		except StandardError:
			#print "Error: Cannot find dataset:", options.dataset, "within fast5 file:", fast5
			#print "Skipping"
			nods+=1
			continue
		#print "Fast5 file", options.fast5, "successfully read"
		success+=1
		print >> output, fq
	print fnf, "files not found"
	print utr, "files were unreadable"
	print nods, "files did not contain 2D basecalling"
	print success, "reads  written to", options.fastq
	print "Done."
	output.close()
	hdf.close()

