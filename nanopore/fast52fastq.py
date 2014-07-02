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
		usage = "usage: %prog [options] <list of fast5 files or directories containing fast5 files>"
		parser = OptionParser(usage=usage)
        
		parser.add_option("-o", "--output", action="store", dest="fastq", help="output fastq file name", type="string", metavar="fastq file", default="")
		parser.add_option("-d", "--2Ddataset", action="store", dest="twoDdataset", help="2D dataset path within fast5 file [default = %default]", metavar="path", default="/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq")
		parser.add_option("-D", "--1Ddataset", action="store", dest="oneDdataset", help="1D dataset path within fast5 file [default = %default]", metavar="path", default="/Analyses/Basecall_1D_000/BaseCalled_template/Fastq")
        
		return parser.parse_args()


################
# Main program #
################                

if __name__ == "__main__":

	(options, args) = main()
 	
 	if len(args)==0:
 		print "Error: No input fastq files specified"
		sys.exit()
 	
	if options.fastq=="":
		print "Error: No name provided for output fastq file"
		sys.exit()
	
	output=open(options.fastq, "w")
	print "Reading", len(args), "fast5 files"
	nods2D=0
	nods1D=0
	fnf=0
	utr=0
	success=0
	
	fast5s=[]
	
	for arg in args:
		if os.path.isfile(arg):
			fast5s.append(arg)
		elif os.path.isdir(arg):
			for file in os.listdir(arg):
			    if file.endswith(".fast5"):
			        fast5s.append(arg+"/"+file)
		else:
			fnf+=1
	
	for fast5 in fast5s:
		if not os.path.isfile(fast5):
			fnf+=1
			continue
		try:
			hdf = h5py.File(fast5, 'r')
		except StandardError:
			utr+=1
			continue
		try:
			fq = hdf[options.twoDdataset][()]
		except StandardError:
			#print "Error: Cannot find dataset:", options.dataset, "within fast5 file:", fast5
			#print "Skipping"
			nods2D+=1
		
			try:
				fq = hdf[options.oneDdataset][()]
			except StandardError:
				#print "Error: Cannot find dataset:", options.dataset, "within fast5 file:", fast5
				#print "Skipping"
				nods1D+=1
				continue
		
		#print "Fast5 file", options.fast5, "successfully read"
		
		success+=1
		print >> output, fq
	print fnf, "files not found"
	print utr, "files were unreadable"
	print nods2D, "files did not contain 2D basecalling"
	print nods1D, "files did not contain 1D basecalling"
	print success, "reads  written to", options.fastq
	print "Done."
	output.close()
	hdf.close()

