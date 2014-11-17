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
		parser.add_option("-t", "--trim", action="store", dest="trim", help="Bases to trim from each end of read [default = %default]", metavar="trim", default=0, type="int")
        
		return parser.parse_args()


################
# Main program #
################                

if __name__ == "__main__":

	(options, args) = main()
 	
 	if len(args)==0:
 		print "Error: No input fast5 files specified"
		sys.exit()
 	
	if options.fastq=="":
		print "Error: No name provided for output fastq file"
		sys.exit()
	
	if options.trim<0:
		print "Error: Trim length must be 0 or greater"
		sys.exit()
	
	output=open(options.fastq, "w")
	print "Reading fast5 files"
	nods2D=0
	nodstemplate2D=0
	nods1D=0
	fnf=0
	utr=0
	tooshort=0
	success=0
	
	fast5s=[]
	hopen=False
	for arg in args:
		if os.path.isfile(arg):
			fast5s.append(arg)
		elif os.path.isdir(arg):
			for file in os.listdir(arg):
			    if file.endswith(".fast5"):
			        fast5s.append(os.path.join(arg, file))
		else:
			fnf+=1
	
	for fast5 in fast5s:
		t=""
		if not os.path.isfile(fast5):
			fnf+=1
			continue
		try:
			hdf = h5py.File(fast5, 'r')
			hopen=True
		except StandardError:
			utr+=1
			if hopen:
				hdf.close()
				hopen=False
			continue
#		try:
#			print hdf["Analyses/Basecall_2D_000/BaseCalled_template/Fastq"]
#		except KeyError:
#			continue
#		continue
		
		try:
			fq = hdf[options.twoDdataset][()]
			t="_BaseCalled_2D"
		except StandardError:
			#print "Error: Cannot find dataset:", options.dataset, "within fast5 file:", fast5
			#print "Skipping"
			nods2D+=1
			
			try:
				fq = hdf["Analyses/Basecall_2D_000/BaseCalled_template/Fastq"][()]
				t="_BaseCalled_template"
			except StandardError:
				#print "Error: Cannot find dataset:", options.dataset, "within fast5 file:", fast5
				#print "Skipping"
				nodstemplate2D+=1
			
				try:
					fq = hdf[options.oneDdataset][()]
					t="_BaseCalled_1D"
				except StandardError:
					#print "Error: Cannot find dataset:", options.dataset, "within fast5 file:", fast5
					print "Skipping", fast5
					print hdf["Analyses/"].keys()
					nods1D+=1
					if hopen:
						hdf.close()
						hopen=False
					continue
		if hopen:
			hdf.close()
			hopen=False
		#print "Fast5 file", options.fast5, "successfully read"
		
		fqstring=fq.data.__str__()
		fqlines=fqstring.split("\n")
		if len(fqlines[1].strip())>options.trim*2:
			print >> output, fqlines[0].strip()+t
			print >> output, fqlines[1].strip()[options.trim:len(fqlines[1].strip())-options.trim]
			print >> output, fqlines[2].strip()
			print >> output, fqlines[3].strip()[options.trim:len(fqlines[3].strip())-options.trim]
		else:
			tooshort+=1
			continue
		success+=1
		
	print fnf, "files not found"
	print utr, "files were unreadable"
	print nods2D, "files did not contain 2D basecalling"
	print nodstemplate2D, "2D files did not contain template basecalling"
	print nods1D, "files did not contain 1D basecalling"
	if options.trim>0:
		print tooshort, "reads too short to trim by", options.trim, "bases"
	print success, "reads  written to", options.fastq
	print "Done."
	output.close()
	

