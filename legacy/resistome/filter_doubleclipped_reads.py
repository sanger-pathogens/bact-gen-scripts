#!/usr/bin/env python

#################################
# Import some necessary modules #
#################################
import string, re
import os, sys
from optparse import OptionParser, OptionGroup
import pysam



################################
# Get the command line options #
################################


def main():
	usage = "usage: %prog [options] args"
	parser = OptionParser(usage=usage)
	
	parser.add_option("-b", "--bam", action="store", dest="bamfile", help="Input bam file", type="string", metavar="FILE", default="")
	parser.add_option("-o", "--output", action="store", dest="output", help="Output bam file", type="string", metavar="FILE", default="")
	parser.add_option("-u", "--unmapped", action="store_true", dest="unmapped", help="keep unmapped reads", default=False)
	
	return parser.parse_args()


################
# Main program #
################		

if __name__ == "__main__":

	#starttime=time.clock()

	#Get command line arguments

	(options, args) = main()
	
	filename=options.bamfile
	
	if not os.path.isfile(filename):
		print 'Cannot find file', filename
		sys.exit()
	
	if options.output=="":
		print "No output file name provided"
		sys.exit()
	

	print "Removing reads not in proper pair that are clipped at both ends", filename
	sys.stdout.flush()
	
	try:
		samfile = pysam.Samfile( filename, "rb" )
	except StandardError:
		print 'Failed to open '+filename+'. Is it in bam format?'
		sys.exit()
	
	refs=samfile.references
	lengths=samfile.lengths

	output=pysam.Samfile(options.output, mode='wb', referencenames=refs, referencelengths=lengths)
	
	unmappedcount=0
	removed=0
	duplicates=0
	readstart=0
	reads={}

	for read in samfile:
		if readstart!=read.pos:
			readstart=read.pos
			reads={}
		
		if read.is_read1:
			fr='f'
		else:
			fr='r'
		
		name=read.qname+"_"+fr
		dup=False
		if not name in reads:
			reads[name]=read
		elif read.pos==reads[name].pos and read.cigar==reads[name].cigar:
			dup=True
		
		if not dup and not read.is_unmapped and (read.is_proper_pair or read.cigar[0][0]!=4 or read.cigar[-1][0]!=4):
			output.write(read)
		elif dup:
			duplicates+=1
		elif read.is_unmapped:
			if options.unmapped:
				output.write(read)
			else:
				unmappedcount+=1
		else:
			removed+=1
	
	print "Removed", removed, "reads that were clipped at both ends and not in proper pairs"
	print "Removed", duplicates, "duplicates"
	if not options.unmapped:
		print "Removed", unmappedcount, "unmapped reads"
	
	samfile.close()
	output.close()
	