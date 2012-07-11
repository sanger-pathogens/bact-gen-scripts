#!/usr/bin/env python
#/usr/bin/env python
import string
import os, sys
import pysam
from optparse import OptionParser
from numpy import mean, std, max


#################################
# Simple Error Printing Funtion #
#################################

def DoError(ErrorString):
	print "!!!Error:", ErrorString,"!!!"
	sys.exit()



##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-b", "--bam", action="store", dest="bam", help="input file in sam or bam format (must have suffix .sam or .bam respectively)", default="", metavar="FILE")
	parser.add_option("-i", "--inserttype", action="store", dest="inserttype", type="choice", choices=["max", "95%", "manual"], help="insert size to search for pairs mapping around start and end of contigs (choose from 'max', '95%', manual). Max = maximum insert sized observed in the bam. 95% = mean insert size in bam +2*std of insert sizes in bam, manual = set value manually with the -I option [default=%default]", default="max")
	parser.add_option("-I", "--insert", action="store", dest="insert", type="int", help="value for manual insert size to search for pairs mapping around start and end of contigs [default=%default]", default=500)
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):

	if options.bam=='':
		DoError('No bam/sam input file selected')
	elif options.bam.split('.')[-1] not in ['sam', 'bam']:
		DoError('Input file '+options.bam+' must have the suffix .sam or .bam')
	elif not os.path.isfile(options.bam):
		DoError('Cannot find file '+options.bam)
	if options.inserttype=="manual" and (options.insert<1 or options.insert>1000000):
		DoError('Maximum insert size (-i) must be between 1 and 1000000')

	
	return



########
# Main #
########

if __name__ == "__main__":
	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)
	
	#open the bam/sam file
	
	if options.bam.split(".")[-1]=="bam":
		print "Reading bam file"
		format="bam"
		try:
			samfile = pysam.Samfile( options.bam, "rb" )
		except StandardError:
			DoError('Failed to open '+options.bam+'. Is it in bam format?')
	elif options.bam.split(".")[-1]=="sam":
		print "Reading sam file"
		format="sam"
		try:
			samfile = pysam.Samfile( options.bam, "r" )
		except StandardError:
			DoError('Failed to open '+options.bam+'. Is it in sam format?')
	
	#get reference names and lengths from the sam file header
	refs=samfile.references
	lengths=samfile.lengths


	if options.inserttype=="manual":
		insert=options.insert
		
	else:
		print "Calculating insert size info"
		
	
		inserts=[]
		for read in samfile:
			if not read.is_reverse and read.is_proper_pair:
				inserts.append(read.isize)
		
		meaninsert=mean(inserts)
		stdinsert=std(inserts)
		maxinsert=max(inserts)
		
		print "\tMean insert size =", meaninsert
		print "\tStd insert size =", stdinsert
		print "\tMax insert size =", maxinsert
		
		if options.inserttype=="95%":
			insert=meaninsert+(2*stdinsert)
		else:
			insert=maxinsert
	
	print "Setting search insert size to", insert, "("+options.inserttype+")"	

	

	
	#Iterate through input file and fix appropriate lines
	for ref in refs:
		count=0
		ccount=0
		for read in samfile.fetch(ref, 0, insert):
		
			if not read.is_unmapped and not read.mate_is_unmapped and not read.is_proper_pair:
			
				if read.is_reverse:
					count+=1
			
					if not read.mate_is_reverse and (read.mpos+read.rlen)>(lengths[read.rname]-(insert-read.pos)):
						ccount+=1
					#print read.pos, read.mpos, read.is_proper_pair, read.insertsize

		print ref+":", ccount/2, "("+str((float(ccount)/count)*100)+"% of relevant mapped reads) pairs provide evidence for circularity." 
		    	
			