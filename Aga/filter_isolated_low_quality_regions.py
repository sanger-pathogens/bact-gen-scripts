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
	
	
	parser.add_option("-o", "--output", action="store", dest="outputfile", help="output plot file name [default= %default]", type="string", metavar="FILE", default="")
	parser.add_option("-i", "--input", action="store", dest="inputfile", help="Input plot file name", type="string", metavar="FILE", default="")
	
	return parser.parse_args()


################
# Main program #
################		

if __name__ == "__main__":

	#starttime=time.clock()

	#Get command line arguments

	(options, args) = main()
	
	filename=options.inputfile
	
	if not os.path.isfile(filename):
		print 'Cannot find file', filename
		sys.exit()
	
	if options.outputfile=="":
		print 'No output file name given'
		sys.exit()
	
	output=open(options.outputfile,"w")

	print "Removing isolated low quality regions", filename
	sys.stdout.flush()
	
	inredblock=False
	blocklines=[]
	redcount=0
	lastlinenum=0
	print >> output, "#Base\tCoverage"
	for line in open(filename, "rU"):
		line=line.strip()
		if line[0]=="#":
			continue
		else:
			words=line.split()
			
			if len(words)!=3:
				print "Unexpected number of columns in row:", line
				sys.exit()
			
			try:
				values=map(float,words)
			except ValueError:
				print "Invalid value in row"
				sys.exit()
			
			
			if values[0]!=lastlinenum+1:
				if inredblock:
					for blockline in blocklines:
						print >> output, str(blockline[0])+"\t"+str(blockline[1]+blockline[2])
				inredblock=False
				redcount=0
				blocklines=[]
			elif values[1]>0:
				blocklines.append(values)
				redcount+=1
				if redcount==100:
					inredblock=True
			elif values[2]>0:
				blocklines.append(values)
				redcount=0
			else:
				if inredblock:
					for blockline in blocklines:
						print >> output, str(blockline[0])+"\t"+str(blockline[1]+blockline[2])
				inredblock=False
				redcount=0
				blocklines=[]
				
			lastlinenum=values[0]
		
			
	if inredblock:
		for blockline in blocklines:
			print >> output, str(blockline[0])+"\t"+str(blockline[1]+blockline[2])
	else:
		print >> output, str(values[0])+"\t0"
			
				
		
		