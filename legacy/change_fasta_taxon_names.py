#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math
from optparse import OptionParser, OptionGroup
sys.path.extend(map(os.path.abspath, ['/nfs/pathogen/sh16_scripts/modules/']))
from Si_general import *

	
##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-f", "--fasta", action="store", dest="fasta", help="Input fasta file", default="", metavar="FILE")
	parser.add_option("-i", "--info", action="store", dest="info", help="Input information (csv) file", default="", metavar="FILE")
	parser.add_option("-c", "--columns", action="store", dest="columnsstring", help="Column numbers to include in the new names", default="", metavar="LIST")
	parser.add_option("-o", "--output", action="store", dest="outfile", help="Name for output fasta file", default="")

	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):
	
	if len(options.columnsstring)==0:
		DoError('You need to choose at least one column (-c option)')
	
	try: options.columns=map(int,options.columnsstring.split(","))
	except ValueError:
		DoError('column numbers must be integers separated by commas')
	
	if options.fasta=='':
		DoError('No fasta file selected')
	elif not os.path.isfile(options.fasta):
		DoError('Cannot find file '+options.fasta)
	
	if options.info=='':
		DoError('No info file selected')
	elif not os.path.isfile(options.info):
		DoError('Cannot find file '+options.info)

	
	
	if len(options.outfile)==0:
		DoError('You need to give an output fasta file name (-o option)')



################
# Main program #
################		

if __name__ == "__main__":

	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)

	fastadata=open(options.fasta,'rU').read()
	infolines=open(options.info,'rU').readlines()
	for line in infolines:
		newname=''
		words=line.strip().split(',')
		
		count=1
		for column in options.columns:
			
			if len(words)<(int(column)) or len(words[int(column)-1].strip())==0:
				continue
			if count==1:
				newname=words[int(column)-1].strip().replace(' ','_')
				count=count+1
			else:
				newname=newname+'_'+words[int(column)-1].strip().replace(' ','_')
		
		if newname!='':
			while newname[-1]=='_':
				newname=newname[:-1]
			
			result = re.search(re.compile("^>"+words[0].strip()+"\\b", re.MULTILINE), fastadata)
			
			if result == None:
				print "No match for", words[0].strip()

			#print substring
			#print ">"+words[0].strip(), ">"+newname
			fastadata=re.sub(re.compile("^>"+words[0].strip()+"\\b", re.MULTILINE), ">"+newname, fastadata)
			
			#fastadata=fastadata.replace('>'+words[0].strip(), '>'+newname)
	
	output=open(options.outfile,'w')
	print >> output, fastadata
	output.close()
        
        
        
        
	