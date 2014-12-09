#!/usr/bin/env python
import string, re, gzip
import os, sys
from Bio import SeqIO
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_SeqIO import *
from optparse import OptionParser


def DoError(message):
	print "Error:", message
	sys.exit()

##########################################
# Function to Get command line arguments #
##########################################


def main():

        usage = "usage: %prog [options]"
        parser = OptionParser(usage=usage)
        parser.add_option("-m", "--map", action="store", dest="map", help="map file name", default="")
        parser.add_option("-l", "--list", action="store", dest="baselist", help="File containing list of bases to remove", default="tmp_bases")
	parser.add_option("-o", "--output", action="store", dest="output", help="output map file name", default="")
	
        return parser.parse_args()

################
# Main program #
################                

if __name__ == "__main__":

        #Get command line arguments
        
        (options, args) = main()

        #Do some checking of the input files
        if options.baselist=="":
                DoError("No list file specified")
		
	if not os.path.isfile(options.baselist):
                DoError('Cannot find file '+options.baselist)
	try:
		bases=[]
		for line in open(options.baselist, "rU"):
			bases.append(int(line.strip().split()[0]))
	except StandardError:
		DoError("Failed to read list file")
	
	
	
	if options.map=="":
                DoError("No map file specified")
		
	if not os.path.isfile(options.map):
                DoError('Cannot find file '+options.map)
	
	output=open(options.output, "w")
	for x, line in enumerate(open(options.map, "rU")):
		if not x in bases:
			print >> output, line.strip()
	output.close()
