#!/usr/bin/env python3

import os, sys
from optparse import OptionParser

##########################################
# Function to Get command line arguments #
##########################################

def main():

	usage = "usage: %prog [options] <input file>"
	parser = OptionParser(usage=usage)
	parser.add_option("-s", "--separator", action="store", dest="separator", help="Separator for columns in input file [default= %default]", default=",")
	return parser.parse_args()
	
	
################
# Main program #
################                

if __name__ == "__main__":
	(options, args) = main()
	
	infile=args[0]
	if  not os.path.isfile(infile):
		print("Cannot find file", infile)
		sys.exit()
	
	for line in open(infile, "rU", encoding='ISO-8859-1'):
		words=line.strip().split(options.separator)
		for x, word in enumerate(words):
			print(x+1, word)
		sys.exit()

		
