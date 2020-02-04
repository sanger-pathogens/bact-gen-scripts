#!/usr/bin/env python

#Compares mummer snp output files.  NOTE: must be from the same reference sequence



#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------

import string, re
import os, sys, getopt, random, math


def Usage():
	print 'Concatenate_aligns.py Usage:'
	print 'Concatenate_aligns.py [output alignment file name] [input alignments]'
	print 'Copyright Simon R Harris, Wellcome Trust Sanger Institute, UK. 2008'

	

#------------------------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------------------------


def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "h", ["help"])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	outfile=''
	inputdirs=[]

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
	
	if len(args)<2:
		Usage()
		sys.exit()

	outfile=args[0]
	inputfiles=args[1:]
	


	if inputfiles==[]:
		print 'Error: No input files selected!'
		Usage()
		sys.exit()
	if outfile=='':
		outfile=ref.split('.')[0]

	return inputfiles, outfile



if __name__ == "__main__":
        argv=sys.argv[1:]
        inputfiles, outfile=getOptions(argv)

		
sequences={}

for file in inputfiles:
	lines=open(file,"rU").read().split(">")[1:]
#	print lines
	for line in lines:
		words=line.split('\n')
		if not sequences.has_key(words[0].split()[0]):
			sequences[words[0].split()[0]]=''
		sequences[words[0].split()[0]]=sequences[words[0].split()[0]]+''.join(words[1:])



seqlen=len(sequences[sequences.keys()[0]])
#print seqlen, len(sequences["Arch"])
for sequence in sequences.keys():
#	print sequence, len(sequences[sequence])
	if len(sequences[sequence])!=seqlen:
		print "Error! Concatenated sequences are of different lengths. Check corresponding sequences are present in each alignment and have the same name."
		sys.exit()

output=open(outfile,"w")

for sequence in sequences.keys():
	print >> output, ">"+sequence
	print >> output, sequences[sequence]

output.close()










