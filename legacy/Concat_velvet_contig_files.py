#!/usr/bin/env python

#Compares mummer snp output files.  NOTE: must be from the same reference sequence



#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------

import string, re
import os, sys, getopt, random, math


def Usage():
	print 'compare_mummer_out.py Usage:'
	print 'compare_mummer_out.py -r=[reference sequence] [input alignment(s)] > [output file pool] {-h}'
	print 'Copyright Simon R Harris, Wellcome Trust Sanger Institute, UK. 2008'

	

#------------------------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------------------------


def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "ho:", ["help", "out="])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)

	inputdirs=[]
	outfile=''
	
	for opt, arg in opts:
	
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-o", "--out"):
			outfile=arg

	inputdirs=args
	

	if inputdirs==[]:
		print 'Error: No input files selected!'
		Usage()
		sys.exit()
	if outfile=='':
		outfile="Velvet_contigs.fasta"

	return inputdirs, outfile




if __name__ == "__main__":
        argv=sys.argv[1:]
        inputdirs, outdir=getOptions(argv)



pools=[]



for pool in inputdirs:
	pool=pool.split('/')[-1].split('.')[0]
	pools.append(pool)



#Running Maq and Velvet where required
concat=open("Velvet_contigs.fasta", "w")
for pool in pools:
	
	#Changing unmap txt file to fastq format
	
	if not os.path.isfile(pool+"/velvet/contigs.fa"):
		continue
	print pool
	lines=open(pool+"/velvet/contigs.fa", "rU").readlines()
	for line in lines:
		line=line.strip()
		if len(line)>0 and line[0]==">":
			print >> concat, ">"+pool+"_"+line[1:]
		else:
			print >> concat, line
concat.close()
		
	