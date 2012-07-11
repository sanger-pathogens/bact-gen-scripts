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
	outdir=''
	
	for opt, arg in opts:
	
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-o", "--out"):
			outdir=arg

	inputdirs=args
	

	if inputdirs==[]:
		print 'Error: No input files selected!'
		Usage()
		sys.exit()
	if outdir[-1]!='/':
		outdir=outdir+'/'

	return inputdirs, outdir




if __name__ == "__main__":
        argv=sys.argv[1:]
        inputdirs, outdir=getOptions(argv)



pools=[]



for pool in inputdirs:
	pool=pool.split('/')[-1].split('.')[0]
	pools.append(pool)



#Running Maq and Velvet where required

for pool in pools:

	#Creating single unmap file where there are multiple
	os.system("cat "+pool+"/unmap*-*.txt > "+pool+"/unmap.txt")
	
	#Changing unmap txt file to fastq format
	unmapout=open(outdir+pool+".fastq", "w")
	lines=open(pool+"/unmap.txt", "rU").readlines()
	for line in lines:
		words=line.split()
		print >> unmapout, '@'+words[0]
		print >> unmapout, words[2]
		print >> unmapout, '+'
		print >> unmapout, words[3]
	unmapout.close()
		
	