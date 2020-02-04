#!/usr/bin/env python

#Compares mummer snp output files.  NOTE: must be from the same reference sequence



#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------

import string, re, math
import os, sys, getopt


def Usage():
	print "\nUsage:"
	print "\nalign_Chlamydia_loci_from_fastas.py <fasta files>"
	print '\nWritten by: Simon R Harris, Wellcome Trust Sanger Institute, UK. 2009\n'



#------------------------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------------------------

def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "h", ["help", "ref=", "xmfa="])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	files=[]

	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
	
	files=args
	if args==[]:
		Usage()
		sys.exit()
	

	return files

def rev(sequence):
	rev=sequence[::-1]
	
	return rev


if __name__ == "__main__":
	argv=sys.argv[1:]
	files=getOptions(argv)
	
	loci={}
	for file in files:
		lines=open(file, 'rU').read().split(">")[1:]
		for line in lines:
			locus=line.split()[0].split("_")[1]
			bug=line.split()[0].split("_")[0]
			sequence=''.join(line.split('\n')[1:])
			if not loci.has_key(locus):
				loci[locus]={}
			oldbug=bug
			count=1
			while loci[locus].has_key(bug):
				count=count+1
				bug=oldbug+"_"+str(count)
			loci[locus][bug]=sequence




for locus in loci.keys():
	output=open(locus+".fna", "w")
	bugs = loci[locus].keys()
	bugs.sort()
	for bug in bugs:
		print >> output, ">"+bug
		print >> output, loci[locus][bug]
	output.close()	
	os.system("muscle -stable -in "+locus+".fna -out "+locus+".aln")
	

	
	
	
	