#!/usr/bin/env python

#Compares mummer snp output files.  NOTE: must be from the same reference sequence



#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------

import string, re, math
import os, sys, getopt


def Usage():
	print '\nClonalFrame2ART.py\n'
	print '\nConverts ClonalFrame output to Artemis readable graph files showing recombination and substitution probabilities for each node on the ClonalFrame tree\n'
	print '\nClonalFrame2ART.py Usage:\n'
	print '\t/nfs/team81/sh16/scripts/mauve2clonalframe.py [options]\nOptions:\n\t-i\tInput file name [ClonalFrame output file]\n\t-o\tPrefix for output file names\n\t-h\tPrint this help'
	print '\nWritten by: Simon R Harris, Wellcome Trust Sanger Institute, UK. 2009\n'



#------------------------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------------------------

def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "hi:o:", ["help", "out=", "infile="])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)

	outfile=''
	infile=''
	mauve='n'
	
	for opt, arg in opts:
		
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-o", "--out"):
			outfile=arg
		elif opt in ("-i", "--infile"):
			infile=arg
		
	
	
	
	if infile=='':
		print 'no input file selected'
		Usage()
		sys.exit()
	if outfile=='':
		outfile=infile+".out"
	
	return outfile, mauve, infile
	
	


if __name__ == "__main__":
	argv=sys.argv[1:]
	outfile, mauve, infile=getOptions(argv)
	
	print "Reading file: "+infile
	data=open(infile, "rU").read().split('#')
	
	names=[]
	
	consevents={}
	
	poly=data[7].split('\n')[1:-1]
	
	conseventsbloc=data[3].split('\n')[1:-1]
	
	for x, name in enumerate(conseventsbloc):
		names.append(str(x+1))
		consevents[str(x+1)]=conseventsbloc[x].replace('0.000000e+00','0').split()
	
	recombination={}
	substitution={}
	
	for name in names:
		recombination[name]=[]
		substitution[name]=[]
		for y in range(0,len(consevents[name]),2):
			recombination[name].append(consevents[name][y])
		for y in range(1,len(consevents[name]),2):
			substitution[name].append(consevents[name][y])
	
	alnlength=len(consevents[names[0]])
	
	for name in names:
		
		print "Printing outputfile for node "+name+" to "+outfile+name+'_CFdata.txt'
		
		output=open(outfile+name+'_CFdata.txt','w')
		
		curposn=1
		currecomb=0.0
		for count, posn in enumerate(poly):
			if posn=='0':
				continue
			posn=int(posn)
			while posn>curposn:
				print >> output, str(currecomb)+' 0'
				curposn=curposn+1

			print >> output, recombination[name][count]+' '+substitution[name][count]
			currecomb=recombination[name][count]
			curposn=curposn+1
		
		output.close()
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
