#!/usr/bin/env python

#Extracts orthologs from tab files and 



#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------

import string, re
import os, sys, commands, getopt, gzip


def Usage():
	print 'extractindels.py Usage:'
	print 'extractindels.py -i [input alignment] -o [output file name] {-h}'
	print 'or'
	print 'extractindels.py --in [input alignment] --out [output file name] --help'
	print 'Copyright Simon R Harris, Wellcome Trust Sanger Institute, UK. 2008'



#------------------------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------------------------

def getOptions(arg):


	try:
		opts, args = getopt.getopt(argv, "hi:o:", ["help", "in=", "out="])
	except getopt.GetoptError:
		Usage()
		sys.exit(2)

	inputfile=''
	outfile=''

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-i", "--in"):
			inputfile=arg
                elif opt in ("-o", "--out"):
			outfile=arg


	if inputfile=='':
		print 'no input fasta file selected'
		Usage()
		sys.exit()
	if outfile=='':
		outfile=inputfile+"_indels.aln"
	
	return inputfile, outfile




if __name__ == "__main__":
	argv=sys.argv[1:]
	inputfile, outfile=getOptions(argv)

print "Removing uninformative sites from "+inputfile

#Open the input file

lines=open(inputfile, "rU").readlines()

sequs={}
currseq=''
seqorder=[]

#Read in each sequence

for line in lines:
	if len(line)>0 and line[0]==">":
		sequs[line.strip().split()[0][1:]]=''
		currseq=line.strip().split()[0][1:]
		seqorder.append(line.strip().split()[0][1:])
	elif len(line)>0:
		sequs[currseq]=sequs[currseq]+line.strip()

#Check all sequences are the same length

seqlen=len(sequs[currseq])

for key in sequs.keys():
	if len(sequs[key])!=seqlen:
		print "Sequences of different lengths! Exciting..."
		sys.exit()

#Find sites that contain a gaps
gapsites=[]
indel=0

for i in range(0,seqlen):
	states={}
	gapsite='n'
	contindel='n'
	for key in sequs.keys():
		if sequs[key][i]=='-':
			gapsite='y'
			if i>0 and sequs[key][i-1]=='-':
				contindel='y'
	if contindel=='n':
		indel=indel+1
	if gapsite=='y':
		gapsites.append([i,indel])

#Print sites with gaps to output file

output=open(outfile, 'w')
for key in seqorder:
	print >> output, '>'+key
	newseq=''
	lastindel=-1
	for i in gapsites:
		if i[1]!=lastindel:
			lastindel=i[1]
			newseq=newseq+'-'
		newseq=newseq+sequs[key][i[0]]
	print >> output, newseq
output.close()

print "Done."
