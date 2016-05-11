#!/usr/bin/env python

#Extracts orthologs from tab files and 



#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------

import string, re
import os, sys, commands, getopt, gzip


def Usage():
	print 'removeuninfsites.py Usage:'
	print 'removeuninfsites.py -i [input alignment] -o [output file name] -t [type of uninformative sites to remove (ML, parsimony or constant, m/p/c Default=)] {-h}'
	print 'or'
	print 'removeuninfsites.py --in [input alignment] --out [output file name] --type [type of uninformative sites to remove (ML, parsimony or constant, m/p/c Default=c)] --help'
	print 'Copyright Simon R Harris, Wellcome Trust Sanger Institute, UK. 2008'



#------------------------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------------------------

def getOptions(arg):


	try:
		opts, args = getopt.getopt(argv, "hi:t:o:", ["help", "in=", "type=", "out="])
	except getopt.GetoptError:
		Usage()
		sys.exit(2)

	inputfile=''
	inftype='c'
	outfile=''

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-i", "--in"):
			inputfile=arg
		elif opt in ("-t", "--type"):
                        inftype=arg
		elif opt in ("-o", "--out"):
                        outfile=arg


	if inputfile=='':
		print 'no input fasta file selected'
		Usage()
		sys.exit()
	if inftype not in ['m','p','c']:
		print 'invalid type of uninformative site to remove selected. Must be m or p (maximum likelihood or parsimony)'
                Usage()
                sys.exit()
	if outfile=='':
		outfile=inputfile+"-uninf.aln"
	
	return inputfile, inftype, outfile




if __name__ == "__main__":
	argv=sys.argv[1:]
	inputfile, inftype, outfile=getOptions(argv)

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
		print "Sequences of different lengths! Exiting..."
		sys.exit()

#Find sites that are informative
	
infsites=[]

for i in range(0,seqlen):
	states={}
	for key in sequs.keys():
		if sequs[key][i].lower() not in ['-','x','?', 'n'] and sequs[key][i].lower() in states:
			states[sequs[key][i].lower()]+=1
		elif sequs[key][i].lower() not in ['-','x','?', 'n']:
			states[sequs[key][i].lower()]=1

	#count number of states that are in more than one taxon

	nostates=len(states.keys())
	nostatesoverone=0

	for key in states.keys():
		if states[key]>1:
			nostatesoverone=nostatesoverone+1
	
	
	if inftype=='m':
		if nostatesoverone>1 or nostates>2:
			infsites.append(i)
	elif inftype=='p':
		if nostatesoverone>1:
			infsites.append(i)
	else:
		if nostates>1:
			infsites.append(i)

#Print informative sites to output file

output=open(outfile, 'w')
for key in seqorder:
	print >> output, '>'+key
	newseq=''
	for i in infsites:
		newseq=newseq+sequs[key][i]
	print >> output, newseq
output.close()

print "Done."
