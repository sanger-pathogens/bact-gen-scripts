#!/usr/bin/env python
import string, re, gzip
import os, sys, getopt, random, math


def Usage():
	print 'join_alns.py Usage:'
	print 'Joins a set of alignments that contain the same reference sequence'
	print 'join_alns.py [options] <input plot files>'
	print 'Options:'
	print '-r <name>\treference name'
	print '-o <filename>\toutput file name'
	print '-h\t\tshow this help'
	print 'Written by Simon R Harris, Wellcome Trust Sanger Institute, UK. 2009'



##############################
# Get command line arguments #
##############################


def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "ho:r:", ["ref=", "out="])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	outfile=''
	inputdirs=[]
	ref=''

	for opt, arg in opts:
	
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-o", "--out"):
			outfile=arg
		elif opt in ("-r", "--ref"):
			ref=arg


	inputdirs=args
	
	if inputdirs==[]:
		print 'Error: No input files selected!'
		Usage()
		sys.exit()
	elif outfile=='':
		print 'Error: No output file specified!'
		Usage()
		sys.exit()
	elif ref=='':
		print 'Error: No reference specified!'
		Usage()
		sys.exit()
		
	overwrite='n'
	
	while overwrite!='y' and os.path.isfile(outfile):
		overwrite=raw_input(outfile+' exists! Overwrite? (y/n): ')
		overwrite.lower()
		if overwrite=='n':
			outfile=raw_input('Enter a new output file name: ')
	
	
	
	return inputdirs, outfile, ref




if __name__ == "__main__":
	argv=sys.argv[1:]
	inputdirs, outfile, ref=getOptions(argv)


output=open(outfile,'w')


alnfiles=[]
for y, alnfile in enumerate(inputdirs):
	sequences={}
	currseq=''

	if os.path.getsize(alnfile)<2000000000:
		lines=open(alnfile, "rU").read().split('>')[1:]
	else:
		lines=[]
		count=-1
		for linea in open(alnfile, "rU"):
			if linea[0]==">":
				count=count+1
				lines.append(linea.split()[0][1:]+'\n')
			else:	
				lines[count]=lines[count]+linea
		linesa=[]
	
	print len(lines)
	for line in lines:
		words=line.strip().split('\n')
		sequences[words[0]]=''.join(words[1:])
	print sequences.keys()
	if not sequences.has_key(ref):
		print "error: file", alnfile, "does not contain", ref
		continue
	
	
	
	refseq=sequences[ref]
	
	if y==0:
		print >> output, ">"+ref
		print >> output, refseq.replace('-','').replace('N','')
	
	for seq in sequences.keys():
		if seq==ref:
			continue
		
		subjectnew=''
		
		for x in range(len(refseq)):
			if refseq[x]!='-' and refseq[x]!='N':
				subjectnew=subjectnew+sequences[seq][x]
		
		print >> output, ">"+seq.split('.')[0]
		print >> output, subjectnew




output.close()
		