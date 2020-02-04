#!/usr/bin/env python

#Compares mummer snp output files.  NOTE: must be from the same reference sequence



##################
# Import modules #
##################

import string, re
import os, sys, getopt


#########
# Usage #
#########

def Usage():
	print '\nMLST_concatenator.py Usage:\n'
	print '\nMLST_concatenator.py -m <tab delimited ST list file (see below)> -o <output file name> <MLST fasta files>'
	print '\nST list file must have a header with gene names and ST as the last column:'
	print '\naroC	dnaN	hemD	hisD	purE	sucA	thrA	ST'
	print '1	1	1	1	1	1	5	1'

	print '\nCopyright Simon R Harris, Wellcome Trust Sanger Institute, UK. 2009\n'

##############################
# Get command line arguments #
##############################


def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "hm:o:", ["mlst=", "output="])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	mlstfile=""
	outfile=""
	fastafiles=[]

	for opt, arg in opts:
	
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-m", "--mlst"):
			mlstfile=arg
		elif opt in ("-o", "--output"):
			outfile=arg
	
	fastafiles=args
	
	if mlstfile=='':
		print 'Error: Missing MLST file'
		Usage()
		sys.exit()
	elif fastafiles==[]:
		print 'Error: No fasta files given'
		Usage()
		sys.exit()
	elif outfile=='':
		print 'Error: Missing output file'
		Usage()
		sys.exit()
	
	
	return mlstfile, outfile, fastafiles




if __name__ == "__main__":
	argv=sys.argv[1:]
	mlstfile, outfile, fastafiles=getOptions(argv)
	
	
	seqs={}
	for f in fastafiles:
		name=f.split('.')[0].upper()

		if seqs.has_key(name):
			print "Error: Two of your fasta files have the same prefix:", name
			sys.exit()
		else:
			seqs[name]={}
		lines=open(f,'rU').read().split(">")[1:]
		for line in lines:
			seqs[name][line.split('\n')[0].split()[0].replace(name,'')]=''.join(line.split('\n')[1:])
		

	lines=open(mlstfile, 'rU').readlines()
	
	output=open(outfile,"w")
	
	header=lines[0]
	
	genes=header.strip().split()[:-1]
	
	lines=lines[1:]
	for line in lines:
		if len(line.strip().split())<len(header.strip().split()):
			continue
		words=line.strip().split()[:-1]
		currst=line.strip().split()[-1]
		
		print >> output, ">ST"+currst
		
		for x, word in enumerate(words):
			if not seqs.has_key(genes[x].upper()):
				print "Error: No fasta for gene in list:", name
				sys.exit()
			elif not seqs[genes[x].upper()].has_key(word):
				print "Error: type", word, "missing from fasta file for gene", genes[x]
				sys.exit()	
			
			print >> output, seqs[genes[x].upper()][word]
		
		
	
	
	output.close()
	
	
	
	
	
	


