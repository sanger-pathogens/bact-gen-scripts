#!/usr/bin/env python
import string, re, gzip
import os, sys, getopt, random, math
from Bio import SeqIO
from Bio.Seq import Seq
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_SeqIO import *
from Si_general import *

def revcomp(sequence):
	rev=sequence[::-1]
	revcomp=''
	d={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'g':'c', 'c':'g'}
	for i in rev:
		if d.has_key(i):
			revcomp=revcomp+d[i]
		else:
			revcomp=revcomp+i
	
	return revcomp


def Usage():
	print 'remove_blocks_from_aln.py Usage:'
	print 'Removes regions from a fasta alignment based on a tab file input. You can choose to remove or keep the regions in the tab file'
	print 'remove_blocks_from_aln.py [options]'
	print 'Options:'
	print '-a <file name>\talignment file name'
	print '-k\t\tkeep only taxa in text file (default is to remove them)'
	print '-o <file name>\toutput file name'
	print '-l <file name>\tinput file name (containing list of taxa to keep/remove)'
	print '-h\t\tshow this help'
	print 'Written by Simon R Harris, Wellcome Trust Sanger Institute, UK. 2011'



##############################
# Get command line arguments #
##############################


def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "ho:a:l:k", ["align=", "out=", "list=", "keep"])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	outfile=''
	alnfile=''
	infile=''
	keepremove='r'
	

	for opt, arg in opts:
	
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-o", "--out"):
			outfile=arg
		elif opt in ("-a", "--align"):
			alnfile=arg
		elif opt in ("-l", "--list"):
			infile=arg
		elif opt in ("-k", "--keep"):
			keepremove='k'
		

	
	if alnfile=='' or not os.path.isfile(alnfile):
		print 'Error: Alignment file not found!'
		Usage()
		sys.exit()
	elif infile=='' or not os.path.isfile(infile):
		print 'Error: Input file not found!'
		Usage()
		sys.exit()
	elif outfile=='':
		print 'Error: No output file specified!'
		Usage()
		sys.exit()
	elif keepremove not in ['k','r']:
		print "What the???"
		sys.exit()
	
		
	
	
	
	return alnfile, outfile, infile, keepremove




if __name__ == "__main__":
	argv=sys.argv[1:]
	
	alnfile, outfile, infile, keepremove=getOptions(argv)
	
	
	
	taxa=[]
	
	for line in open(infile, 'rU'):
		if len(line.strip().split()[0])>0:
		       taxa.append(line.strip().split()[0])
	
	print "Found", len(taxa), "taxa in list"
	
	sequences={}
	seqorder=[]
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
	
	for line in lines:
		words=line.strip().split('\n')
		sequences[words[0].split()[0]]=''.join(words[1:])
		seqorder.append(words[0].split()[0])
	print "Found", len(sequences.keys()), "sequences"
	
	

	
	newsequences={}
	newseqorder=[]
	
	for sequence in seqorder:

		if keepremove=="k" and sequence in taxa:
			newsequences[sequence]=sequences[sequence]
			newseqorder.append(sequence)
		elif keepremove=="r" and not sequence in taxa:
			newsequences[sequence]=sequences[sequence]
			newseqorder.append(sequence)
	
	
				
	alnout=open(outfile,"w")
	for sequence in newseqorder:
		print >> alnout, ">"+sequence
		print >> alnout, newsequences[sequence]
			
	alnout.close()
	
	
	
	print "Original number of taxa:", len(sequences), "New number of taxa:", len(newsequences)
	print "Done."
