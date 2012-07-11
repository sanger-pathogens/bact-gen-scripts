#!/usr/bin/env python

#Compares mummer snp output files.  NOTE: must be from the same reference sequence



#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------

import string, re
import os, sys, getopt, random, math


def Usage():
	print 'Nuc_to_aa.py Usage:'
	print 'Nuc_to_aa.py -i [input nucleotide fasta file] -o [output protein file name]'
	print 'Copyright Simon R Harris, Wellcome Trust Sanger Institute, UK. 2008'

	

#------------------------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------------------------


def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "hi:o:", ["help", "in=", "out="])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	outfile=''
	infile=''

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-i", "--in"):
			infile=arg
		elif opt in ("-o", "--out"):
			outfile=arg


	


	if infile=='':
		print 'Error: Missing input file!'
		Usage()
		sys.exit()
	if outfile=='':
		print 'Error: Missing output file!'
		Usage()
		sys.exit()

	return outfile, infile



if __name__ == "__main__":
        argv=sys.argv[1:]
        outfile, infile=getOptions(argv)

		
nucsequences={}
nucnames=[]

lines=open(infile,"rU").read().split(">")[1:]
for line in lines:
	words=line.split('\n')
	nucsequences[words[0]]=''.join(words[1:]).replace('-','').upper()
	nucnames.append(words[0])

output=open(outfile,"w")

geneticcode={'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT': 'S', 'TCC': 'S','TCA': 'S','TCG': 'S', 'TAT': 'Y','TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CTT': 'L','CTC': 'L','CTA': 'L','CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

protsequences={}

for sequence in nucsequences.keys():
	protsequences[sequence]=''
	
	for x in range(0,len(nucsequences[sequence]),3):
		codon=nucsequences[sequence][x:x+3]
		if len(codon)<3:
			break

		if codon[0] in ["A", "C", "G", "T"] and codon[1] in ["A", "C", "G", "T"] and codon[2] in ["A", "C", "G", "T"]:
			protsequences[sequence]=protsequences[sequence]+geneticcode[codon]
		else:
			protsequences[sequence]=protsequences[sequence]+"X"
		
		

		
		
for sequence in nucnames:	
	print >> output, ">"+sequence
	print >> output, protsequences[sequence]

output.close()










