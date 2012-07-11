#!/usr/bin/env python

#Compares mummer snp output files.  NOTE: must be from the same reference sequence



#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------

import string, re
import os, sys, getopt, random, math


def Usage():
	print 'Protein_align_to_nuc.py Usage:'
	print 'Protein_align_to_nuc.py -n [input nucleotide fasta file] -p [input protein alignment] -o [output nucleotide alignment file name]'
	print 'Written by: Simon R Harris, Wellcome Trust Sanger Institute, UK. 2009'

	

#------------------------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------------------------


def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "hn:p:o:", ["help", "nuc=", "prot=", "out="])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	outfile=''
	protfile=''
	nucfile=''

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-n", "--nuc"):
			nucfile=arg
		elif opt in ("-p", "--prot"):
			protfile=arg
		elif opt in ("-o", "--out"):
			outfile=arg


	


	if nucfile=='' or protfile=='':
		print 'Error: Missing input file!'
		Usage()
		sys.exit()
	if outfile=='':
		print 'Error: Missing output file!'
		Usage()
		sys.exit()

	return outfile, nucfile, protfile



if __name__ == "__main__":
        argv=sys.argv[1:]
        outfile, nucfile, protfile=getOptions(argv)

		
nucsequences={}

lines=open(nucfile,"rU").read().split(">")[1:]
for line in lines:
	words=line.split('\n')
	nucsequences[words[0]]=''.join(words[1:]).replace('-','').upper()

protsequences={}

lines=open(protfile,"rU").read().split(">")[1:]
for line in lines:
	words=line.split('\n')
	protsequences[words[0]]=''.join(words[1:]).upper()

output=open(outfile,"w")

geneticcode={'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT': 'S', 'TCC': 'S','TCA': 'S','TCG': 'S', 'TAT': 'Y','TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CTT': 'L','CTC': 'L','CTA': 'L','CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}
geneticcode={'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT': 'S', 'TCC': 'S','TCA': 'S','TCG': 'S', 'TAT': 'Y','TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CTT': 'L','CTC': 'L','CTA': 'L','CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

nucalign={}

for sequence in protsequences.keys():
	nucalign[sequence]=''
	if not nucsequences.has_key(sequence):
		print "Error: Cannot find sequence "+sequence+" in nucleotide file"
	
	x=0
	for aa in protsequences[sequence]:
		if aa in ['G','A','L','M','F','W','K','Q','E','S','P','V','I','C','T','H','R','N','D','T','Y','*']:
			codon=nucsequences[sequence][x:x+3]
			if aa==geneticcode[codon]:
				nucalign[sequence]=nucalign[sequence]+codon
				x=x+3
			else:
				print "Error:", aa, geneticcode[codon]
		elif aa in ['-', '?']:
			nucalign[sequence]=nucalign[sequence]+aa*3
			
		elif aa=="X":
			nucalign[sequence]=nucalign[sequence]+aa*3
			x=x+3
		
		
		
for sequence in nucalign.keys():	
	print >> output, ">"+sequence
	print >> output, nucalign[sequence]

output.close()










