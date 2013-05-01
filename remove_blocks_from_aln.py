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
	print '-k\t\tkeep regions in tab file (default is to remove them)'
	print '-o <file name>\toutput file name'
	print '-t <file name>\ttab file name (containing regions to keep/remove)'
	print '-r <name>\treference name (optional, but required if there are gaps in the reference sequence relative to the tab file)'
	print '-R\t\tDo not remove blocks from reference sequence (default is to remove from all sequences)'
	print '-s <char>\tSymbol to use for removed regions (default = N)'
	print '-h\t\tshow this help'
	print 'Written by Simon R Harris, Wellcome Trust Sanger Institute, UK. 2010'



##############################
# Get command line arguments #
##############################


def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "ho:a:t:kr:Rs:c", ["align=", "out=", "=tab", "keep", "reference=", "refrem", "symbol=", "cut"])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	outfile=''
	alnfile=''
	tabfile=''
	keepremove='r'
	reference=''
	refrem=True
	symbol="N"

	for opt, arg in opts:
	
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-o", "--out"):
			outfile=arg
		elif opt in ("-a", "--align"):
			alnfile=arg
		elif opt in ("-t", "--tab"):
			tabfile=arg
		elif opt in ("-k", "--keep"):
			keepremove='k'
		elif opt in ("-c", "--cut"):
			keepremove='c'
		elif opt in ("-r", "--reference"):
			reference=arg
		elif opt in ("-R", "--refrem"):
			refrem=False
		elif opt in ("-s", "--symbol"):
			symbol=arg

	
	if alnfile=='' or not os.path.isfile(alnfile):
		print 'Error: Alignment file not found!'
		Usage()
		sys.exit()
	elif tabfile=='' or not os.path.isfile(tabfile):
		print 'Error: Tab file not found!'
		Usage()
		sys.exit()
	elif outfile=='':
		print 'Error: No output file specified!'
		Usage()
		sys.exit()
	elif keepremove not in ['k','r', 'c']:
		print "What the???"
		sys.exit()
	symbol=symbol.upper()
	if symbol not in ["N", "X", "?", "-"]:
		print 'Error: Symbol must be N, X, ? or -!'
		Usage()
		sys.exit()
		
	
	
	
	return alnfile, outfile, tabfile, keepremove, reference, refrem, symbol




if __name__ == "__main__":
	argv=sys.argv[1:]
	
	alnfile, outfile, tabfile, keepremove, reference, refrem, symbol=getOptions(argv)
	
	if keepremove=='c':
		refrem=True
	
	regions=[]
	
	for line in open(tabfile, 'rU'):
		if len(line.split())>2 and line.split()[0]=="FT" and (line.split()[1].lower() in ["misc_feature", "cds", "mobile_element", "fasta_record"]):
			if len(line.split()[2].split('..'))==1:
				start=int(line.split()[2])
				end=start
				regions.append([start-1,end-1,"f"])
			elif line.split()[2][:10]=="complement":
				line=line.replace(")","").replace("complement(", "")
				
				if int(line.split()[2].split('..')[0])<int(line.split()[2].split('..')[1]):
					start=int(line.split()[2].split('..')[0])
					end=int(line.split()[2].split('..')[1])
				else:
					start=int(line.split()[2].split('..')[1])
					end=int(line.split()[2].split('..')[2])
				
				regions.append([start-1,end-1,"r"])
			else:
				if int(line.split()[2].split('..')[0])<int(line.split()[2].split('..')[1]):
					start=int(line.split()[2].split('..')[0])
					end=int(line.split()[2].split('..')[1])
				else:
					start=int(line.split()[2].split('..')[1])
					end=int(line.split()[2].split('..')[2])
				regions.append([start-1,end-1,"f"])
	
	print "Found", len(regions), "regions"
	
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
	
	for line in lines:
		words=line.strip().split('\n')
		if words[0].split()[0] in sequences:
			DoError("Found duplicate sequence name: "+words[0].split()[0])
			
		sequences[words[0].split()[0]]=''.join(words[1:])
	print "Found", len(sequences.keys()), "sequences"
	
	
	if reference!="":
		if not reference in sequences.keys():
			DoError("Cannot find reference in alignment")
		print "Adjusting region locations to alignment positions"
		reftoaln={}
		refnum=0
		print reference
		for alnnum, base in enumerate(sequences[reference]):
			if base!="-" and base!="N":
				reftoaln[refnum]=alnnum
				refnum+=1
		for region in regions:
			region[0]=reftoaln[region[0]]
			region[1]=reftoaln[region[1]]
			
	
	
	regions.sort()
	
	newsequences={}
	
	for sequence in sequences.keys():
		if not refrem and sequence==reference:
			newsequences[sequence]=sequences[sequence]
			continue
		newsequences[sequence]=''
		for x, region in enumerate(regions):
			if keepremove=='k':
				if region[2]=="f":
					newsequences[sequence]=newsequences[sequence]+sequences[sequence][region[0]:region[1]+1]
				else:
					newsequences[sequence]=newsequences[sequence]+revcomp(sequences[sequence][region[0]:region[1]+1])
			elif keepremove=='c':
				if x==0:
					newsequences[sequence]=newsequences[sequence]+sequences[sequence][:region[0]]
				else:
					newsequences[sequence]=newsequences[sequence]+sequences[sequence][regions[x-1][1]+1:region[0]]
				if x==len(regions)-1:
					newsequences[sequence]=newsequences[sequence]+sequences[sequence][region[1]+1:]
			else:
				if x==0:
					newsequences[sequence]=newsequences[sequence]+sequences[sequence][:region[0]]+symbol*(region[1]+1-region[0])
				else:
					newsequences[sequence]=newsequences[sequence]+sequences[sequence][regions[x-1][1]+1:region[0]]+symbol*(region[1]+1-region[0])
				if x==len(regions)-1:
					newsequences[sequence]=newsequences[sequence]+sequences[sequence][region[1]+1:]
	
				
	alnout=open(outfile,"w")
	for sequence in newsequences:
		print >> alnout, ">"+sequence
		print >> alnout, newsequences[sequence]
			
	alnout.close()
	
	if len(newsequences[sequence])!=len(sequences[sequence]) and keepremove!='c':
		print "ERROR: Output and input sequences of different lengths. Do you have overlapping features in your inputfile?"
	
	print "Original alignment length:", len(sequences[sequence]), "New alignment length:", len(newsequences[sequence])
	print "Done."
