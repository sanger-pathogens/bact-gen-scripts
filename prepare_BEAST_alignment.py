#!/usr/bin/env python

#/usr/bin/env python


##################
# Import modules #
##################

import string, re
import os, sys, getopt, math
from random import *
#sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/lib/python2.5/site-packages/']))
#from scipy.stats import chi2
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, Gapped
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import Generic
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_nexus import *


from optparse import OptionParser, OptionGroup



####################
# Set some globals #
####################

RAXML_PATH=""


##########################
# Error message function #
##########################

def DoError(errorstring):
	print "\nError:", errorstring
	print "\nFor help use -h or --help\n"
	sys.exit()
	

##########################################
# Function to Get command line arguments #
##########################################


def get_user_options():

	usage = "usage: %prog [options]"
	version="%prog 1.0. Written by Simon Harris, Wellcome Trust Sanger Institute, 2010"
	parser = OptionParser(usage=usage, version=version)

	group = OptionGroup(parser, "Required")
	group.add_option("-i", "--input", action="store", dest="inputfile", help="Input file name", default="")
	group.add_option("-o", "--output", action="store", dest="outfile", help="Output file name", default="")
	group.add_option("-n", "--nons", action="store_true", dest="nons", help="Exclude sites which are constant other than Ns", default=False)
	group.add_option("-g", "--gaps", action="store_true", dest="gaps", help="Gaps (-) are real", default=False)
	parser.add_option_group(group)
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


	if options.inputfile=='':
		DoError('No input file selected')
	elif not os.path.isfile(options.inputfile):
		DoError('Cannot find file '+options.inputfile)
	
	
	if options.outfile=='':
		options.outfile=options.ref.split("/")[-1].split(".")[0]


	while os.path.isfile(options.outfile) and options.overwrite==False:
		outopt=""
		outopt=raw_input('\nOutput files with chosen prefix already exist.\n\nWould you like to overwrite (o), choose a new output file prefix (n) or quit (Q): ')
		if outopt=='Q':
			sys.exit()
		elif outopt=="o":
			break
		elif outopt=="n":
			options.outfile=raw_input('Enter a new output file prefix: ')
		
		
	return



if __name__ == "__main__":
		#argv=sys.argv[1:]
		#ref, inputfile, outfile, tabfile, align, embl, raxml, graphs, bootstrap, model, chisquare, recomb=getOptions(argv)
		
	
	(options, args)=get_user_options()
	
	check_input_validity(options, args)
	
	
	print '\nReading input alignment...',
	sys.stdout.flush()
	
	sequences={}


	currseq=''

	#Read the alignment. If it's bigger than 2Gb read it line by line. Else read it all at once (faster)
	
	try:
		open(options.inputfile, "rU")
	except IOError:
		DoError('Cannot open alignment file '+options.inputfile)
	
	if os.path.getsize(options.inputfile)<2000000000:
			lines=open(options.inputfile, "rU").read().split('>')[1:]

	else:
		lines=[]
		count=-1
		for linea in open(options.inputfile, "rU"):
			if linea[0]==">":
				count=count+1
				lines.append(linea.split()[0][1:]+'\n')
			else:	
				lines[count]=lines[count]+linea
		linesa=[]
		sequences={}
	
	
	for line in lines:
		words=line.strip().split('\n')
		sequences[words[0].split()[0]]=''.join(words[1:])


	alnlen=len(sequences[sequences.keys()[0]])

	for sequence in sequences.keys():
		if len(sequences[sequence])!=alnlen:
			print "\nERROR!: sequences are not all of the same length!!!\n"
			sys.exit()
#		if sequence!=ref:
#			sequences[sequence]=sequences[sequence].upper().replace('N','-')
#		else:
#			sequences[sequence]=sequences[sequence].upper()
		if not options.gaps:
			sequences[sequence]=sequences[sequence].upper().replace('-','N')
		else:
			sequences[sequence]=sequences[sequence].upper()
		
		#replace uncertainties in each sequence
		for x in ["R", "S", "B", "Y", "W", "D", "K", "H", "M", "V"]:
			sequences[sequence]=sequences[sequence].replace(x,"N")

	
	#print sequences.keys()
	print "Found", len(sequences.keys()), "sequences of length", alnlen
	sys.stdout.flush()

	
	
	
	snplocations=[]
	Ncount=0
	snpsequence={}
			
	#Identify snps in the alignment
	
	print "\nIdentifying SNPs...",
	sys.stdout.flush()
	
	constants={"A":0, "C":0, "G":0, "T":0}
	
	for x in range(alnlen):
		numbases=0
		foundbases={}
		
		for key in sequences.keys():
			if not key in snpsequence:
				snpsequence[key]=[]
			base=sequences[key][x].upper()
			if base not in foundbases.keys() and ((options.nons and base not in ['-', "N"]) or (not options.nons and base!="-")):
				foundbases[base]=1
				numbases=numbases+1
#			elif base!='-':
#				foundbases[base]+=1
		if numbases>1:
			if numbases==2 and "N" in foundbases.keys():
				Ncount+=1
			snplocations.append(x)
			for key in sequences.keys():
				snpsequence[key].append(sequences[key][x].upper())
		elif numbases==1:
			if foundbases.keys()[0] in constants:
				constants[foundbases.keys()[0]]+=1
	
	
	print "Done"
	print "Found", len(snplocations)-Ncount, "sites with a SNP"
	if not options.nons:
		print "Found a further ", Ncount, "sites with a no SNP, but at least one N"
	print "Constant bases:"

	for base in ["A","C","G","T"]:
		print base+":", constants[base]
	print
	sys.stdout.flush()
	
	print "Replace the patterns block in your BEAST xml file with this:"
	
	print '\t<mergePatterns id="patterns">'
	print '\t\t<patterns from="1" every="1">'
	print '\t\t\t<alignment idref="alignment"/>'
	print '\t\t</patterns>'
	print '\t\t<constantPatterns>'
	print '\t\t\t<alignment idref="alignment"/>'
	print '\t\t\t<counts>'
	print '\t\t\t\t<parameter value="', constants['A'], constants['C'], constants['G'], constants['T'],'"/>'
	print '\t\t\t</counts>'
	print '\t\t</constantPatterns>'
	print '\t</mergePatterns>'


		
	alignment = Generic.Alignment(Gapped(IUPAC.unambiguous_dna, "-"))
	for name in snpsequence:
		 
#		if len(''.join(snpsequence[name]).replace("-","").replace("N",""))>float(len(snpsequence[name]))*(float(options.exclude)/100):
#			alignment.add_sequence(name, ''.join(snpsequence[name]))
#		else:
#			print name, "excluded from snp alignment as it is < "+str(options.exclude)+"% mapped"
		
		alignment.add_sequence(name, ''.join(snpsequence[name]))
	
	
	AlignIO.write([alignment], open(options.outfile, 'w'), "fasta")
		

	

