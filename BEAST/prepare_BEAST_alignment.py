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
sys.path.extend(map(os.path.abspath, ['/nfs/pathogen/sh16_scripts/modules/']))
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
	group.add_option("-w", "--overwrite", action="store_true", dest="overwrite", help="Force overwrite", default=False)
	group.add_option("-d", "--dates", action="store", dest="dates", help="Dates csv file name (Name,date)", default="")
	group.add_option("-n", "--nons", action="store_true", dest="nons", help="Exclude sites which are constant other than Ns", default=False)
	group.add_option("-g", "--gaps", action="store_true", dest="gaps", help="Gaps (-) are real", default=False)
	group.add_option("-x", "--exclude", action="store_true", dest="exclude", help="Exclude isolates without dates from output", default=False)
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
	
	if options.dates!='' and not os.path.isfile(options.dates):
		DoError('Cannot find file '+options.dates)
	
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



def read_metadata(filehandle, name_column=1, unit_column=2, value_column=3, header=False, name_heading="", unit_heading="", value_heading="", split_value="\t"):
	
	metadata={}
	for x, line in enumerate(filehandle):
		line=line.strip()
		if x==0 and header:
			headings=line.split(split_value)
			for colcount, colhead in enumearate(headings):
				if unit_heading!="" and colhead==unit_heading:
					unit_column=colcount+1
				elif value_heading!="" and colhead==value_heading:
					value_column=colcount+1
				elif name_heading!="" and colhead==name_heading:
					name_column=colcount+1
		
		else:
			columns=line.split()
			if len(columns)<max([unit_column, value_column]):
				print 'metadata has rows without correct number of columns...skipping:'
				print line
			else:
				name=columns[name_column]
				if name in metadata:
					DoError("Repeated names in metadata")
				metadata["name"]={}
				
				unit=columns[unit_column].lower()
				if not unit in ["days", "months", "years"]:
					DoError("Time unit must be days, months or years!")
				
				try:
					value=float(columns[value_column])
				except ValueError:
					DoError("Time value column must be a float!")
	



if __name__ == "__main__":
		#argv=sys.argv[1:]
		#ref, inputfile, outfile, tabfile, align, embl, raxml, graphs, bootstrap, model, chisquare, recomb=getOptions(argv)
		
	
	(options, args)=get_user_options()
	
	check_input_validity(options, args)
	
	
	dates={}
	
	if options.dates!="":
		
		for line in open(options.dates):
			words=line.strip().split(',')
			if len(words)<2:
				continue
				
				
			try:
				dates[words[0]]=int(words[1])
			except ValueError:
				continue
	
	
	print '\nReading input alignment...'
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
	
	exccount=0
	for line in lines:
		words=line.strip().split('\n')
		if options.exclude and not words[0].split()[0] in dates:
			exccount+=1
		else:
			sequences[words[0].split()[0]]=''.join(words[1:])
	
	if exccount>0:
		print "Removed", exccount, "sequences with no dates"
	if len(sequences)==0:
		print "\nNo sequences found"
		sys.exit()
	

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
			else:
				constants[foundbases.keys()[0]]=0
				constants[foundbases.keys()[0]]+=1
	
	
	print "Done"
	print "Found", len(snplocations)-Ncount, "sites with a SNP"
	if not options.nons:
		print "Found a further ", Ncount, "sites with no SNP, but at least one N"
	print "Constant bases:"
	
	sortbase=constants.keys()
	sortbase.sort()
	
	for base in sortbase:
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
	
	print '\nOr use /nfs/pathogen/sh16_scripts/BEAST/replace_BEAST_blocks.py and provide the file', options.outfile+".patterns", "with the -p flag"
	
	output=open(options.outfile+".patterns","w")
	print >> output, ' '.join(map(str,[constants['A'], constants['C'], constants['G'], constants['T']]))
	output.close()


		
	alignment = Generic.Alignment(Gapped(IUPAC.unambiguous_dna, "-"))
	for name in snpsequence:
		 
#		if len(''.join(snpsequence[name]).replace("-","").replace("N",""))>float(len(snpsequence[name]))*(float(options.exclude)/100):
#			alignment.add_sequence(name, ''.join(snpsequence[name]))
#		else:
#			print name, "excluded from snp alignment as it is < "+str(options.exclude)+"% mapped"
		if name in dates:
			alignment.add_sequence(name+"_"+str(dates[name]), ''.join(snpsequence[name]))
		else:
			alignment.add_sequence(name, ''.join(snpsequence[name]))
			
	
	
	AlignIO.write([alignment], open(options.outfile, 'w'), "fasta")
		

	

