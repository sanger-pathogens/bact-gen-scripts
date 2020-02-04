#!/usr/bin/env python
#/usr/bin/env python

##################
# Import modules #
##################

import string, re
import os, sys, getopt, random, math
from random import *
from Bio import SeqIO
from Bio.Seq import Seq
from optparse import OptionParser, OptionGroup



####################
# Set some globals #
####################

RAXML_DIR=""


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
	group.add_option("-r", "--reference", action="store", dest="ref", help="Name of reference strain", default="")
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "Optional")
	group.add_option("-o", "--output", action="store", dest="outfile", help="Output file prefix", default="")
	group.add_option("-a", "--align", action="store_true", dest="align", help="Produce snp alignment file (in phylip format)", default=False)
	group.add_option("-p", "--phylogeny", action="store_true", dest="raxml", help="Run phylogeny with RAxML", default=False)
	group.add_option("-m", "--model", action="store", dest="model", help="Model of evolution to use. [Default= %default]", default="GTRGAMMA", type="choice", choices=["GTRGAMMA","GTRGAMMAI", "GTRCAT", "GTRMIX", "GTRMIXI"])
	group.add_option("-b", "--bootstrap", action="store", dest="bootstrap", help="Number of bootstrap replicates (0 = do not run bootstrap). [Default= %default]", default=100, type="int", metavar="int")
	parser.add_option_group(group)
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


	if options.ref=='':
		DoError('No reference selected!')
	elif options.inputfile=='':
		DoError('No input file selected!')
	elif not os.path.isfile(options.inputfile):
		DoError('Cannot find file '+options.inputfile+'!')
	elif options.bootstrap>10000 or options.bootstrap<0:
		DoError('Number of bootstrap replicates (-b) must be between 0 (do not run bootstrap) and 10,000!')
		
	if options.outfile=='':
		options.outfile=options.ref.split("/")[-1].split(".")[0]


	while os.path.isfile(options.outfile+".aln") and options.overwrite==False:
		outopt=""
		outopt=raw_input('\nOutput files with chosen prefix already exist.\n\nWould you like to overwrite (o), choose a new output file prefix (n) or quit (Q): ')
		if outopt=='Q':
			sys.exit()
		elif outopt=="o":
			break
		elif outopt=="n":
			options.outfile=raw_input('Enter a new output file prefix: ')
		
		
	return





class SNPanalysis:
	def __init__(self, fastq='', directory='', mapped={}, runmaq='n', CDSseq=''):
		self.fastq=fastq
		self.directory=directory
		self.mapped=mapped
		self.goodlen=0
		self.CDSseq=''
		self.snpsummary={}
		self.nummapped=0
		self.percentmapped=0.0
		self.SNPs={}



if __name__ == "__main__":
		#argv=sys.argv[1:]
		#ref, inputfile, outfile, tabfile, align, embl, raxml, graphs, bootstrap, model, chisquare, recomb=getOptions(argv)
		
	
	(options, args)=get_user_options()
	
	check_input_validity(options, args)
	
	snps={}
	
	refbases={}
	bases={}#0=A 1=C 2=G 3=T
	nostates=0
	converter={}
	convertback={}
	snpbases={}
	
	snpstructs=[]
	count=0
	namesort=[]
	runpileups='n'
	
	
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
	snpstructs.append(SNPanalysis())	
	
	
	ref=options.ref

	reflen=len(sequences[sequences.keys()[0]])

	for sequence in sequences.keys():
		if len(sequences[sequence])!=reflen:
			print "\nERROR!: sequences are not all of the same length!!!\n"
			sys.exit()
		sequences[sequence]=sequences[sequence].upper().replace('N','-')

	
	#print sequences.keys()
	print "Found", len(sequences.keys()), "sequences of length", reflen
	sys.stdout.flush()
	
	if not ref in sequences.keys():
		print "Error!!!: Reference ("+ref+") not in alignment"
		sys.exit()

	
	count=0
	if ref!='':
		aln_pos_to_ref_pos={}
		ref_pos_to_aln_pos={}
		#print reflen
		refseq=sequences[ref].replace('-','')
		reflennogaps=len(refseq)
		for x, y in enumerate(sequences[ref]):
			if y!='-':
				
				aln_pos_to_ref_pos[x]=count#use aln_pos_to_ref_pos to convert alignment positions to positions in ref sequence
				ref_pos_to_aln_pos[count]=x#opposite to aln_pos_to_ref_pos. Use this to convert reference positions to alignment positions
				count=count+1
	
	snplocations=[]
	snpbases=[]
	basecounts=[]

	
	
	#Identify snps in the alignment
	
	print "\nIdentifying SNPs...",
	sys.stdout.flush()
	for x in range(reflen):
		numbases=0
		foundbases={}	
		for key in sequences.keys():
			base=sequences[key][x].upper()
			if base not in foundbases.keys() and base!='-':
				foundbases[base]=1
				numbases=numbases+1
			elif base!='-':
				foundbases[base]=foundbases[base]+1
		if numbases>1:
			snplocations.append(x)
			snpbases.append(foundbases)
			basecounts.append(numbases)
	
	print "Done"
	sys.stdout.flush()
	
	
	
	
	#print CDSposns
	#sys.exit()
	namesort.sort()
	
	#avcoverage={}
	#covcount={}
	
	allsnps={}
	refbases={}
	
	alssnpsummary={}
	
				
	
	tempsnps={}
	tempsnpsummary={'A':{'C':0, 'G':0, 'T':0}, 'C':{'A':0, 'G':0, 'T':0}, 'G':{'C':0, 'A':0, 'T':0}, 'T':{'C':0, 'G':0, 'A':0}}
	
	olddone=-1
	


	
	#print summary file and alignment file
	
	print "\nWriting output file(s)...",
	sys.stdout.flush()
	
	output=open(options.outfile+'.out','w')
	
	#if len(allsnps.keys())>1:
	#	print >> output, 'Ref_sequence\tPosition_in_ref',
	#else:
	print >> output, 'Position_in_alignment',
	if ref!='':
		 print >> output, '\tPosition_in_'+ref,
	print >> output, '\tRef_base\tSNP_base\tTotal',
	
	for name in sequences.keys():
		print >> output, '\t'+name,
	print >> output, '\n',
	
	if options.align==True:
		aln={}
		alnfile=open(options.outfile+'_snps.aln', 'w')
	
	#minquality=10
	#minqualityperdepth=4
	
	snpsequence={}
	for name in sequences.keys():
		snpsequence[name]=''
	
	
	refsequence=''
	if ref=='':
		ref=sequences.keys()[0]
	
	
	poscount=1
	
	
	for x, j in enumerate(snplocations):
		
		outstring=''
		tabstring=''
		
		outstring=outstring+str(j+1)
		snpcolour='1'
		if sequences[ref][j]!='-':
			outstring=outstring+'\t'+str(aln_pos_to_ref_pos[j]+1)
		else:
			outstring=outstring+'\t-'
			snpcolour='1'
			if options.embl!='':
				outstring=outstring+'\t-\t-\t-'
		
		
		outstring=outstring+'\t'+sequences[ref][j]
		refbase=sequences[ref][j]
		refsequence=refsequence+refbase
		
			
		snpbase=''
		snpbasecount=0
		for base in snpbases[x].keys():
			if base!=refbase and snpbases[x][base]>snpbasecount:
				snpbasecount=snpbases[x][base]
				snpbase=base
			elif base!=refbase and  snpbases[x][base]==snpbasecount:
				snpbase=snpbase+','+base
	
		outstring=outstring+'\t'+snpbase+'\t'+str(snpbasecount)# need to think how to do this

		
		numsnpbases=1
		sitesnpbases=[]
		for name in sequences.keys():
			if name!=ref:
				if sequences[name][j]!='-':# and name.SNPs[key][j][1]>8:# and name.SNPs[key][j][2]>(name.SNPs[key][j][3]*4):
					if sequences[name][j]!=sequences[ref][j]:
						
						outstring=outstring+'\t'+sequences[name][j]#
						
					else:
						outstring=outstring+'\t.'
					
				else:
					outstring=outstring+'\t-'
				
			snpsequence[name]=snpsequence[name]+sequences[name][j]
			
		
		print >> output, outstring

	
	
	if options.align==True:
		
		count=0
		
		print >> alnfile, len(snpsequence.keys()), len(refsequence)
			
		for name in snpsequence.keys():
			#print >> alnfile, name.replace('pool','').replace('.snps','')+' '*(10-len(name))+snpsequence[name]
			print >> alnfile, name.replace('pool','').replace('.snps','')+' '+snpsequence[name]	
					
		alnfile.close()	

	output.close()
	
	#sys.exit()
	if options.raxml==True:
		outlen=0
		userinput='x'
		if '/' in options.outfile:
			outlen=-1*len(options.outfile.split('/')[-1])
		if os.path.isfile('RAxML_info.'+options.outfile.split('/')[-1]):
			print '\nRAxML files with extension '+options.outfile.split('/')[-1]+' already exist!'
			while userinput not in ['y','n']:
				userinput=raw_input('Overwrite? (y/n): ')
				userinput.lower()
			if userinput=='y':
				os.system('rm RAxML_*'+options.outfile.split('/')[-1])
		if options.bootstrap>0:
			print "Running RAxML phylogeny with "+options.model+" model of evolution and "+str(options.bootstrap)+" bootstrap replicates..."
			os.system(RAXML_DIR+"RAxML -f a -x "+str(randrange(1,99999))+" -p "+str(randrange(1,99999))+" -# "+str(options.bootstrap)+" -m "+options.model+" -s "+options.outfile+"_snps.aln -n "+options.outfile.split('/')[-1])
		else:
			print "Running RAxML phylogeny with "+options.model+" model of evolution"
			os.system(RAXML_DIR+"RAxML -f d -m "+options.model+" -s "+options.outfile+"_snps.aln -n "+options.outfile.split('/')[-1])
		
	
	
	
	
	print "All finished.\n"
