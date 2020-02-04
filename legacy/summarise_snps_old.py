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
from Si_SeqIO import *
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/lib/python2.7/site-packages/fisher-0.1.4-py2.7-linux-x86_64.egg']))
#from scipy import stats
import fisher
#from guppy import hpy


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
	group.add_option("-r", "--reference", action="store", dest="ref", help="Name of reference strain", default="")
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "Optional")
	group.add_option("-e", "--embl", action="store", dest="embl", help="Embl/genbank annotation file for reference strain (for dN/dS etc.)", default="")
	group.add_option("-o", "--output", action="store", dest="outfile", help="Output file prefix", default="")
	group.add_option("-t", "--tab", action="store_true", dest="tabfile", help="Produce tab files of snp locations", default=False)
	group.add_option("-a", "--align", action="store_true", dest="align", help="Produce snp alignment file (in phylip format)", default=False)
	group.add_option("-E", "--exclude", action="store", dest="exclude", help="Exclude sequences from SNP alignment if they are less thanINT% mapped [Default= %default]", default=50, type="float", metavar="int")
	group.add_option("-p", "--phylogeny", action="store_true", dest="raxml", help="Run phylogeny with RAxML", default=False)
	group.add_option("-v", "--asrv", action="store", dest="asrv", help="Method of correction for among site rate variation (optional) [Choices = GAMMA, CAT, CAT_GAMMA, MIX] [Default = %default]", default="GAMMA", type="choice", choices=["GAMMA","CAT", "CAT_GAMMA", "MIX"])
	group.add_option("-I", "--pinvar", action="store_true", dest="pinvar", help="Use correction for proportion of invariant sites", default=False)
	group.add_option("-F", "--frequencies", action="store_true", dest="f", help="Use empirical base frequencies (protein models only)", default=False)
	group.add_option("-w", "--overwrite", action="store_true", dest="overwrite", help="Overwrite old RAxML files without warning", default=False)
	group.add_option("-b", "--bootstrap", action="store", dest="bootstrap", help="Number of bootstrap replicates (0 = do not run bootstrap). [Default= %default]", default=100, type="int", metavar="int")
	group.add_option("-f", "--fast", action="store_true", dest="fast", help="Use fast bootstrap method", default=False)
	group.add_option("-l", "--LSF", action="store_true", dest="LSF", help="Run bootstrap replicates over LSF (does not apply with fast bootstrap method)", default=False)
	group.add_option("-g", "--gaps", action="store_true", dest="gaps", help="Gaps (-) are real", default=False)
	group.add_option("-G", "--genetic_code", action="store", dest="genetic_code_number", help="Genetic code number to use. Choose from 1: Standard, 4: Mycoplasma. Default is 1", default=False, type="choice", choices=[1, 4])
	parser.add_option_group(group)
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


	if options.ref=='':
		DoError('No reference selected')
	elif options.inputfile=='':
		DoError('No input file selected')
	elif not os.path.isfile(options.inputfile):
		DoError('Cannot find file '+options.inputfile)
	elif options.embl!='' and not os.path.isfile(options.embl):
		DoError('Cannot find file '+options.embl)
	elif options.bootstrap>10000 or options.bootstrap<0:
		DoError('Number of bootstrap replicates (-b) must be between 0 (do not run bootstrap) and 10,000')
	elif options.exclude>=100 or options.exclude<0:
		DoError('Exclude percentage must be an float between 0 and 100')
	
	
	if options.outfile=='':
		options.outfile=options.ref.split("/")[-1].split(".")[0]
	
	if options.raxml and not options.align:
		options.align=True

	while os.path.isfile(options.outfile+".out") and options.overwrite==False:
		outopt=""
		outopt=raw_input('\nOutput files with chosen prefix already exist.\n\nWould you like to overwrite (o), choose a new output file prefix (n) or quit (Q): ')
		if outopt=='Q':
			sys.exit()
		elif outopt=="o":
			break
		elif outopt=="n":
			options.outfile=raw_input('Enter a new output file prefix: ')
		
		
	return



def split_len(seq, length):
    return [seq[i:i+length] for i in range(0, len(seq), length)]

#############################################
# Function to create and embl style ID line #
#############################################

def embl_style_id(sequence):
	
	emblstring="ID   XX000001; SV 1; circular; genomic DNA; STD; PRO; "+str(len(sequence.strip()))+" BP.\nXX"

	return emblstring



#############################################
# Function to create and embl style ID line #
#############################################

def embl_style_header(sequence):

	emblstring="FH   Key             Location/Qualifiers"
	emblstring=emblstring+"\nFH"
	emblstring=emblstring+"\nFT   source          1.."+str(len(sequence.strip()))

	return emblstring


####################################################
# Function to create an embl style sequence format #
####################################################

def embl_style_sequence(sequence):
	
	emblstring="XX\n"
	sequence=sequence.strip().lower()
	basehash={"a":0,"c":0,"g":0,"t":0,"Other":0}
	
	for x in sequence:
		if basehash.has_key(x):
			basehash[x]=basehash[x]+1
		else:
			basehash["Other"]=basehash["Other"]+1
	
	emblstring=emblstring+"SQ   Sequence "+str(len(sequence))+" BP; "+str(basehash["a"])+" A; "+str(basehash["c"])+" C; "+str(basehash["g"])+" G; "+str(basehash["t"])+" T; "+str(basehash["Other"])+" other;\n"
	
	currentposition=0
	for x in range(0,len(sequence), 60):
		emblstring=emblstring+"     "
		charactersadded=0
		
		for y in range(0,60,10):
			if x+y+10<len(sequence):
				emblstring=emblstring+sequence[x+y:x+y+10]+" "
				charactersadded=charactersadded+11
				currentposition=currentposition+10
			elif x+y<len(sequence):
				emblstring=emblstring+sequence[x+y:]+" "
				charactersadded=charactersadded+(1+len(sequence))-(x+y)
				currentposition=currentposition+len(sequence)-(x+y)
		emblstring=emblstring+" "*(75-(charactersadded+len(str(currentposition))))+str(currentposition)+"\n"
	
	emblstring=emblstring+"//"

	return emblstring


#def embl_style_sequence(sequence):
#	
#	emblstring="XX\n"
#	sequence=sequence.strip().upper()
#	basehash={"A":0,"C":0,"G":0,"T":0,"Other":0}
#	
#	for x in sequence:
#		if basehash.has_key(x):
#			basehash[x]=basehash[x]+1
#		else:
#			basehash["Other"]=basehash["Other"]+1
#	
#	
#	split_seq=split_len(sequence, 10)
#	
#	emblstring=emblstring+"SQ   Sequence "+str(len(sequence))+" BP; "+str(basehash["A"])+" A; "+str(basehash["C"])+" C; "+str(basehash["G"])+" G; "+str(basehash["T"])+" T; "+str(basehash["Other"])+" other;\n"
#	
#	currentposition=0
#	for x in range(0,len(split_seq), 6):
#		emblstring=emblstring+"     "
#		charactersadded=0
#		
#		for y, subseq in enumerate(split_seq[x:x+6]):
#			if len(subseq)==10:
#				emblstring=emblstring+subseq+" "
#				charactersadded=charactersadded+11
#				currentposition=currentposition+10
#			else:
#				emblstring=emblstring+subseq+" "
#				charactersadded=charactersadded+(1+len(sequence))-(x+y)
#				currentposition=currentposition+len(sequence)-(x+y)
#		emblstring=emblstring+" "*(68-charactersadded)+str(currentposition)+"\n"
#	
#	emblstring=emblstring+"\\\\"
#
#	testout=open("test.embl","w")
#	print >> testout, emblstring


##################################################################
# Function to count the minimum number of changes between codons #
##################################################################

def countcodonchanges(codon, SNPcodon, geneticcode, sd, nd, loopsd=0, loopnd=0, pathcount=0):
	
	for x in range(3):
		if codon[x]!=SNPcodon[x]:
			newSNPcodon=SNPcodon[:x]+codon[x]+SNPcodon[x+1:]
			
			#print  SNPcodon, newSNPcodon, geneticcode[SNPcodon], geneticcode[newSNPcodon]
			
			if geneticcode[newSNPcodon]=='*':
				continue
			elif geneticcode[SNPcodon]==geneticcode[newSNPcodon]:
				newloopnd=loopnd
				newloopsd=loopsd+1
			else:
				newloopnd=loopnd+1
				newloopsd=loopsd
			
			
			#print SNPcodon, newSNPcodon, codon, sd, nd, newloopsd, newloopnd, pathcount
			if newSNPcodon!=codon:
				sd, nd, pathcount=countcodonchanges(codon, newSNPcodon, geneticcode, sd, nd, newloopsd, newloopnd, pathcount)
			
			else:
				sd=sd+newloopsd
				nd=nd+newloopnd
				pathcount=pathcount+1
			
			
	return sd, nd, pathcount


#####################################################
# Function to calculate dN/dS between two sequences #
#####################################################


def dnbyds(CDS, SNPseq, CDSbasenumbers, genetic_code_number=1):
	
	#standard genetic code
	geneticcode_1={'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT': 'S', 'TCC': 'S','TCA': 'S','TCG': 'S', 'TAT': 'Y','TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CTT': 'L','CTC': 'L','CTA': 'L','CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}
	#mycoplasma genetic code
	geneticcode_4={'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT': 'S', 'TCC': 'S','TCA': 'S','TCG': 'S', 'TAT': 'Y','TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': 'W', 'TGG': 'W', 'CTT': 'L','CTC': 'L','CTA': 'L','CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}
	
	
	
	geneticcodes=[geneticcode_1,geneticcode_1,geneticcode_1,geneticcode_4]
	geneticcode=geneticcodes[genetic_code_number-1]
	
	codonsynonyms={}
	
	for codon in geneticcode.keys():
		thiscodon={}
		codonsynonyms[codon]=0.0
		
		for x in range(3):
			numsyn=0.0
			numnotstop=3
			for y in ['A', 'C', 'G', 'T']:
				if codon[x]!=y:
					newcodon=codon[:x]+y+codon[x+1:]
					if geneticcode[newcodon]==geneticcode[codon]:
						numsyn=numsyn+1
					elif geneticcode[newcodon]=='*':
						numnotstop=numnotstop-1
			codonsynonyms[codon]=codonsynonyms[codon]+(numsyn/numnotstop)
	
	S=0.0
	N=0.0
	S1=0.0
	N1=0.0
	S2=0.0
	N2=0.0
	Sd=0.0
	Nd=0.0
	pS=0.0
	pN=0.0
	gapcount=0
	numcodons=0
	varianceS=0.0
	varianceN=0.0
	z=0.0
	dN=0.0
	dS=0.0
	SNPtype={}
	AAfromtype={}
	AAtotype={}
	Fisher=-0
	
	if len(CDS)!=len(SNPseq):
		print "Error: sequences must be the same length to calculate dN/dS!"
		[N, S, dN, dS, pN, pS, varianceS, varianceN, z, (len(CDS)-gapcount), Nd, Sd, Fisher], SNPtype, AAfromtype, AAtotype
	
	
	
	
	for x in range(0,len(CDS),3):
		if len(CDS)<x+3:
			continue
		numcodons=numcodons+1
		codon=CDS[x:x+3]
		SNPcodon=SNPseq[x:x+3]
		codonposn=[CDSbasenumbers[x], CDSbasenumbers[x+1], CDSbasenumbers[x+2]]
		
		if codon!=SNPcodon:
			
			for y,z in enumerate(codon):
				if SNPcodon[y]!=z:
					newSNPcodon=codon[:y]+SNPcodon[y]+codon[y+1:]
					if 'N' in codon or 'N' in newSNPcodon or '-' in codon or '-' in newSNPcodon:
						SNPtype[codonposn[y]]='-'
					elif geneticcode[newSNPcodon]=='*':
						SNPtype[codonposn[y]]='2'
					elif geneticcode[newSNPcodon]=='*':
						SNPtype[codonposn[y]]='3'
					elif geneticcode[newSNPcodon] == geneticcode[codon]:
						SNPtype[codonposn[y]]='S'
					else:
						SNPtype[codonposn[y]]='N'
				if not 'N' in codon and not 'N' in SNPcodon and geneticcode[codon]!=geneticcode[SNPcodon]:
					AAfromtype[codonposn[y]]=geneticcode[codon]
					AAtotype[codonposn[y]]=geneticcode[SNPcodon]
		
		if 'N' in codon or 'N' in SNPcodon or '-' in codon or '-' in SNPcodon:
			gapcount=gapcount+3
			continue
		
#		if geneticcode[codon]=='*' or geneticcode[SNPcodon]=='*':
#			continue
#			print codon, x
		
		
		#s=float(codonsynonyms[codon])/3
		#n=float(3-s)
		S1=S1+(float(codonsynonyms[codon]))
		S2=S2+(float(codonsynonyms[SNPcodon]))
		
		
		
		sd=0.0
		nd=0.0
		
		pathcount=0
		if codon!=SNPcodon:
					
			sd, nd, pathcount=countcodonchanges(codon, SNPcodon, geneticcode, sd, nd)
		
		#if nochanges==2:
		#	sd=float(sd)/2
		#	nd=float(nd)/2
		#elif nochanges==3:
		#	sd=float(sd)/6
		#	nd=float(nd)/6
		if pathcount>0:
			sd=float(sd)/pathcount
			nd=float(nd)/pathcount
		#elif codon!= SNPcodon:
		#	print codon, SNPcodon, "All paths lead to early stop codon in SNP sequence."
			#return
		
		
		Sd=Sd+sd
		Nd=Nd+nd
		#pS=pS+(sd/s)
		#pN=pN+(nd/n)
		
		#print codon, SNPcodon
	
	S=(S1+S2)/2
	N=len(CDS)-gapcount-S
	
	#pNb=pN/numcodons
	#pSb=pS/numcodons
	
	if N!=0:
		pN=Nd/N
	if S!=0:
		pS=Sd/S
	
	
	if pS==0:
		#print "No sites are synonymous."
		return [N, S, dN, dS, pN, pS, varianceS, varianceN, z, (len(CDS)-gapcount), Nd, Sd, Fisher], SNPtype, AAfromtype, AAtotype
	if pS<0.75 and pN<0.75:
		dS=(-3*(math.log(1-((pS*4)/3))))/4
		dN=(-3*(math.log(1-((pN*4)/3))))/4
		varianceS=(9 * pS * (1 -pS))/(((3 - 4 *pS) **2) * (len(CDS)-gapcount));
		varianceN=(9 * pN * (1 -pN))/(((3 - 4 *pN) **2) * (len(CDS)-gapcount));
		z=(dN - dS) / math.sqrt(varianceS + varianceN)
		#Fisher exact test of S and N
		#oddsratio, Fisher = stats.fisher_exact([[N, S], [Nd, Sd]])
		Fisher=fisher.pvalue(N,S,Nd,Sd).two_tail
		
	else:
		#print "Too divergent for JC! Using pN/pS instead."
		dS=pS
		dN=pN
		
	
	
	#print dN, dS, S, N, S+N, Sd, Nd, pS, pN#, pSb, pNb, dS, dN, pN/pS, dN/dS
	#print "N =", N
	#print "S =", S
	#print "dN/dS =", dN/dS
	
	return [N, S, dN, dS, pN, pS, varianceS, varianceN, z, (len(CDS)-gapcount), Nd, Sd, Fisher], SNPtype, AAfromtype, AAtotype
	





#########################################################
# Function to concatenate CDSs for one of the sequences #
#########################################################

def concatenate_CDS_sequences(record, sequence, ref):
	
	def add_feature_seq(feature):
		start=ref_pos_to_aln_pos[feature.location.start.position]
		end=ref_pos_to_aln_pos[feature.location.end.position-1]
		
		my_seq= Seq(sequences[sequence][start : end+1])
		if feature.strand==-1:
			my_seq=my_seq.reverse_complement()
			#print len(my_seq),
		my_seq=str(my_seq)
		#my_seq=str(feature.extract(sequences[sequence]))
			
		my_ref= Seq(sequences[ref][start : end+1])
		if feature.strand==-1:
			my_ref=my_ref.reverse_complement()
			#print len(my_ref)
		my_ref=str(my_ref)
		
		#my_ref=str(feature.extract(sequences[ref]))
		
		for base in xrange(len(my_ref)):
			if my_ref[base]=="-" and my_seq[base]=="-":
				continue
			elif my_ref[base]=="-" or my_seq[base]=="-":
				return

		if len(str(my_seq.replace("-",""))) % 3 or len(str(my_ref.replace("-",""))) % 3:
			#print float(len(str(my_seq)))/3, float(len(str(my_ref))/3)
			#print feature.location.start.position, feature.location.end.position, feature.location.end.position-feature.location.start.position
			#print feature
			#print feature.sub_features
			return
				

			
		if feature.strand==-1:
			for x, y in enumerate(range(end,start-1,-1)):
				if my_ref[x]!="-":
					CDSbasenumbers.append(y)
		else:
			for x, y in enumerate(range(start, end+1)):
				if my_ref[x]!="-":
					CDSbasenumbers.append(y)
		
		my_ref=my_ref.replace("-","")
		my_seq=my_seq.replace("-","")
		
		
		concatenated_sequence.append(my_seq)
		concatenated_ref.append(my_ref)
			
	
	def iterate_subfeatures(feature):
		if len(feature.sub_features)>0:
			if feature.strand==-1:
				for subfeature in feature.sub_features[::-1]:
					iterate_subfeatures(subfeature)
			else:
				for subfeature in feature.sub_features:
					iterate_subfeatures(subfeature)
		else:
			add_feature_seq(feature)
	
	
	concatenated_sequence=[]
	concatenated_ref=[]
	CDSbasenumbers=[]
	
	for feature in record.features:
		if feature.type=="CDS" and not "pseudo" in feature.qualifiers:
			add_feature_seq(feature)
		
	print len(''.join(concatenated_sequence)), len(''.join(concatenated_ref)), len(CDSbasenumbers)
	return ''.join(concatenated_sequence), ''.join(concatenated_ref), CDSbasenumbers




########################################
# Function to get a name for a feature #
########################################

def get_best_feature_name(feature):
	
	name_types=["gene", "primary_name", "systematic_id", "locus_tag"]
	
	for name in name_types:
		if feature.qualifiers.has_key(name):
			return feature.qualifiers[name][0]

	return ""


###########################################
# Function to get a product for a feature #
###########################################

def get_feature_product(feature):
	
	if feature.qualifiers.has_key("product"):
		return feature.qualifiers["product"][0]

	return ""


class SNPanalysis:
	def __init__(self, fastq='', directory='', mapped={}, runmaq='n', CDSseq=''):
		self.fastq=fastq
		self.directory=directory
		self.mapped=mapped
		self.runmaq=runmaq
		self.CDSseq=CDSseq
		self.N=0.0
		self.S=0.0
		self.dN=0.0
		self.dS=0.0
		self.pN=0.0
		self.pS=0.0
		self.varianceN=0.0
		self.varianceS=0.0
		self.z=0.0
		self.Nd=0.0
		self.Ns=0.0
		self.goodlen=0
		self.CDSseq=''
		self.snpsummary={}
		self.nummapped=0
		self.percentmapped=0.0
		self.intragenic=0
		self.intergenic=0
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
	runpileups='n'
	
	
	print '\nReading input alignment...',
	sys.stdout.flush()
	
	#else:
	lines=[]
	count=-1
	sequences={}
	names=set([])
	curseqlist=[]
	append=curseqlist.append
	for linea in open(options.inputfile, "rU"):
		linea=linea.strip()
		if len(linea)==0:
			continue
		if linea[0]==">":
			if count>-1:
				sequence=''.join(curseqlist)
				sequences[name]=sequence
			count=count+1
			curseqlist=[]
			append=curseqlist.append
			name=linea.split()[0][1:]
			
			lines.append(linea.split()[0][1:]+'\n')
		else:	
			append(linea)
	
	if count>-1:
		sequence=''.join(curseqlist)
		sequences[name]=sequence
		
	
	
	#for line in lines:
	#	words=line.strip().split('\n')
	#	sequences[words[0].split()[0]]=''.join(words[1:])
	snpstructs.append(SNPanalysis())	
	
	
	ref=options.ref

	reflen=len(sequences[sequences.keys()[0]])

	for sequence in sequences.keys():
		if len(sequences[sequence])!=reflen:
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

		
	

	###################
	# Open embl files #
	###################
	
	if options.embl!='':
		
		print "\nReading EMBL file(s)...",
		
		try:
			record = open_annotation(options.embl)
		except StandardError:
			print "\nCannot open annotation file"
			sys.exit()
#		
#		readok=False
#		try:
#			record = SeqIO.read(open(options.embl), "embl")
#			
#		except StandardError:
#			print "\nCannot open annotation file as embl"
#		else:
#			readok=True
#		if readok==False:
#		
#			print "Trying to open file as genbank instead"
#		
#			try:
#				record = SeqIO.read(open(options.embl), "gb")
#			except StandardError:
#				print "Cannot open annotation file as genbank"
#				
#		if readok==False:
#			emblfile=open(options.embl, "rU")
#			found_id_line=False
#			found_header=False
#			found_sequence=False
#			newembl=[]
#			
#			for line in emblfile:
#				if len(line.split())>0:
#					if line.split()[0]=="ID":
#						found_id_line=True
#						if len(line.strip().split(";"))>=7 and int(line.strip().split(";")[6].split()[0])==len(sequences[options.ref]):
#							newembl.append(line)
#						else:
#							print "Found invalid ID line. Will be replaced."
#							found_id_line=False
#					elif line.split()[0]=="FH":
#						found_header=True
#						newembl.append(line)
#					elif line.split()[0]=="SQ":
#						found_sequence=True
#						newembl.append(line)
#					else:
#						newembl.append(line)
#			
#			emblfile.close()
#			
#			if found_id_line==False or found_header==False or found_sequence==False:
#				
#				print "Last try. Trying to convert file into readable embl format"
#				
#				chars = string.ascii_letters + string.digits
#				tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))
#				stringfile = open(tmpname+".embl", "w")
#				
#				if options.gaps:
#					refseqforembl=sequences[options.ref].replace("-","")
#				else:
#					refseqforembl=sequences[options.ref]
#				print len(sequences[options.ref]), len(refseqforembl)
#				
#				if found_id_line==False:
#					print "Adding ID line"
#					print >> stringfile, embl_style_id(refseqforembl)
#				
#				if found_header==False:
#					print "Adding header line"
#					print >> stringfile, embl_style_header(refseqforembl)
#				
#				print >> stringfile, ''.join(newembl)
#				
#				if found_sequence==False:
#					print "Appending sequence"
#					print >> stringfile, embl_style_sequence(refseqforembl)
#				
#				stringfile.close()
#				
#				try:
#					record = SeqIO.read(open(tmpname+".embl"), "embl")
#					os.system("rm "+tmpname+".embl")
#				except StandardError:
#					os.system("rm "+tmpname+".embl")
#					DoError("Cannot read annotation file (is it in embl/genbank format?)")
#			else:
#				DoError("Cannot read annotation file (is it in embl/genbank format?)")

		print "Done"
		if len(record.seq)!=reflennogaps:
			DoError("embl file ("+str(len(record.seq))+") and alignment ("+str(reflennogaps)+") have different reference lengths")
			sys.exit()
		
	
	
	
	
	#Identify snps in the alignment
	
	print "\nIdentifying SNPs...",
	sys.stdout.flush()
	
	constants={"A":0, "C":0, "G":0, "T":0}
	
	for x in range(reflen):
		numbases=0
		foundbases={}
		
		if sequences[ref][x].upper()=="N":
			continue
		
		for key in sequences.keys():
			base=sequences[key][x].upper()
			if base not in foundbases.keys() and base not in ['-', "N"]:
				foundbases[base]=1
				numbases=numbases+1
			elif base not in ['-', "N"]:
				foundbases[base]+=1
		if numbases>1:
			snplocations.append(x)
			snpbases.append(foundbases)
			basecounts.append(numbases)
		elif numbases==1:
			constants[foundbases.keys()[0]]+=1
	
	
	print "Done"
	print "Found", len(snplocations), "sites with a SNP"
	print "Constant bases:"
	for base in ["A","C","G","T"]:
		print base+":", constants[base]
	print
	sys.stdout.flush()
	
	
	seqsort=sequences.keys()
	seqsort.sort()
	
	allsnps={}
	refbases={}
	
	alssnpsummary={}
	
				
	
	tempsnps={}
	tempsnpsummary={'A':{'C':0, 'G':0, 'T':0}, 'C':{'A':0, 'G':0, 'T':0}, 'G':{'C':0, 'A':0, 'T':0}, 'T':{'C':0, 'G':0, 'A':0}}
	
	olddone=-1
	

	if options.embl!='' and ref!='':
	
		def add_feature_pos_to_basenums(feature, pseudo, featurenum):
			start=ref_pos_to_aln_pos[feature.location.start.position]
			end=ref_pos_to_aln_pos[feature.location.end.position-1]
			
			if feature.type in ['CDS', 'rRNA', 'tRNA']:
				#print feature.id, feature.type, feature.location, start,end+1, len(embldata)
				for x in range(start,end+1):
					embldata[x]=featurenum
					
		def add_subfeature_positions(feature, pseudo, featurenum):
			if len(feature.sub_features)>0:
				for subfeature in feature.sub_features:
					add_subfeature_positions(subfeature, pseudo, featurenum)
				
			else:
				add_feature_pos_to_basenums(feature, pseudo, featurenum)
	
	
	
	
		print "\nCalculating dN/dS values...",# NEEDS FIXING!!!
		sys.stdout.flush()
		#dndsout=open('dNdS_by_CDS.out','w')
		
		sys.stdout.flush()
		dnbydsstats={}
		
		comp={'A':'T','T':'A','G':'C','C':'G', 'N':'N', '-':'-'}
		
		embldata=[]

		
		for i in range(0,reflen):
			embldata.append(-1)
		
		
		
		for y, feature in enumerate(record.features):
			if feature.type in ['CDS', 'rRNA', 'tRNA']:
				if "pseudo" in feature.qualifiers:
					pseudo=True
				else:
					pseudo=False
				add_subfeature_positions(feature, pseudo, y)
					
		tmpout=open("tmpold.txt","w")		
		print >> tmpout, embldata
		tmpout.close()
		#refCDSseq=concatenate_CDS_sequences(record, ref)
		
		
		
		print 'Done'
		sys.stdout.flush()
		tmpseqs={}
		snptypes={}
		AAfromtypes={}
		AAtotypes={}
		for sequence in seqsort:
			if sequence==ref:
				continue
			
#			h = hpy()
#			print h.heap()

			print sequence+'...',
			sys.stdout.flush()
			
			tmpCDSseq, refCDSseq, CDSbasenumbers=concatenate_CDS_sequences(record, sequence, ref)
			#print len(tmpCDSseq), len(refCDSseq), len(CDSbasenumbers)
			
			dnbydsstats[sequence], snptypes[sequence], AAfromtypes[sequence], AAtotypes[sequence]=dnbyds(refCDSseq, tmpCDSseq, CDSbasenumbers, genetic_code_number=options.genetic_code_number)
			#print dnbydsstats[sequence]
			#N, S, dN, dS, pN, pS, varianceS, varianceN, z, (len(CDS)-gapcount), Nd, Sd
			
			if dnbydsstats[sequence][3]!=0:
				print "dNdS = %.2f,  p-value = %.2f" % (dnbydsstats[sequence][2]/dnbydsstats[sequence][3], dnbydsstats[sequence][12])
				#print >> dndsout, '\t'+str(dnbydsstats[sequence][2]/dnbydsstats[sequence][3]),
				sys.stdout.flush()
			else:
				print "-"			
		
		
		#dndsout.close()









	
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
	if options.embl!='':
			print >> output, '\tCDS/rRNA/tRNA/Intergenic\tstrand\tCDS_name\tproduct\tSynonymous/Non-synonymous',
	print >> output, '\tRef_base\tSNP_base\tTotal',
	
	for name in seqsort:
		if name!=ref:
			print >> output, '\t'+name,
	print >> output, '\n',
	
#
	if options.tabfile==True:
			tabout=open(options.outfile+'_snps.tab','w')
		#covgraph=open(outfile+'_cov.txt','w')
		#covcountgraph=open(outfile+'_covcount.txt','w')
	
	#minquality=10
	#minqualityperdepth=4
	
	snpsequence={}
	for name in seqsort:
		snpsequence[name]=''
	
	
	refsequence=''
	if ref=='':
		ref=sequences.keys()[0]
	
	if options.tabfile==True:
		print >> tabout, 'ID   SNP'
	
	poscount=1
	
	#for key in sequences.keys():
		
		
	
		#snpsort=allsnps[key].keys()
		#snpsort.sort()
	
	strainsnpsummary={}
	for name in seqsort:
		if name!=ref:
			strainsnpsummary[name]={}
			strainsnpsummary[name]['count']=0
			for x in ['A', 'C', 'G', 'T']:
				strainsnpsummary[name][x]={}
				for y in ['A', 'C', 'G', 'T']:
					strainsnpsummary[name][x][y]=0
	

	
	
	for x, j in enumerate(snplocations):
		
		outstring=''
		tabstring=''
		
		#if len(allsnps.keys())>1:
		#	outstring=outstring+key+'\t'+str(j)
		#else:
		outstring=outstring+str(j+1)
		snpcolour='1'
		
		
		
		if sequences[ref][j]!='-':
			outstring=outstring+'\t'+str(aln_pos_to_ref_pos[j]+1)
			if options.embl!='':
				
#				if embldata[aln_pos_to_ref_pos[j]]<0:
#					outstring=outstring+'\tIntergenic\t-\t-\t'
#				else:
#					outstring=outstring+'\t'+record.features[embldata[aln_pos_to_ref_pos[j]]].type+'\t'+str(record.features[embldata[aln_pos_to_ref_pos[j]]].strand)+'\t'+get_best_feature_name(record.features[embldata[aln_pos_to_ref_pos[j]]])+'\t'+get_feature_product(record.features[embldata[aln_pos_to_ref_pos[j]]])
#				snptype='-'
				if embldata[j]<0:
					outstring=outstring+'\tIntergenic\t-\t-\t-'
				else:
					outstring=outstring+'\t'+record.features[embldata[j]].type+'\t'+str(record.features[embldata[j]].strand)+'\t'+get_best_feature_name(record.features[embldata[j]])+'\t'+get_feature_product(record.features[embldata[j]])
				snptype='-'
				for name in seqsort:
					if name==ref:
						continue
				#	if snptypes[name].has_key(aln_pos_to_ref_pos[j]):
					#	if snptypes[name][aln_pos_to_ref_pos[j]]=='S':
					#		snpcolour='3'
					#	elif snptypes[name][aln_pos_to_ref_pos[j]]=='N':
					#		snpcolour='2'
					#	elif snptypes[name][aln_pos_to_ref_pos[j]]=='2':
					#		snpcolour='4'
					#	elif snptypes[name][aln_pos_to_ref_pos[j]]=='3':
					#		snpcolour='4'
					#	if snptype=='-' and snptypes[name][aln_pos_to_ref_pos[j]]!='-':
					#		snptype=snptypes[name][aln_pos_to_ref_pos[j]]
					#	elif snptypes[name][aln_pos_to_ref_pos[j]] not in snptype and snptypes[name][aln_pos_to_ref_pos[j]]!='-':
					#		snptype=snptype+'/'+snptypes[name][aln_pos_to_ref_pos[j]]
					#	if record.features[embldata[aln_pos_to_ref_pos[j]]].qualifiers.has_key("pseudo"):
					#		snpcolour='11'
					#		snptypes[name][aln_pos_to_ref_pos[j]]="P"
					if snptypes[name].has_key(j):
						if snptypes[name][j]=='S':
							snpcolour='3'
						elif snptypes[name][j]=='N':
							snpcolour='2'
						elif snptypes[name][j]=='2':
							snpcolour='4'
						elif snptypes[name][j]=='3':
							snpcolour='4'
						if snptype=='-' and snptypes[name][j]!='-':
							snptype=snptypes[name][j]
							if snptype=='2':
								snptype='SNOP'
							elif snptype=='3':
								snptype='STIP'
						elif snptypes[name][j] not in snptype and snptypes[name][j]!='-':
							snptmp=snptypes[name][j]
							if snptmp=='2':
								snptype='SNOP'
							elif snptmp=='3':
								snptype='STIP'
							snptype=snptype+'/'+snptmp
						if record.features[embldata[j]].qualifiers.has_key("pseudo"):
							snpcolour='11'
							snptypes[name][j]="P"
						
									
				outstring=outstring+'\t'+snptype
		else:
			outstring=outstring+'\t-'
			snpcolour='1'
			if options.embl!='':
				outstring=outstring+'\t-\t-\t-\t-\t-'
		
		
		outstring=outstring+'\t'+sequences[ref][j]
		refbase=sequences[ref][j]
		refsequence=refsequence+refbase
		
			
		snpbase=''
		snpbasecount=0
		snpbaselist=[]
		for base in snpbases[x].keys():
			if base!=refbase and len(snpbaselist)==0 and not base in ["N", "-"]:
				snpbasecount=str(snpbases[x][base])
				snpbase=base
				snpbaselist.append(base)
			elif base!=refbase and not base in snpbaselist and not base in ["N", "-"]:
				snpbase=snpbase+','+base
				snpbasecount=snpbasecount+','+str(snpbases[x][base])
				snpbaselist.append(base)
	
		outstring=outstring+'\t'+snpbase+'\t'+snpbasecount# need to think how to do this
		
		if options.tabfile==True and sequences[ref][j]!='-':
			tabstring=tabstring+'FT   SNP             '+str(aln_pos_to_ref_pos[j]+1)+'\n'
			tabstring=tabstring+'FT                   /note="refAllele: '+sequences[ref][j]
			#tabstring=tabstring+'FT				   /SNPAllele="'+snpbase+'"\n'
			tabstring=tabstring+' SNPstrains: '
			

		

		numsnpbases=1
		sitesnpbases=[]
		taxastring='\nFT                   /taxa="'
		for name in seqsort:
			if name!=ref:
				if sequences[name][j]!='-' and sequences[name][j]!='N':# and name.SNPs[key][j][1]>8:# and name.SNPs[key][j][2]>(name.SNPs[key][j][3]*4):
					if sequences[name][j]!=sequences[ref][j]:
						if sequences[ref][j]!='-':
							strainsnpsummary[name][sequences[ref][j]][sequences[name][j]]+=1
						if options.tabfile==True and sequences[ref][j]!='-':
							tabstring=tabstring+name+'='+sequences[name][j]+' '
							taxastring=taxastring+name+" "
							if options.embl!='':
								if snptypes[name].has_key(j):
									if snptypes[name][j]=='S':
										tabstring=tabstring+'(synonymous) '
									elif snptypes[name][j]=='N':
										tabstring=tabstring+'(non-synonymous) '
									elif snptypes[name][j]=='-':
										tabstring=tabstring+'(gap in SNP codon) '
									elif snptypes[name][j]=='2':
										tabstring=tabstring+'(SNP codon is STOP) '
									elif snptypes[name][j]=='3':
										tabstring=tabstring+'(ref codon is STOP) '
									if name in AAfromtypes and name in AAtotypes:
										if j in AAfromtypes[name] and j in AAtotypes[name]:
											tabstring=tabstring+'(AA '+AAfromtypes[name][j]+'->'+AAtotypes[name][j]+') '
								
						outstring=outstring+'\t'+sequences[name][j]#+str(name.mapped[key][j][1])
						
						#if sequences[name][j] not in sitesnpbases:
						#	sitesnpbases.append(sequences[name][j])
						#	numsnpbases=numsnpbases+1
						#name.SNPs[j]=''
						
					else:
						#outstring=outstring+'\t'+sequences[name][j]
						outstring=outstring+'\t.'
					
				else:
					#snpsequence[name]=snpsequence[name]+'-'
					outstring=outstring+'\t'+sequences[name][j]
				
			snpsequence[name]=snpsequence[name]+sequences[name][j]
		taxastring=taxastring+'"'
		print >> output, outstring
		if options.tabfile==True and sequences[ref][j]!='-':
			print >> tabout, tabstring+'"\nFT                   /colour='+snpcolour+taxastring

	
	
	
	#Creat SNP alignment file
	
#	
#	if options.align==True:
#		
#		count=0
#		
#		#print >> alnfile, len(snpsequence.keys()), len(refsequence)
#		
#		phylipnames={}
#		
#		for name in snpsequence.keys():
#			#print >> alnfile, name.replace('pool','').replace('.snps','')+' '*(10-len(name))+snpsequence[name]
#			#print >> alnfile, name.replace('pool','').replace('.snps','')+' '+snpsequence[name]
#			phylipname=name.replace('pool','').replace('.snps','')
#			if len(name)>9:
#				phylipname=phylipname[:9]
#			count=1
#			while phylipname in phylipnames.keys():
#				phylipname=phylipname[:(len(str(count))+1)*-1]+"_"+str(count)
#			
#			phylipnames[name]=phylipname
#		
#		#print phylipnames
#		
#		alignment = Generic.Alignment(Gapped(IUPAC.unambiguous_dna, "-"))
#		for name in snpsequence.keys():
#			alignment.add_sequence(phylipnames[name], snpsequence[name])
#		
#		
#		AlignIO.write([alignment], open(options.outfile+'_snps.aln', 'w'), "phylip")
		
						
		#alnfile.close()
	
	if options.align==True:
		
		alignment = Generic.Alignment(Gapped(IUPAC.unambiguous_dna, "-"))
		for name in seqsort:
			 
			if len(snpsequence[name].replace("-","").replace("N",""))>float(len(snpsequence[name]))*(float(options.exclude)/100):
				alignment.add_sequence(name, snpsequence[name])
			else:
				print name, "excluded from snp alignment as it is < "+str(options.exclude)+"% mapped"
		
		
		AlignIO.write([alignment], open(options.outfile+'_snps.aln', 'w'), "fasta")
		
						
		#alnfile.close()
	
	
	
	
	
	#Create tab file for each sequence vs ref (extended to length of total alignment to allow insertions)
	#print ref_pos_to_aln_pos
	if options.tabfile==True and reflennogaps<reflen:
		extrefout=open('Padded_'+ref+'.fna', 'w')
		print >> extrefout, '>'+ref
		count=0
		outstring=''
		for x in range(len(sequences[ref])):
			outstring=outstring+sequences[ref][x]
			count=count+1
			if count==60:
				outstring=outstring+'\n'
				count=0
		print >> extrefout, outstring
		outstring=''
		#print reflen, len(sequences[ref]), len(sequences[ref].replace('-',''))
		extrefout.close()
		if options.embl!='':
			lines=open(options.embl,'rU').readlines()
			extrefout=open('Padded_'+ref+'.embl', 'w')
			for line in lines:
				if len(line.strip().split())==3 and line.strip().split()[1]=='source':
					#print '..'+str(len(sequences[ref])), '..'+str(len(sequences[ref].replace('-','')))
					print >> extrefout, line.replace( '..'+str(len(sequences[ref].replace('-',''))), '..'+str(len(sequences[ref]))),
				elif len(line.strip().split())==3 and len(line.strip().split()[2].split('..'))==2:
					positions=line.strip().split()[2]
					#print positions
					num_list = re.findall(r'.[0-9]+.', positions)
					#print num_list
					for numa in num_list:
					
						num = re.findall(r'[0-9]+', numa)[0]
						
						#print numa, num,
						
						numb=numa.replace(num,str(ref_pos_to_aln_pos[int(num)-1]+1))
						#print numb
						
						positions=positions.replace(numa, numb)
						 
						#print num, num[0]+str(ref_pos_to_aln_pos[int(num[1:-1])-1]+1)+num[-1], ref_pos_to_aln_pos[int(num[1:-1])-1]
					
					
					#print positions
					#start=ref_pos_to_aln_pos[int(line.strip().split()[2].split('..')[0])-1]
					#end=ref_pos_to_aln_pos[int(line.strip().split()[2].split('..')[1])-1]
					#posn=line[:22]+str(start+1)+'..'+str(end+1)
					print >> extrefout, line[:21]+positions
				elif len(line.strip().split())>1 and line.strip().split()[0]=='SQ':
					A=0
					C=0
					G=0
					T=0
					Others=0
					
					for x in sequences[ref]:
						if x=='A':
							A=A+1
						elif x=='C':
							C=C+1
						elif x=='G':
							G=G+1
						elif x=='T':
							T=T+1
						else:
							Others=Others+1
					print >> extrefout, "SQ   Sequence "+str(len(sequences[ref]))+" BP; "+str(A)+" A; "+str(C)+" C; "+str(G)+" G; "+str(T)+" T; "+str(Others)+" other;"
					
					count=0
					countb=0
					totalsofar=0
					outstring=''
					for x in range(len(sequences[ref])):
						if count==0:
							outstring=outstring+'	 '
						outstring=outstring+sequences[ref][x].lower()
						count=count+1
						totalsofar=totalsofar+1
						countb=countb+1
						
						if count==60:
							outstring=outstring+'	   '+str(totalsofar)+'\n'
							count=0
							countb=0
						if countb==10:
							outstring=outstring+' '
							countb=0
					
					if count!=0:
						while count!=0:
							outstring=outstring+' '
							count=count+1
							countb=countb+1
							if count==60:
								outstring=outstring+'	   '+str(totalsofar)+'\n//'
								count=0
								countb=0
							if countb==10:
								outstring=outstring+' '
								countb=0
	
					
					print >> extrefout, outstring,
					outstring=''
					
					break
					
					
				else:
					print >> extrefout, line.strip()
			extrefout.close()
#		for name in sequences.keys():
#			indtabfile=open(name+'.tab', 'w')
#			print >> indtabfile, 'ID   SNP'
#			delstart=0
#			insstart=0
#			for i,j in enumerate(sequences[name]):
#				if j=='-':
#					if i==0 or sequences[name][i-1]!='-':
#						delstart=i
#					if (i+1)==len(sequences[ref]) or sequences[name][i+1]!='-':
#						if delstart!=i:
#							print >> indtabfile, 'FT   misc_feature	'+str(delstart+1)+'..'+str(i+1)
#						else:
#							print >> indtabfile, 'FT   misc_feature	'+str(i+1)
#						print >> indtabfile, 'FT				   /colour=3'#Deletion = grey
#						
##				#elif sequences[ref][i]=='-' and j!=sequences[ref][i]:
##				#	if i==0 or sequences[ref][i-1]!='-' or sequences[sequence][i-1]=='-':
##				#		insstart=i
#					#if (i+1)==len(sequences[ref]) or sequences[ref][i+1]!='-' or sequences[sequence][i+1]=='-':
#						#if insstart!=i:
#						#	print >> indtabfile, 'FT   misc_feature	'+str(insstart+1)+'..'+str(i+1)
#						#else:
#						#	print >> indtabfile, 'FT   misc_feature	'+str(i+1)
#						#print >> indtabfile, 'FT				   /colour=10'#Insertion = orange
#						
#				elif j!=sequences[ref][i] and sequences[ref][i]!='-':
#					print >> indtabfile, 'FT   misc_feature	'+str(i+1)
#					print >> indtabfile, 'FT                   /refBase="'+sequences[ref][i]+'"'
#					#tabstring=tabstring+'FT                   /SNPAllele="'+snpbase+'"\n'
#					print >> indtabfile, 'FT                   /SNPBase="'+j+'"'
#					if embl!='' and snptypes[name].has_key(aln_pos_to_ref_pos[i]):
#						if snptypes[name][aln_pos_to_ref_pos[i]]=='S':
#							print >> indtabfile, 'FT				   /colour=2'#Synonymous=red
#						elif snptypes[name][aln_pos_to_ref_pos[i]]=='N':
#							print >> indtabfile, 'FT				   /colour=4'#Non-synonymous=blue
#						elif snptypes[name][aln_pos_to_ref_pos[i]]=='1':
#							print >> indtabfile, 'FT                   /colour=1'#add synonymous stuff here
#						else:
#							print >> indtabfile, 'FT                   /colour=1'#add synonymous stuff here
#					elif embl!='':
#						print >> indtabfile, 'FT                   /colour=1'#add synonymous stuff here
#					else:
#						print >> indtabfile, 'FT                   /colour=14'#add synonymous stuff here
#
#			indtabfile.close()
#	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	#[N, S, dN, dS, pN, pS, varianceS, varianceN, z, (len(CDS)-gapcount), Nd, Sd]
	
	#sys.exit()
	summary='y'
	if summary=='y':
		summaryfile=open(options.outfile+'_summary.out','w')
		
		if options.embl!="":
			headingline='Strain\tPercent Mapped (100-percent Ns)\tdN/dS\t2-tailed Fisher\'s exact test p-value'
		else:
			headingline='Strain\tPercent Mapped'
		
		for x in ['A', 'C', 'G', 'T']:
			for y in ['A', 'C', 'G', 'T']:
				if x!=y:
					headingline=headingline+'\t'+x+'->'+y
		
		headingline=headingline+'\tTotal SNPs'
		
		print >> summaryfile, headingline
		
		for name in seqsort:
			if name==ref:
				continue
			intragenic=0
			intragenic_mapped=0
			
			if options.embl!="" and dnbydsstats[name][3]!=0:
				strainline=name+'\t'+str(100*(float(len(sequences[name].replace("N","")))/reflennogaps))+'\t'+str(dnbydsstats[name][2]/dnbydsstats[name][3])+'\t'+str(dnbydsstats[name][12])
			elif options.embl!="":
				strainline=name+'\t'+str(100*(float(len(sequences[name].replace("N","").replace("-","")))/len(sequences[name].replace("-",""))))+'\t-'
			else:
				strainline=name+'\t'+str(100*(float(len(sequences[name].replace("N","").replace("-","")))/len(sequences[name].replace("-",""))))
			count=0
			for x in ['A', 'C', 'G', 'T']:
				for y in ['A', 'C', 'G', 'T']:
					if x!=y:
						strainline=strainline+'\t'+str(strainsnpsummary[name][x][y])
						count=count+strainsnpsummary[name][x][y]

			
			
			strainline=strainline+"\t"+str(count)
			print >> summaryfile, strainline
		

		print >> summaryfile, '\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t'+str(len(snplocations))
		summaryfile.close()
		
		print 'Done'
	
	if options.tabfile==True:
		tabout.close()
	

	output.close()
	
	#sys.exit()
	if options.raxml==True:
		outlen=0
		userinput='x'
		if '/' in options.outfile:
			outlen=-1*len(options.outfile.split('/')[-1])
		if os.path.isfile('RAxML_info.'+options.outfile.split('/')[-1]) and options.overwrite==False:
			print '\nRAxML files with extension '+options.outfile.split('/')[-1]+' already exist!'
			while userinput not in ['y','n']:
				userinput=raw_input('Overwrite? (y/n): ')
				userinput.lower()
			if userinput=='y':
				print 'RAxML files with extension '+options.outfile.split('/')[-1]+' will be overwritten'
				os.system('rm RAxML_*'+options.outfile.split('/')[-1])
				
		RAxMLcommand="/nfs/pathogen/sh16_scripts/run_RAxML.py -a "+options.outfile+"_snps.aln -t DNA -v "+options.asrv+" -w -o "+options.outfile.split('/')[-1]
		
		if options.pinvar:
			RAxMLcommand=RAxMLcommand+" -i"
		if options.f:
			RAxMLcommand=RAxMLcommand+" -F"
		if options.bootstrap>0:
			RAxMLcommand=RAxMLcommand+" -b "+str(options.bootstrap)
			if options.fast:
				RAxMLcommand=RAxMLcommand+" -f"
		else:
			RAxMLcommand=RAxMLcommand+" -b 0"
		print RAxMLcommand
		os.system(RAxMLcommand)
		
#		if options.bootstrap>0:
#			os.system("/nfs/pathogen/sh16_scripts/run_RAxML.py -a "+options.outfile+"_snps.aln -t DNA -v "+options.asrv+" -l -w -o "+options.outfile.split('/')[-1])
#			#print "Running RAxML phylogeny with "+options.model+" model of evolution and "+str(options.bootstrap)+" bootstrap replicates..."
#			#os.system(RAXML_PATH+"RAxML -f a -x "+str(randrange(1,99999))+" -p "+str(randrange(1,99999))+" -# "+str(options.bootstrap)+" -m "+options.model+" -s "+options.outfile+"_snps.aln -n "+options.outfile.split('/')[-1])
#		else:
#			os.system("/nfs/pathogen/sh16_scripts/run_RAxML.py -a "+options.outfile+"_snps.aln -t DNA -v "+options.asrv+" -b 0 -w -o "+options.outfile.split('/')[-1])
#			#print "Running RAxML phylogeny with "+options.model+" model of evolution"
#			#os.system(RAXML_PATH+"RAxML -f d -m "+options.model+" -s "+options.outfile+"_snps.aln -n "+options.outfile.split('/')[-1])
#		
#	
	
	
	
	print "All finished.\n"
