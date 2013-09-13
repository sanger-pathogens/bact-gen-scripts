#!/usr/bin/env python

SMALT_DIR=""
SAMTOOLS_DIR=""
BCFTOOLS_DIR=""

##################
# Import modules #
##################

import os, sys, string
from random import randint, choice
from optparse import OptionParser
import pysam
from numpy import min, max, median, mean, std
from scipy.stats import mannwhitneyu, ttest_ind
from math import sqrt, pow

##############################################
## Function to reverse complement a sequence #
##############################################

def revcomp(sequence):
	rev=sequence[::-1]
	revcomp=''
	d={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'g':'c', 'c':'g', "n":"n", "N":"N"}
	for i in rev:
		if d.has_key(i):
			revcomp=revcomp+d[i]
		else:
			revcomp=revcomp+i
	
	return revcomp



##############################
# Get command line arguments #
##############################

def get_user_options():
	usage = "usage: %prog [options]"
	version="%prog 1.0. Written by Simon Harris, Wellcome Trust Sanger Institute, 2012"
	parser = OptionParser(usage=usage, version=version)
	
	#do not allow arguments to be interspersed. e.g. -a -b arg1 agr2. MUST be -a arg1 -b arg2.
	parser.disable_interspersed_args()
	
	parser.add_option("-C", "--crunch", action="store", dest="crunch", help="Crunch file", default="", metavar="FILE")
	parser.add_option("-c", "--contigs", action="store", dest="contigs", help="File containing contigs", default="", metavar="FILE")
	parser.add_option("-s", "--scaffolds", action="store", dest="scaffols", help="File containing scaffolds", default="", metavar="FILE")
	

	return parser.parse_args()

(options, args)=get_user_options()

if not os.path.isfile(options.contigs):
	print "Could not find contigs file"
	sys.exit()
if not os.path.isfile(options.crunch):
	print "Could not find crunch file"
	sys.exit()
try:
	contigsfile=open(options.contigs, "rU").read()
except StandardError:
	print "Could not open contigs file"
	sys.exit()
try:
	crunchfile=open(options.crunch, "rU")
except StandardError:
	print "Could not open crunch file"
	sys.exit()
try:
	scaffoldsfile=open(options.scaffolds, "rU").read()
except StandardError:
	print "Could not open scaffolds file"
	sys.exit()


contigs={}
contiglist=[]
contigseqs={}
total=0
for line in contigsfile.split(">")[1:]:
	bits=line.split("\n")
	contigs[bits[0].split()[0]]=len(''.join(bits[1:]).upper())
	contigseqs[bits[0].split()[0]]=''.join(bits[1:]).upper()
	contiglist.append([bits[0].split()[0],total,total+len(''.join(bits[1:]).upper())])
	total+=len(''.join(bits[1:]).upper())

scaffolds={}
scaffoldlist=[]
scaffoldseqs={}
total=0
for line in scaffoldsfile.split(">")[1:]:
	bits=line.split("\n")
	scaffolds[bits[0].split()[0]]=len(''.join(bits[1:]).upper())
	scaffoldseqs[bits[0].split()[0]]=''.join(bits[1:]).upper()
	scaffoldlist.append([bits[0].split()[0],total,total+len(''.join(bits[1:]).upper())])
	total+=len(''.join(bits[1:]).upper())

included=[]
excluded=[]
for line in crunchfile:
	words=line.split()
	print words[1], (int(words[3])+1)-int(words[2]), contigs[words[4]], ((int(words[3])+1)-int(words[2]))==contigs[words[4]]
	if ((int(words[3])+1)-int(words[2]))==contigs[words[4]]:
		included.append(words[4])
output=open("out.tab", "w")
for contig in contiglist:
	if contig[0] in included:
		print contig[0], "included"
	else:
		print contig[0], "excluded"
		print >> output, "FT   misc_feature    "+str(contig[1]+1)+".."+str(contig[2])
	
output.close()



	
	
	
