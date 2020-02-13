#!/usr/bin/env python
#/usr/bin/env python

import os, sys

##################
# Import modules #
##################

from Bio import AlignIO
from Bio.Align.Generic import Alignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, SeqFeature
from optparse import OptionParser

from modules.Si_nexus import get_ref_to_alignment_translations
from modules.Si_general import DoError
from modules.Si_SeqIO import read_alignment, open_annotation
from modules.Si_SNPs_temp import revcomp



##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options] fasta/multifasta input files"
	parser = OptionParser(usage=usage)

	parser.add_option("-a", "--alignment", action="store", dest="alignment", help="Input alignment file", default="", metavar="FILE")
	parser.add_option("-c", "--CDSlist", action="store", dest="CDSlist", help="File containing list of CDSs to extract", default="", metavar="FILE")
	parser.add_option("-e", "--embl", action="store", dest="embl", help="Emble annotation file", default="", metavar="FILE")
	parser.add_option("-r", "--reference", action="store", dest="reference", help="Name of referene in alignment file", default="")
	parser.add_option("-o", "--prefix", action="store", dest="prefix", help="output file prefix", default="")

	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):

	if options.alignment=='':
		DoError('No alignment file selected')
	elif not os.path.isfile(options.alignment):
		DoError('Cannot find file '+options.alignment)
	
	if options.CDSlist=='':
		DoError('No CDS list file selected')
	elif not os.path.isfile(options.CDSlist):
		DoError('Cannot find file '+options.CDSlist)
	
	if options.alignment=='':
		DoError('No embl file selected')
	elif not os.path.isfile(options.embl):
		DoError('Cannot find file '+options.embl)
	
	if options.reference=='':
		DoError('No reference name selected')
	
	if options.prefix=='':
		DoError('No output file prefix selected')

		
		
	return




if __name__ == "__main__":
	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)

	
	myset=set()
	seqorder=[]
	for line in open(options.CDSlist,"rU"):
		myset.add(line.strip())
		seqorder.append(line.strip())
	
	#Read the alignment file
		
	try:
		alignment=read_alignment(options.alignment)
	except StandardError:
		DoError("Cannot open alignment file")
	
	concatseqs={}
	foundref=False
	for sequence in alignment:
		concatseqs[sequence.name]=[]
		if sequence.name==options.reference:
			foundref=True
	
	if not foundref:
		DoError("Can't find reference in alignment")
	
	ref_to_alignment, alignment_to_ref=get_ref_to_alignment_translations(options.reference, alignment)
	
	
	try:
		emblrecord=open_annotation(options.embl)
	except StandardError:
		DoError("Cannot open annotation file "+options.embl+" please check the format")
	
	print seqorder
	
	count=0
	for feature in emblrecord.features:
		if feature.type=="CDS":
			
			foundname=False
			for nametype in ["gene", "locus_tag", "systematic_id"]:
				if nametype in feature.qualifiers and feature.qualifiers[nametype][0] in myset:
					foundname=True
					break
			if not foundname:
				continue
		
			#print feature
			count+=1
			start=feature.location.nofuzzy_start
			end=feature.location.nofuzzy_end
			
			alignment_start=ref_to_alignment[start]
			alignment_end=ref_to_alignment[end]
			
			print "Extracting", feature.qualifiers[nametype][0]
			
			outfile=open(options.prefix+"_"+feature.qualifiers[nametype][0]+".fasta", "w")
			for sequence in alignment:
				if feature.strand==1:
					print >> outfile, ">"+sequence.name
					print >> outfile, str(sequence.seq[alignment_start:alignment_end])
					concatseqs[sequence.name].append(str(sequence.seq[alignment_start:alignment_end]))
				elif feature.strand==-1:
					print >> outfile, ">"+sequence.name
					print >> outfile, revcomp(str(sequence.seq[alignment_start:alignment_end]))
					concatseqs[sequence.name].append(revcomp(str(sequence.seq[alignment_start:alignment_end])))
			outfile.close()
					
			
	outfile=open(options.prefix+"_MLST.mfa", "w")
	for sequence in concatseqs:
		
		print >> outfile, ">"+sequence
		print >> outfile, '\n'.join(concatseqs[sequence])
	outfile.close()
			
			
			
