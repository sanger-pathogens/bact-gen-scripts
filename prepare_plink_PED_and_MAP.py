#!/usr/bin/env python
import string, re, gzip
import os, sys
from Bio import SeqIO
sys.path.extend(map(os.path.abspath, ['/nfs/pathogen/sh16_scripts/modules/']))
from Si_SeqIO import *
from optparse import OptionParser


def DoError(message):
	print "Error:", message
	sys.exit()

##########################################
# Function to Get command line arguments #
##########################################


def main():

        usage = "usage: %prog [options]"
        parser = OptionParser(usage=usage)

        parser.add_option("-a", "--alignment", action="store", dest="alignment", help="alignment file name", default="")
        parser.add_option("-m", "--metadata", action="store", dest="metadata", help="metadata csv file name", default="")
        parser.add_option("-c", "--column", action="store", dest="column", help="metadata column number containing phenotype", type="int", default=2)
        parser.add_option("-s", "--separator", action="store", dest="separator", help="if column separator character in metadata file is not tab, specify it here",default="\t")
	parser.add_option("-o", "--output", action="store", dest="output", help="output file name prefix", default="")
	
        return parser.parse_args()

################
# Main program #
################                

if __name__ == "__main__":

        #Get command line arguments
        
        (options, args) = main()

        #Do some checking of the input files
        if options.metadata=="":
                DoError("No metadata file specified")
		
	if not os.path.isfile(options.metadata):
                DoError('Cannot find file '+options.metadata)
		
	
	metadata={}
	
	for line in open(options.metadata, "rU"):
		words=line.strip().split(options.separator)
		
		if len(words)<options.column or words[0]=="":
			continue
		
		if words[0] in metadata:
			print "WARNING: found duplicate entry in metadata for", words[0]

		else:
			metadata[words[0]]=words[options.column-1]
	
	
	if options.output=="":
                DoError("No output file specified")
		
        if options.alignment=="":
                DoError("No alignment file specified")
		
	if not os.path.isfile(options.alignment):
                DoError('Cannot find file '+options.alignment)
        
        try:
                alignment=read_alignment(options.alignment)
                #alignment = AlignIO.read(open(options.alignment), "fasta")
        except StandardError:
                DoError("Cannot open alignment file "+options.alignment+". Is it in the correct format?")
	
	
	
	biallelic_bases=[]
	excluded=[]
	for basenum in xrange(0,len(str(alignment[0].seq))):
		alleles=["?", "N", "-"]
		added_alleles=0
		bases=alignment[:,basenum]
		for base in bases:
			if not base in alleles:
				added_alleles+=1
				alleles.append(base)
		if added_alleles==2:
			biallelic_bases.append(basenum)
		elif added_alleles>2:
			excluded.append(basenum+1)
		
	print len(excluded), "bases excluded as not biallelic"
	
	
	output=open(options.output+".map", "w")
	#map output y	count	0	base
	
	for i, basenum in enumerate(biallelic_bases):
		print >> output, '\t'.join(["y", str(i+1), "0", str(basenum+1)])
		
	
	output.close()
	
	
	#output=open("tmp_bases", "w")
	#for base in excluded:
	#	print >> output, base
	#output.close()
	
		
	output=open(options.output+".ped", "w")
		
	for record in alignment:
		if not record.id in metadata:
			print "Sample in alignment not in metadata:", record.id
		else:
			outseq=[]
			sequence=str(record.seq)
			for basenum in biallelic_bases:
				base=sequence[basenum]
				if base in ["-", "N", "?"]:
					base="0"
				outseq.append(base)
				outseq.append(base)
					
			print >> output, record.id,record.id, 0, 0, 1, metadata[record.id], ' '.join(outseq)
	
	output.close()
