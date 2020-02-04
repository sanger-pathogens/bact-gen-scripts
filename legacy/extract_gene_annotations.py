#!/usr/bin/env python


##################
# Import modules #
##################

import string, re, os, sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, Gapped
from Bio.SeqRecord import SeqRecord
sys.path.extend(map(os.path.abspath, ['/nfs/pathogen/sh16_scripts/modules/']))
from Si_general import *
from Si_SeqIO import *
from optparse import OptionParser


################################
# Get the command line options #
################################


def main():
	usage = "usage: %prog [options] args"
	parser = OptionParser(usage=usage)

	parser.add_option("-e", "--embl", action="store", dest="embl", help="embl file", type="string", metavar="FILE", default="")
	parser.add_option("-l", "--list", action="store", dest="list", help="file containing list of genes", type="string", metavar="FILE", default="")
	
	return parser.parse_args()
	
################
# Main program #
################		

if __name__ == "__main__":

	(options, args) = main()
	
	if options.embl=="":
		print "No embl file specified"
		sys.exit()
	elif not os.path.isfile(options.embl):
		print "Cannot find embl file:", options.embl
		sys.exit()
	
	if options.list=="":
		print "No list file specified"
		sys.exit()
	elif not os.path.isfile(options.list):
		print "Cannot find embl file:", options.list
		sys.exit()
	
	try:
		emblrecord=open_annotation(options.embl, quiet=True)
	except StandardError:
		print "Cannot open embl file:", options.embl
		sys.exit()
	
	try:
		genelist=open(options.list, "rU").readlines()
	except StandardError:
		print "Cannot open list file:", options.list
		sys.exit()
	
	searchgenes=set([])
	
	for gene in genelist:
		searchgenes.add(gene.strip().split()[0])
	
	pseudodict={}
	print "\t".join(["Name", "Product", "Length"])
	for feature in emblrecord.features:
		if feature.type=="CDS":
			featureseq=feature.extract(emblrecord.seq)
			length=len(featureseq)
			
			name=""
			for nametype in ["systematic_id", "gene", "locus_tag"]:
				if nametype in feature.qualifiers:
					name=feature.qualifiers[nametype][0]
					break
			if "product" in feature.qualifiers:
				product=feature.qualifiers["product"][0]
			else:
				product=""
			
			if name in searchgenes:
			
				print "\t".join([name, product, str(length)])
	
