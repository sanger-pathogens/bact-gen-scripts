#!/usr/bin/env python


##################
# Import modules #
##################

import string, re, os, sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, Gapped
from Bio.SeqRecord import SeqRecord
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
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
	parser.add_option("-o", "--output", action="store", dest="outfile", help="name for output file", type="string", metavar="FILE", default="protein.fasta")
	parser.add_option("-s", "--include_stop", action="store_true", dest="inc_stop", help="include stop codons in protein sequences (if this option is not chosen the translation will stop at the first stop codon found, so pesudogenes will be truncated)", metavar="FILE", default=False)
	
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
	
	
	try:
		emblrecord=open_annotation(options.embl, quiet=True)
	except StandardError:
		print "Cannot open embl file:", options.embl
		sys.exit()
	
	output=open(options.outfile,"w")
	for feature in emblrecord.features:
		if feature.type=="CDS":
			name=""
			featureseq=feature.extract(emblrecord.seq)
			if options.inc_stop:
				proteinseq=feature.extract(emblrecord.seq).translate()
			else:
				proteinseq=feature.extract(emblrecord.seq).translate(to_stop=True)
			
			for qualifier in ["systematic_id", "locus_tag"]:
				if qualifier in feature.qualifiers:
					name=feature.qualifiers[qualifier][0]
					break
			if name!="":
				print >> output, ">"+name
				print >> output, proteinseq
			else:
				print "Cannot find name for feature"
				print feature
				sys.exit()
	
	output.close()
	
