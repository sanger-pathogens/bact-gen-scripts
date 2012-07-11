#!/usr/bin/env python

#/usr/bin/env python


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


emblrecord=open_annotation(sys.argv[1], quiet=True)
pseudodict={}
print "\t".join(["Name", "Product", "Artemis colour", "longer/shorter", "Notes"])
for feature in emblrecord.features:
	if feature.type=="CDS":
		featureseq=feature.extract(emblrecord.seq)
		aaseq= featureseq.translate()
		stopcount=list(aaseq[:-1]).count("*")
		if stopcount>0:
			name=""
			for nametype in ["systematic_id", "gene", "locus_tag"]:
				if nametype in feature.qualifiers:
					name=feature.qualifiers[nametype][0]
					break
			if "product" in feature.qualifiers:
				product=feature.qualifiers["product"][0]
			else:
				product=""
			colour=""
			for colourtype in ["color", "colour"]:
				if colourtype in feature.qualifiers:
					colour=feature.qualifiers[colourtype][0]
					break
			print "\t".join([name, product, colour, "shorter", "contains "+str(stopcount)+" unexpected stop codon(s)"])
			pseudodict[name]=featureseq
		elif aaseq[-1] != "*":
			name=""
			for nametype in ["systematic_id", "gene", "locus_tag"]:
				if nametype in feature.qualifiers:
					name=feature.qualifiers[nametype][0]
					break
			if "product" in feature.qualifiers:
				product=feature.qualifiers["product"][0]
			else:
				product=""
			colour=""
			for colourtype in ["color", "colour"]:
				if colourtype in feature.qualifiers:
					colour=feature.qualifiers[colourtype][0]
					break
			
			print "\t".join([name, product, colour, "longer", "last codon is not a stop codon"])
			pseudodict[name]=featureseq

sys.exit()

emblrecord=open_annotation(sys.argv[2], quiet=False)
tmpname="temp"
for feature in emblrecord.features:
	if feature.type=="CDS":
		name=""
		for nametype in ["systematic_id", "gene", "locus_tag"]:
				if nametype in feature.qualifiers:
					name=feature.qualifiers[nametype][0]
					break
		if name in pseudodict:
			referenceseq=feature.extract(emblrecord.seq)
			
			output=open("temp.fasta", "w")
			print >> output, ">ref"
			print >> output, referenceseq
			print >> output, ">new"
			print >> output, pseudodict[name]
			output.close()
			
			os.system("muscle -in "+tmpname+".fasta -out "+tmpname+".aln > /dev/null 2>&1")
				
			muscleout=open(tmpname+".aln", "rU")
			os.system("seaview "+tmpname+".aln")
			curseq=""
			for line in muscleout:
				words=line.strip().split()
				if len(words)>0 and len(words[0])>0:
					if words[0][0]==">":
						curseq=words[0][1:]
						seqs[curseq]=""
					elif curseq in seqs:
						seqs[curseq]+="".join(words)
		