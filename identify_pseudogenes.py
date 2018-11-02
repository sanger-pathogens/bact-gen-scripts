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
from optparse import OptionParser, OptionGroup
sys.path.extend(map(os.path.abspath, ['/nfs/pathogen/sh16_scripts/modules/']))
from Si_general import *
from Si_SeqIO import *
from random import *



##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	
	parser.add_option("-r", "--reference", action="store", dest="ref", help="Reference embl file name", default="", metavar="FILE")
	parser.add_option("-q", "--query", action="store", dest="query", help="Query embl file name", default="", metavar="FILE")
	parser.add_option("-o", "--output", action="store", dest="output", help="Prefix for output files", default="", metavar="STRING")
	parser.add_option("-s", "--snps", action="store_true", dest="SNPs", help="Print non-pseudogenes with SNPs to separate file", default=False, metavar="STRING")
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="Print more output to screen", default=False, metavar="STRING")
	
	
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


#	if options.processors>cpu_count():
#		DoError('You do not have '+str(options.processors)+' processors')
#
	if options.ref=='':
		DoError('No reference embl file selected (-r)')
	elif not os.path.isfile(options.ref):
		DoError('Cannot find file '+options.ref)
	
	if options.query=='':
		DoError('No query embl file selected (-q)')
	elif not os.path.isfile(options.query):
		DoError('Cannot find file '+options.query)
	
	if options.output=="":
		DoError('No output prefix selected (-o)')
	
	
	return


################
# Main program #
################		

if __name__ == "__main__":

	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)


	
	refrecord=open_annotation(options.ref, quiet=True)
	refdict={}
	reffound={}
	reforder=[]
	for feature in refrecord.features:
		if feature.type=="CDS":
			name=""
			for nametype in [ "locus_tag", "systematic_id", "gene"]:
				if nametype in feature.qualifiers:
					name=feature.qualifiers[nametype][0]
					break
			if not name in refdict:
				refdict[name]=feature
				reffound[name]=False
				reforder.append(name)
			else:
				print "Repeated name:", name, "Please note that this may cause errors in the output of the script"
				
	
	
	
	print len(refdict), "CDSs found in reference embl file"
	
	pseudofile=open(options.output+"_pseudo.tsv", "w")
	if options.SNPs:
		SNPsfile=open(options.output+"_SNPs.tsv", "w")
	
	chars = string.ascii_letters + string.digits
	tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))
	bases=["A", "C", "G", "T"]
	missing=["N", "?"]
	
	emblrecord=open_annotation(options.query, quiet=True)
	pseudodict={}
	print >> pseudofile, "\t".join(["Name", "locus_tag", "Product", "Artemis colour", "pseudo qualifier in ref", "joins in ref", "joins in query", "ref length", "query length", "SNPs", "Stop codons in ref", "Stop codons in query", "Insertions", "deletions", "Notes"])
	for feature in emblrecord.features:
		if feature.type=="CDS":
			featureseq=feature.extract(emblrecord.seq)
			aaseq= featureseq.translate()
			stopcount=list(aaseq[:-1]).count("*")
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
			lname=""
			for nametype in [ "locus_tag", "systematic_id", "gene"]:
				if nametype in feature.qualifiers:
					lname=feature.qualifiers[nametype][0]
					break
			if lname in reffound:
				reffound[lname]=True
			else:
				print "Name not found in reference annotation!!", lname
				continue
			
			refjoins=len(refdict[lname].sub_features)
			if refjoins>0:
				refjoins-=1
			
			if "pseudo" in refdict[lname].qualifiers:
				refpseudo=True
			else:
				refpseudo=False
			
			
			joins=len(feature.sub_features)
			if joins>0:
				joins-=1
			
			reffeatureseq=refdict[lname].extract(refrecord.seq)
			refaaseq= reffeatureseq.translate()
			refstopcount=list(refaaseq[:-1]).count("*")
			
			
			if id(featureseq) != id(reffeatureseq):
				
				output=open(tmpname+".fasta", "w")
				print >> output, ">ref"
				print >> output, reffeatureseq
				print >> output, ">new"
				print >> output, featureseq
				output.close()
				
				os.system("muscle -in "+tmpname+".fasta -out "+tmpname+".aln > /dev/null 2>&1")
					
				muscleout=open(tmpname+".aln", "rU")
	#			os.system("seaview "+tmpname+".aln")
				curseq=""
				seqs={}
				for line in muscleout:
					words=line.strip().split()
					if len(words)>0 and len(words[0])>0:
						if words[0][0]==">":
							curseq=words[0][1:]
							seqs[curseq]=""
						elif curseq in seqs:
							seqs[curseq]+="".join(words)
				
				

				SNPs=0
				insertions=0
				deletions=0
				SNOPs=0
				STIPs=0
				indel=False
				inins=False
				for x, refbase in enumerate(seqs["ref"]):
					subbase=seqs["new"][x]
					if refbase!=subbase:
						if refbase.upper() in bases and subbase.upper() in bases:
							SNPs+=1
							indel=False
							inins=False
						elif refbase.upper() in bases and subbase.upper()=="-" and not indel:
							indel=True
							deletions+=1
							inins=False
						elif refbase.upper()=="-" and subbase.upper() in bases and not inins:
							inins=True
							insertions+=1
							indel=False
						elif (refbase.upper() not in bases and refbase.upper()!="-") or (subbase.upper() not in bases and subbase.upper()!="-"):
							print refbase.upper(), subbase.upper()
					else:
						indel=False
						inins=False
	#			print SNPs, insertions, deletions, stopcount, refstopcount
	
				notes=[]
				outlist=[name, lname, product, colour, refpseudo, refjoins, joins, len(refdict[lname]), len(feature), SNPs, refstopcount, stopcount, insertions, deletions]
				pseudogene=False
				SNPgene=False
				if refjoins<joins:
					notes.append("EXTRA JOIN(S)")
					pseudogene=True
				elif joins>0:
					notes.append("HAS JOIN(S)")
					pseudogene=True
				elif refjoins>joins:
					notes.append("LOST ALL JOIN(S)")
					pseudogene=True
				if len(refdict[lname])!=len(feature):
					notes.append("CHANGED LENGTH")
					pseudogene=True
				if insertions>0 :
					notes.append("INSERTIONS")
					pseudogene=True
				if deletions >0:
					notes.append("DELETIONS")
					pseudogene=True
				if SNPs>0:
					notes.append("SNPS")
					SNPgene=True
				if stopcount>0:
					notes.append("STOP CODON(S) IN CDS")
					pseudogene=True
				if refjoins==joins and joins==0 and len(refdict[lname])==len(feature) and insertions==0 and deletions==0 and (SNPs==0 or not options.SNPs) and stopcount==0:
					notes.append("NOT PSEUDOGENE")
						
					
				outlist.append(', '.join(notes))
								
				if pseudogene:
					print >> pseudofile, "\t".join(map( str, outlist))
					if name==lname:
						print lname, product, 
					else:
						print lname, name, product,
					print ":", ', '.join(notes)
				elif options.SNPs and SNPgene:
					print >> SNPsfile, "\t".join(map( str, outlist))
					if name==lname:
						print lname, product, 
					else:
						print lname, name, product,
					print ":", ', '.join(notes)
				elif options.verbose:
					if name==lname:
						print lname, product, 
					else:
						print lname, name, product,
					print ":", ', '.join(notes)
				
			else:
				if options.verbose:
					print ": CONSERVED"
				continue
	
	pseudofile.close()
	if options.SNPs:
		SNPsfile.close()
	
	missingfile=open(options.output+"_missing.txt", "w")
	print >> missingfile, "\nCDSs in reference missing in subject:"
	count=0
	for name in reforder:
		if not reffound[name]:
			print >> missingfile, name
			count+=1
	print count, "CDSs in reference annotation missing in query annotation"
	missingfile.close()
	
	os.system("rm -f "+tmpname+"*")
	
	
