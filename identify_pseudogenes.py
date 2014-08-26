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



refrecord=open_annotation(sys.argv[2], quiet=True)
refdict={}
reffound={}
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
		else:
			print "Repeated name:", name, "Please note that this may cause errors in the output of the script"
			



print len(refdict), "CDSs found in reference embl file"

pseudofile=open(sys.argv[3]+"_pseudo.tsv", "w")
missingfile=open(sys.argv[3]+"_missing.tsv", "w")
SNPsfile=open(sys.argv[3]+"_SNPs.tsv", "w")

tmpname="temp12345"
bases=["A", "C", "G", "T"]
missing=["N", "?"]

emblrecord=open_annotation(sys.argv[1], quiet=True)
pseudodict={}
print >> pseudofile, "\t".join(["Name", "Product", "Artemis colour", "pseudo qualifier in ref", "joins in ref", "joins in subject", "ref length", "subject length", "SNPs", "SNOPs", "STIPs", "Insertions", "deletions" "Notes"])
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
		
		if name==lname:
			print lname, product, 
		else:
			print lname, name, product, 
		
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
			
			
			
#			print featureseq
#			print reffeatureseq
#			print name, lname,
#			print len(refdict[lname]), len(feature),
#			print refjoins, refpseudo, joins,
#			print feature.sub_features,
#			print refdict[lname].sub_features,
#			print feature
#			print refdict[lname]
#			print seqs
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
			if refjoins!=joins or len(refdict[lname])!=len(feature) or insertions>0 or deletions >0:
				print >> pseudofile, "\t".join(map( str, [name, product, colour, refpseudo, refjoins, joins, len(refdict[lname]), len(feature), SNPs, stopcount, refstopcount, insertions, deletions]))
				print "pseudogene or changed length"
			elif SNPs>0:
				print >> SNPsfile, "\t".join(map( str, [name, product, colour, refpseudo, refjoins, joins, len(refdict[lname]), len(feature), SNPs, stopcount, refstopcount, insertions, deletions]))
				print "variant"
			else:
				print "conserved"
		else:
			print "conserved"
			continue
		
print "\nCDSs in reference missing in subject:"
for name in reffound:
	if not reffound[name]:
		print >> missingfile, name

pseudofile.close()
SNPsfile.close()
missingfile.close()

#		hasjoin=0
#		
#		if len(feature.sub_features)>0:
#			hasjoin=len(feature.sub_features)
##			for x, join in enumerate(feature.sub_features):
#				
#		
#		
#		
#		if stopcount>0:
#			if hasjoin>0:
#				print "\t".join([name, product, colour, str(refjoins), str(refpseudo), "shorter", "contains "+str(stopcount)+" unexpected stop codon(s) and "+str(hasjoin-1)+" frame shift(s)"])
#			else:
#				print "\t".join([name, product, colour, str(refjoins), str(refpseudo), "shorter", "contains "+str(stopcount)+" unexpected stop codon(s)"])
#			pseudodict[name]=featureseq
#		elif aaseq[-1] != "*":
#			print "\t".join([name, product, colour, str(refjoins), str(refpseudo), "longer", "last codon is not a stop codon"])
#			pseudodict[name]=featureseq
#		elif hasjoin>0:
#			print "\t".join([name, product, colour, str(refjoins), str(refpseudo), "shorter", "contains "+str(hasjoin-1)+" frame shift(s)"])
#		
#
#sys.exit()
#
#emblrecord=open_annotation(sys.argv[2], quiet=False)
#tmpname="temp"
#for feature in emblrecord.features:
#	if feature.type=="CDS":
#		name=""
#		for nametype in ["systematic_id", "gene", "locus_tag"]:
#				if nametype in feature.qualifiers:
#					name=feature.qualifiers[nametype][0]
#					break
#		if name in pseudodict:
#			referenceseq=feature.extract(emblrecord.seq)
#			
#			output=open("temp.fasta", "w")
#			print >> output, ">ref"
#			print >> output, referenceseq
#			print >> output, ">new"
#			print >> output, pseudodict[name]
#			output.close()
#			
#			os.system("muscle -in "+tmpname+".fasta -out "+tmpname+".aln > /dev/null 2>&1")
#				
#			muscleout=open(tmpname+".aln", "rU")
#			os.system("seaview "+tmpname+".aln")
#			curseq=""
#			for line in muscleout:
#				words=line.strip().split()
#				if len(words)>0 and len(words[0])>0:
#					if words[0][0]==">":
#						curseq=words[0][1:]
#						seqs[curseq]=""
#					elif curseq in seqs:
#						seqs[curseq]+="".join(words)
#		