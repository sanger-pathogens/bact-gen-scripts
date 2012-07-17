#!/usr/bin/env python

#Create plot of dissimilarity between pairs of strains in an alignment



#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------
import dendropy
from fisher import pvalue
import string, re
import os, sys
from optparse import OptionParser
from Bio.Align import AlignInfo
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_general import *
from Si_SeqIO import *

##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options] args"
	parser = OptionParser(usage=usage)

	parser.add_option("-a", "--alignment", action="store", dest="alignment", help="alignment file name", default="", metavar="FILE")
	parser.add_option("-t", "--tree", action="store", dest="tree", help="tree file name", default="", metavar="FILE")
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


#	if options.embl!="" and options.reference=='':
#		DoError('No reference selected! If you give an embl file you must specify which sequence it is linked to (i.e. the reference)')
	if options.alignment=='':
		DoError('No alignment file selected!')
	elif not os.path.isfile(options.alignment):
		DoError('Cannot find file '+options.alignment+'!')
	if options.tree=='':
		DoError('No tree file selected!')
	elif not os.path.isfile(options.tree):
		DoError('Cannot find file '+options.tree+'!')

gap_and_missing=set(["-", "N", "?"])
missing_set=([ "N", "?"])

def chisquarepvalue(taxon,taxonb,taxonc):
	
	taxonseq=alignment[taxon]
	taxonbseq=alignment[taxonb]
	taxoncseq=alignment[taxonc]
	
	diffN=0.0
	sameN=0.0
	diffnotN=0.0
	samenotN=0.0
	Ncount=0.0
	
	for x in xrange(len(taxonseq)):
		if taxonseq[x] in missing_set:
			Ncount+=1
		if taxonbseq[x] in gap_and_missing or taxoncseq[x] in gap_and_missing:
			continue
		elif taxonbseq[x]==taxoncseq[x]:
			if taxonseq[x] in missing_set:
				sameN+=1
			elif taxonseq[x]!="-":
				samenotN+=1
		elif taxonbseq[x]!=taxoncseq[x]:
			if taxonseq[x] in missing_set:
				diffN+=1
			elif taxonseq[x]!="-":
				diffnotN+=1
	
	#matrix=[[sameN, samenotN],[diffN, diffnotN]]
	p = pvalue(diffN, diffnotN, sameN, samenotN)
#	if p.right_tail<(0.05/tot):
#		print taxon, taxonb, taxonc, Ncount, sameN, diffN, samenotN, diffnotN, p
#	print p.right_tail
	return p.right_tail	




################
# Main program #
################		

if __name__ == "__main__":
	

	
	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)
	
	
	#Read the alignment file
	
	
	try:
		readalignment=read_alignment(options.alignment)
	except StandardError:
		DoError("Cannot open alignment file")

	seqnames=[]
	
	pairwise_distances={}
	tot=0
	
	alignment={}
	for sequence in readalignment:
		alignment[sequence.name]=str(sequence.seq).upper()
		if not sequence.name in seqnames:
			seqnames.append(sequence.name)
			pairwise_distances[sequence.name]={}
	
	alnlen=len(alignment[sequence.name])
	seqcount=len(seqnames)
	toremove=[]
	cutoff=1
	if float(alnlen)/10>cutoff:
		cutoff=float(alnlen)/10
	for x in xrange(alnlen):
		ncount=0.0
		for seq in seqnames:
			if alignment[seq][x] in gap_and_missing:
				ncount+=1
		if ncount>cutoff:
			toremove.append(x)
	
	for x in toremove[::-1]:
		for seq in seqnames:
			alignment[seq]=alignment[seq][:x]+alignment[seq][x+1:]
	alnlen=len(alignment[sequence.name])
	
	for x, taxon in enumerate(seqnames):
		pairwise_distances[taxon]={}
		taxonseq=""
		print x+1, "\r",
		sys.stdout.flush()
		try: taxonseq=alignment[taxon]
		except StandardError:
			print "Cannot find ", taxon
			continue
#		print x
		for taxonb in seqnames[x+1:]:
						
			if taxon==taxonb:
				continue
			
			try: taxonbseq=alignment[taxonb]
			except StandardError:
				print "Cannot find ", taxonb
				continue
			
			count=0
			ncount=0			
			
			if len(taxonseq)==len(taxonbseq):
			
				for y, base in enumerate(taxonseq):
					if base!=taxonbseq[y] and base not in gap_and_missing and taxonbseq[y] not in gap_and_missing:
						count+=1
						
						#print y+1, base, taxonbseq[y]
					elif base=="N":
						ncount+=1
						
			else:
				print taxon, "and", taxonb, "are not the same length"
			pairwise_distances[taxon][taxonb]=count
			pairwise_distances[taxonb][taxon]=count
			tot+=1
			#print taxon, taxonb, count, ncount

#	print "\n"
	
	significant={}
	for x, taxon in enumerate(seqnames):
		print x+1, "\r",
		sys.stdout.flush()
		for y, taxonb in enumerate(seqnames[x+1:]):
			
			sys.stdout.flush()
			for taxonc in seqnames[x+y+2:]:
				if pairwise_distances[taxon][taxonb]+pairwise_distances[taxon][taxonc]<pairwise_distances[taxonb][taxonc]:
					cspvalue=chisquarepvalue(taxon,taxonb,taxonc)
					if cspvalue<(0.05/tot):
						if not taxon in significant:
							significant[taxon]={}
						if not taxonb in significant[taxon]:
							significant[taxon][taxonb]=0
						significant[taxon][taxonb]+=1
						if not taxonc in significant[taxon]:
							significant[taxon][taxonc]=0
						significant[taxon][taxonc]+=1
				elif pairwise_distances[taxon][taxonb]+pairwise_distances[taxonb][taxonc]<pairwise_distances[taxon][taxonc]:
					cspvalue=chisquarepvalue(taxonb,taxon,taxonc)
					if cspvalue<(0.05/tot):
						if not taxonb in significant:
							significant[taxonb]={}
						if not taxon in significant[taxonb]:
							significant[taxonb][taxon]=0
						significant[taxonb][taxon]+=1
						if not taxonc in significant[taxonb]:
							significant[taxonb][taxonc]=0
						significant[taxonb][taxonc]+=1
				elif pairwise_distances[taxon][taxonc]+pairwise_distances[taxonb][taxonc]<pairwise_distances[taxon][taxonb]:
					cspvalue=chisquarepvalue(taxonc,taxonb,taxon)
					if cspvalue<(0.05/tot):
						if not taxonc in significant:
							significant[taxonc]={}
						if not taxonb in significant[taxonc]:
							significant[taxonc][taxonb]=0
						significant[taxonc][taxonb]+=1
						if not taxon in significant[taxonc]:
							significant[taxonc][taxon]=0
						significant[taxonc][taxon]+=1
#				else:
#					print >> outfile, taxonb+","+"2"
#					print >> outfile, taxonc+","+"2"
	print
	print "Taxon", "Number of significant comparisons"
	for taxon in significant:
		print taxon, len(significant[taxon])
		outfile=open("test.csv","w")
		print >> outfile, "name,type,number:c:1:"+str(len(seqnames)-1)
		print >> outfile, taxon+","+"1,"+"-"
		maxnum=0
		for taxonb in significant[taxon]:
			print >> outfile, taxonb+","+"-,"+str(significant[taxon][taxonb])
			if significant[taxon][taxonb]>maxnum:
				maxnum=significant[taxon][taxonb]
		outfile.close()
		os.system("~/scripts/reportlabtest.py -t "+options.tree+" -m test.csv -C 2,3 -O portrait -M -L left -a 2 -o "+taxon+"_"+str(maxnum)+".pdf")
	sys.exit()
                                    
