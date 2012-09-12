#!/usr/bin/env python

#Create plot of dissimilarity between pairs of strains in an alignment



#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------
import string, re
import os, sys
import dendropy
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/lib/python2.7/site-packages/fisher-0.1.4-py2.7-linux-x86_64.egg/', '/nfs/users/nfs_s/sh16/lib/python2.5/site-packages/fisher-0.1.4-py2.7-linux-x86_64.egg/', '/nfs/users/nfs_s/sh16/scripts/modules/']))
from fisher import pvalue
from optparse import OptionParser
from Bio.Align import AlignInfo
#sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_general import *
from Si_SeqIO import *

##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options] args"
	parser = OptionParser(usage=usage)

	parser.add_option("-a", "--alignment", action="store", dest="alignment", help="alignment file name. For the tab output to fit you need to give the entire alignment, not just SNP sites, which is faster.", default="", metavar="FILE")
	parser.add_option("-t", "--tree", action="store", dest="tree", help="tree file name", default="", metavar="FILE")
	parser.add_option("-p", "--proportion", action="store", dest="proportion", help="maximum proportion of Ns to allow in a column for it to be included in the test (i.e. ignore any sites with > than this proportion of Ns) [default=%default]", default=0.05)
	parser.add_option("-r", "--reference", action="store", dest="reference", help="Name of reference in alignment", default="")
	
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
	if options.proportion>1 or options.proportion<0:
		DoError('Maximum proportion of Ns must be between 0 and 1!')
	if options.reference=="":
		DoError('No reference selected!')

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
	diffNpositions=set([])
	
	for x in xrange(len(taxonseq)):
		if x in toremove:
			continue
		if taxonseq[x] in missing_set:
			Ncount+=1
		if taxonbseq[x] in gap_and_missing or taxoncseq[x] in gap_and_missing:
			continue
		elif taxonbseq[x]==taxoncseq[x]:
			if taxonseq[x] in missing_set:
				sameN+=1
			elif taxonseq[x] not in missing_set:
				samenotN+=1
		elif taxonbseq[x]!=taxoncseq[x]:
			if taxonseq[x] in missing_set:
				diffN+=1
				diffNpositions.add(x)
			elif taxonseq[x] not in missing_set:
				diffnotN+=1
	
	#matrix=[[sameN, samenotN],[diffN, diffnotN]]
	p = pvalue(diffN, diffnotN, sameN, samenotN)
#	if p.right_tail<(0.05/tot):
#		print taxon, taxonb, taxonc, Ncount, sameN, diffN, samenotN, diffnotN, p
#	print p.right_tail
	return p.right_tail	, diffNpositions




################
# Main program #
################		

if __name__ == "__main__":
	

	
	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)
	
	
	#Read the alignment file
	
	
	#first find the reference
	reffound=False
	refseq=[]
	print "Finding reference in alignment"
	for line in open(options.alignment, "rU"):
		line=line.strip()
		if len(line)>0 and line[0]==">":
			if line.split()[0][1:]==options.reference:
				reffound=True
			elif reffound:
				break
		elif reffound:
			refseq.append(''.join(line.split()))
		
	
	if not reffound:
		DoError("Cannot find reference in alignment")
	else:
		print "Reference found"
	refseq=''.join(refseq).upper()
	
	print "Finding SNP sites"
	SNPsites=set([])
	seq=[]
	
	ok_bases=set(['A','G','C','T'])
	
	for line in open(options.alignment, "rU"):
		line=line.strip()
		if len(line)>0 and line[0]==">":
			if line.split()[0][1:]==options.reference:
				seq=[]
				continue
			seq=''.join(seq).upper()
			if len(seq)==len(refseq):
				for base in xrange(len(refseq)):
					if base in SNPsites:
						continue
					if refseq[base]!=seq[base] and refseq[base] in ok_bases and seq[base] in ok_bases:
						SNPsites.add(base)
		
			seq=[]
		elif reffound:
			seq.append(''.join(line.split()))
	if line.split()[0][1:]!=options.reference:
		seq=''.join(seq).upper()
		if len(seq)==len(refseq):
			for base in xrange(len(refseq)):
				if base in SNPsites:
					continue
				if refseq[base]!=seq[base] and refseq[base] in ok_bases and seq[base] in ok_bases:
					SNPsites.add(base)
	
	print "Found", len(SNPsites), "SNP sites"
	
	sorted_snps=list(SNPsites)
	sorted_snps.sort()
	seqnames=[]
	
	pairwise_distances={}
	tot=0
	
	alignment={}
	
	refpos={}
	
	print "Extracting SNP sites and their positions in reference"
	name=""
	seq=[]
	for line in open(options.alignment, "rU"):
		line=line.strip()
		if len(line)>0 and line[0]==">":
			if name!="":
				seqnames.append(name)
				pairwise_distances[name]={}
				seq=''.join(seq).upper()
				snpseq=[]
				for x in sorted_snps:
					snpseq.append(seq[x])
				
				if name==options.reference:
					count=0
					SNPcount=0
					for base in xrange(len(refseq)):
						
						if refseq[base]!="-":
							count+=1
						if base in SNPsites:
							refpos[SNPcount]=count
							SNPcount+=1

				alignment[name]=''.join(snpseq)
				seq=[]
				
			name=line.split()[0][1:]
		elif reffound:
			seq.append(''.join(line.split()))

	if name!="":
		seqnames.append(name)
		pairwise_distances[name]={}
		seq=''.join(seq).upper()
		snpseq=[]
		for x in sorted_snps:
			snpseq.append(seq[x])
		
		if name==options.reference:
			count=0
			SNPcount=0
			for base in xrange(len(refseq)):
				
				if refseq[base]!="-":
					count+=1
				if base in SNPsites:
					refpos[SNPcount]=count
					SNPcount+=1
		alignment[name]=''.join(snpseq)

	
#	print alignment
#	print seqnames
#	print pairwise_distances
#	print refpos
#	sys.exit()	
#	
	
	
#	try:
#		readalignment=read_alignment(options.alignment)
#	except StandardError:
#		DoError("Cannot open alignment file")
#
#	seqnames=[]
#	
#	pairwise_distances={}
#	tot=0
#	
#	alignment={}
#	for sequence in readalignment:
#		alignment[sequence.name]=str(sequence.seq).upper()
#		if not sequence.name in seqnames:
#			seqnames.append(sequence.name)
#			pairwise_distances[sequence.name]={}
	
	alnlen=len(alignment[name])
	taxacount=len(alignment)
	seqcount=len(seqnames)
	toremove=set([])
	cutoff=2
	if float(taxacount)*options.proportion>cutoff:
		cutoff=float(taxacount)*options.proportion
	
	print "Cutoff set at", cutoff	
	
	for x in xrange(alnlen):
		ncount=0.0
		for seq in seqnames:
			if alignment[seq][x] in gap_and_missing:
				ncount+=1
		if ncount>cutoff:
			toremove.add(x)
	
	
	
#	for x in toremove[::-1]:
#		for seq in seqnames:
#			alignment[seq]=alignment[seq][:x]+alignment[seq][x+1:]
#	alnlen=len(alignment[name])
#	
	
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
					if y in toremove:
						continue
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
	Nlocsets={}
	
	for x, taxon in enumerate(seqnames):
		print x+1, "\r",
		sys.stdout.flush()
		for y, taxonb in enumerate(seqnames[x+1:]):
			
			sys.stdout.flush()
			for taxonc in seqnames[x+y+2:]:
				if pairwise_distances[taxon][taxonb]+pairwise_distances[taxon][taxonc]<pairwise_distances[taxonb][taxonc]:
					cspvalue, Nlocset=chisquarepvalue(taxon,taxonb,taxonc)
					if cspvalue<(0.05/tot):
						if not taxon in Nlocsets:
							Nlocsets[taxon]={}
						if not taxonb in Nlocsets[taxon]:
							Nlocsets[taxon][taxonb]={}
						if not taxonc in Nlocsets[taxon][taxonb]:
							Nlocsets[taxon][taxonb][taxonc]=Nlocset
					
						if not taxon in significant:
							significant[taxon]={}
						if not taxonb in significant[taxon]:
							significant[taxon][taxonb]=0
						significant[taxon][taxonb]+=1
						if not taxonc in significant[taxon]:
							significant[taxon][taxonc]=0
						significant[taxon][taxonc]+=1
				elif pairwise_distances[taxon][taxonb]+pairwise_distances[taxonb][taxonc]<pairwise_distances[taxon][taxonc]:
					cspvalue, Nlocset=chisquarepvalue(taxonb,taxon,taxonc)
					if cspvalue<(0.05/tot):
						if not taxonb in Nlocsets:
							Nlocsets[taxonb]={}
						if not taxon in Nlocsets[taxonb]:
							Nlocsets[taxonb][taxon]={}
						if not taxonc in Nlocsets[taxonb][taxon]:
							Nlocsets[taxonb][taxon][taxonc]=Nlocset
							
						if not taxonb in significant:
							significant[taxonb]={}
						if not taxon in significant[taxonb]:
							significant[taxonb][taxon]=0
						significant[taxonb][taxon]+=1
						if not taxonc in significant[taxonb]:
							significant[taxonb][taxonc]=0
						significant[taxonb][taxonc]+=1
				elif pairwise_distances[taxon][taxonc]+pairwise_distances[taxonb][taxonc]<pairwise_distances[taxon][taxonb]:
					cspvalue, Nlocset=chisquarepvalue(taxonc,taxonb,taxon)
					if cspvalue<(0.05/tot):
						if not taxonc in Nlocsets:
							Nlocsets[taxonc]={}
						if not taxonb in Nlocsets[taxonc]:
							Nlocsets[taxonc][taxon]={}
						if not taxonc in Nlocsets[taxonc][taxon]:
							Nlocsets[taxonc][taxon][taxonb]=Nlocset
							
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
	logout=open("heterozygotes.log","w")
	print >> logout, "Taxon", "Number of significant comparisons"
	for taxon in seqnames:
		if taxon in significant:
			
			outfile=open("test.csv","w")
			print >> outfile, "name,type,number:c:1:"+str(len(seqnames)-1)
			print >> outfile, taxon+","+"1,"+"-"
			maxnum=0
			for taxonb in significant[taxon]:
				print >> outfile, taxonb+","+"-,"+str(significant[taxon][taxonb])
				if significant[taxon][taxonb]>maxnum:
					maxnum=significant[taxon][taxonb]

				
			outfile.close()
			print >> logout, taxon, float(maxnum)/(len(seqnames)-1)
			os.system("/nfs/users/nfs_s/sh16/scripts/reportlabtest.py -t "+options.tree+" -m test.csv -C 2,3 -O portrait -M -L left -a 2 -o "+taxon+"_"+str(maxnum)+".pdf")
			os.system("rm -f  test.csv")
		else:
			print >> logout, taxon, 0.0
	
	
	for taxon in Nlocsets:
		Nlocs={}
		for taxonb in Nlocsets[taxon]:
			for taxonc in Nlocsets[taxon][taxonb]:
				for base in Nlocsets[taxon][taxonb][taxonc]:
					if not base in Nlocs:
						Nlocs[base]=[]
					Nlocs[base].append('-'.join([taxonb,taxonc]))
		tabout=open(taxon+"_Ns.tab", "w")
		for Nloc in Nlocs:
			print >> tabout, "FT   misc_feature    "+str(refpos[Nloc])
			print >> tabout, 'FT                   /taxa="'+', '.join(Nlocs[Nloc])+'"'
		tabout.close()
				
	
	sys.exit()
                                    
