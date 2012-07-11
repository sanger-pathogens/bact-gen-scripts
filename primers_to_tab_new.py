#!/usr/bin/env python

#Compares mummer snp taboutput files.  NOTE: must be from the same reference sequence



##################
# Import modules #
##################

import string, re
import os, sys
import tre
from optparse import OptionParser
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_general import *



##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options] fasta/multifasta input files"
	parser = OptionParser(usage=usage)

	parser.add_option("-p", "--primers", action="store", dest="primers", help="text file containing primers. For each region of interest there should be one row. Each row should contain 3 or 5 whitespace-delimited columns: 1) Name of the region, 2) forward primer sequence, 3) reverse primer sequence 4) maximum product length (optional) 5) minimum product length (optional, and can only be included with 4). The file may contain a header row, but if this is the case the -H option must be specified (or strange things might happen)", default="", metavar="FILE")
	parser.add_option("-H", "--header", action="store_true", dest="header", help="primer file has a header row [default=%default]", default=False)
	parser.add_option("-i", "--incprimer", action="store_true", dest="incprimer", help="include primer sequence in fasta output [default=%default]", default=False)
	
	parser.add_option("-m", "--maxcost", action="store", dest="maxcost", help="maximum cost (number of changes in a primer sequences) allowed for a valid match. [Default= %default]", default=3, type="int")
	parser.add_option("-o", "--prefix", action="store", dest="prefix", help="output file prefix", default="")
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):

	if options.primers=='':
		DoError('No primer file selected!')
	elif not os.path.isfile(options.primers):
		DoError('Cannot find file '+options.primers+'!')
	elif options.maxcost>10 or options.maxcost<0:
		DoError('Maximum cost must be between 1 and 10!')

		
		
	return




############################################
# Function to reverse complement sequences #
############################################

def revcomp(sequence):
	rev=sequence[::-1]
	revcomp=''
	d={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'g':'c', 'c':'g'}
	for i in rev:
		if d.has_key(i):
			revcomp=revcomp+d[i]
		else:
			revcomp=revcomp+i
	return revcomp



if __name__ == "__main__":
	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)
	
	fastafiles=args
	seqs={}
	
	for f in fastafiles:
		if not os.path.isfile(f):
			print "Cannot find file", f, "Skipping..."
			continue
		lines=open(f,'rU').read().split(">")[1:]
		for line in lines:
			seq=''.join(line.split('\n')[1:]).upper()#.replace("-","")
			seqname=line.split('\n')[0].split()[0]
			
			if seqs.has_key(seqname):
				print "Error: Two of your sequences have the same name:", name
				sys.exit()
			
			seqs[seqname]=seq

	if len(seqs.keys())==0:
		DoError('No valid sequences found')

	lines=open(options.primers, 'rU').readlines()
	
	
	if options.header:
		header=lines[0]
		lines=lines[1:]
	
	
	patterns={}
	primeroutput={}
	
	for x, line in enumerate(lines):
		words=line.strip().split()
		if len(line.strip().split())<3:
			print "Found invalid line in primer file", line.strip()
			continue
		primername=words[0]
		primer1seq=words[1].upper()
		primer2seq=words[2].upper()
		if patterns.has_key(primername):
			print "Error: Two of your regions have the same name:", name
			print "Skipping..."
			continue
		patterns[primername]=[[],[]]
		patterns[primername][0].append(tre.compile(primer1seq, tre.EXTENDED))
		patterns[primername][0].append(tre.compile(revcomp(primer1seq), tre.EXTENDED))
		patterns[primername][1].append(tre.compile(primer2seq, tre.EXTENDED))
		patterns[primername][1].append(tre.compile(revcomp(primer2seq), tre.EXTENDED))
		if len(line.strip().split())>=4:
			try:
				max_product_len=int(words[3])
			except ValueError:
				print "Invalid product length for region", primername
				max_product_len=10000
		if len(line.strip().split())>=5:
			try:
				min_product_len=int(words[4])
			except ValueError:
				print "Invalid minimum product length for region", primername
				min_product_len=0
			
		else:
			max_product_len=10000
			min_product_len=0
		
		primeroutput[primername]=open(options.prefix+primername+".fasta", "w")
		patterns[primername].append(max_product_len)
		patterns[primername].append(min_product_len)
		
	
	fz = tre.Fuzzyness(maxerr = options.maxcost)
	for seq in seqs.keys():
		taboutput=open(options.prefix+seq+"_pcr_regions.tab",'w')
		print >> taboutput, 'ID   pcr_regions'
		print seq
		for primername in patterns.keys():
			print " Searching for", primername
			
			#search for first primer forward
			matches = patterns[primername][0][0].match(seqs[seq], fz)
			mincost1=4
			primer1matches=[]
			if matches:
				mincost1=matches.cost
				for match in matches.groups():
					primer1matches.append([match,"f"])
			
			#search for first primer reerse
			matches = patterns[primername][0][1].match(seqs[seq], fz)
			if matches and matches.cost<=mincost1:
				if matches.cost<mincost1:
					primer1matches=[]
					mincost1=matches.cost
				for match in matches.groups():
					primer1matches.append([match,"r"])

			
			if len(primer1matches)==0:
				print "  No matches found for primer1. Skipping search for primer2"
				continue
			
			
			#search for second primer forward
			matches = patterns[primername][1][0].match(seqs[seq], fz)
			mincost2=4
			primer2matches=[]
			if matches:
				mincost2=matches.cost
				for match in matches.groups():
					primer2matches.append([match,"f"])
			
			#search for second primer reerse
			matches = patterns[primername][1][1].match(seqs[seq], fz)
			if matches and matches.cost<=mincost2:
				if matches.cost<mincost2:
					primer2matches=[]
					mincost2=matches.cost
				for match in matches.groups():
					primer2matches.append([match,"r"])
			

			if len(primer2matches)==0:
				print "  No matches found for primer2."
				continue
			
			
			max_product_len=patterns[primername][2]
			
			count=0
			for p1 in primer1matches:
				if p1[0][0]>p1[0][1]:
					p1max=p1[0][0]
					p1min=p1[0][1]
				else:
					p1max=p1[0][1]
					p1min=p1[0][0]
					
				for p2 in primer2matches:
					if p2[0][0]>p2[0][1]:
						p2max=p2[0][0]
						p2min=p2[0][1]
					else:
						p2max=p2[0][1]
						p2min=p2[0][0]
					
					
					if p1min>p2max and (p1min-p2max)<max_product_len and (p1min-p2max)>min_product_len and p1[1]!=p2[1]:
						if count==0:
							print >> primeroutput[primername], ">"+seq
						else:
							print >> primeroutput[primername], ">"+seq+"_"+str(count+1)
						count+=1
						if options.incprimer:
							print >> primeroutput[primername], revcomp(seqs[seq][p2min:p1max])
						else:
							print >> primeroutput[primername], revcomp(seqs[seq][p2max:p1min])
							
						print "  Found matching region on reverse strand from", p2max, "to", p1min, "with total cost of", mincost1+mincost2
						print >> taboutput, 'FT   misc_feature    complement('+str(p1min)+".."+str(p1max)+")"
						print >> taboutput, 'FT                   /primer="'+primername+' primer1"'
						print >> taboutput, 'FT                   /match_cost='+str(mincost1)
						if mincost1==0:
							print >> taboutput, 'FT                   /colour=3'
						else:
							print >> taboutput, 'FT                   /colour=1'
						print >> taboutput, 'FT   misc_feature    complement('+str(p2min)+".."+str(p2max)+")"
						print >> taboutput, 'FT                   /primer="'+primername+' primer2"'
						print >> taboutput, 'FT                   /colour=3'
						print >> taboutput, 'FT                   /match_cost='+str(mincost2)
						if mincost2==0:
							print >> taboutput, 'FT                   /colour=3'
						else:
							print >> taboutput, 'FT                   /colour=1'
						print >> taboutput, 'FT   misc_feature    complement('+str(p2max)+".."+str(p1min)+")"
						print >> taboutput, 'FT                   /primer="'+primername+'"'
						print >> taboutput, 'FT                   /colour=2'
					
					elif p2min>p1max and (p2min-p1max)<max_product_len and p1[1]!=p2[1]:
						if count==0:
							print >> primeroutput[primername], ">"+seq
						else:
							print >> primeroutput[primername], ">"+seq+"_"+str(count+1)
						count+=1
						if options.incprimer:
							print >> primeroutput[primername], seqs[seq][p1min:p2max]
						else:
							print >> primeroutput[primername], seqs[seq][p1max:p2min]
						print "  Found matching region on forward strand from", p1max, "to", p2min, "with total cost of", mincost1+mincost2
						print >> taboutput, 'FT   misc_feature    '+str(p1min)+".."+str(p1max)
						print >> taboutput, 'FT                   /primer="'+primername+' primer1"'
						print >> taboutput, 'FT                   /match_cost='+str(mincost1)
						if mincost1==0:
							print >> taboutput, 'FT                   /colour=3'
						else:
							print >> taboutput, 'FT                   /colour=1'
						print >> taboutput, 'FT   misc_feature    '+str(p2min)+".."+str(p2max)
						print >> taboutput, 'FT                   /primer="'+primername+' primer2"'
						print >> taboutput, 'FT                   /colour=3'
						print >> taboutput, 'FT                   /match_cost='+str(mincost2)
						if mincost2==0:
							print >> taboutput, 'FT                   /colour=3'
						else:
							print >> taboutput, 'FT                   /colour=1'
						print >> taboutput, 'FT   misc_feature    '+str(p1max)+".."+str(p2min)
						print >> taboutput, 'FT                   /primer="'+primername+'"'
						print >> taboutput, 'FT                   /colour=2'

		taboutput.close()

	
	for outfile in primeroutput:
		outfile.close()
	
	
	


