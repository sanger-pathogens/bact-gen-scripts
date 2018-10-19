#!/usr/bin/env python

#Compares mummer snp taboutput files.  NOTE: must be from the same reference sequence



##################
# Import modules #
##################

import string, re
import os, sys
from optparse import OptionParser
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_general import *
sys.path.extend(map(os.path.abspath, ['/software/pathogen/external/lib/python2.7/site-packages/']))
import tre



def replace_ambiguity_codes(sequence):
	codes={"A":"A", "C":"C", "G":"G", "T":"T", "R":"[AG]", "S":"[GC]", "B":"[CGT]", "Y":"[CT]", "W":"[AT]", "D":"[AGT]", "K":"[GT]", "N":"[ACGT]", "H":"[ACT]", "M":"[AC]", "V":"[ACG]", "X": "[ACGT]"}
	regexlist=[]
	for base in sequence:
		if not base in codes:
			print "Unidentified nucleotide code in primer:", base
			sys.exit()
		else:
			regexlist.append(codes[base])
	
	return ''.join(regexlist)
		




##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options] fasta/multifasta input files"
	parser = OptionParser(usage=usage)

	parser.add_option("-p", "--primers", action="store", dest="primers", help="text file containing primers. For each region of interest there should be one row. Each row should contain 3 or 5 whitespace-delimited columns: 1) Name of the region, 2) forward primer sequence, 3) reverse primer sequence 4) maximum product length (optional) 5) minimum product length (optional, and can only be included with 4). The file may contain a header row, but if this is the case the -H option must be specified (or strange things might happen)", default="", metavar="FILE")
	parser.add_option("-H", "--header", action="store_true", dest="header", help="primer file has a header row [default=%default]", default=False)
	parser.add_option("-r", "--rcreverse", action="store_true", dest="rcreverse", help="Reverse primers require reverse complementing (use this option if your reerse primers are in the same direction as your forward primers) [default=%default]", default=False)
	parser.add_option("-i", "--incprimer", action="store_true", dest="incprimer", help="include primer sequence in fasta output [default=%default]", default=False)
	
	parser.add_option("-m", "--maxcost", action="store", dest="maxcost", help="maximum cost (number of changes in a primer sequences) allowed for a valid match. [Default= %default]", default=3, type="int")
	parser.add_option("-o", "--prefix", action="store", dest="prefix", help="output file prefix", default="")
	parser.add_option("-e", "--potentials", action="store_true", dest="potentials", help="include potential primer locations in tab files in cases where no pair within selected product size range is found [default=%default]", default=False)
	parser.add_option("-c", "--circular", action="store_true", dest="circular", help="sequences are circular [default=%default]", default=False)
	
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
	if options.prefix!="":
		options.prefix=options.prefix+"_"

		
		
	return




############################################
# Function to reverse complement sequences #
############################################

def revcomp(sequence):
	rev=sequence[::-1]
	revcomp=''
	d={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'g':'c', 'c':'g', "[":"]", "]":"["}
	for i in rev:
		if d.has_key(i):
			revcomp=revcomp+d[i]
		else:
			revcomp=revcomp+i
	return revcomp



if __name__ == "__main__":

	print"WARNING: Please use in_silico_pcr.py instead. This script is no longer supported."

	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)
	
	fastafiles=args
	seqs={}
	seqorder=[]
	filenum=-1
	
	print "Reading sequence files"
	
	for f in fastafiles:
		if not os.path.isfile(f):
			print "Cannot find file", f, "Skipping..."
			continue
		lines=open(f,'rU').read().split(">")[1:]
		seqorder.append([f.split("/")[-1].split(".")[0],[]])
		filenum+=1
		for line in lines:
			seq=''.join(line.split('\n')[1:]).upper()#.replace("-","")
			seqname=line.split('\n')[0].split()[0]
			
			if seqs.has_key(seqname):
				print "Error: Two of your sequences have the same name:", seqname, "Skipping"
				continue
				sys.exit()
			seqorder[filenum][1].append(seqname)
			seqs[seqname]=seq
	
	
	
	if len(seqs.keys())==0:
		DoError('No valid sequences found')
	else:
		print "Found", len(seqs.keys()), "sequences in", len(seqorder), "files"
	lines=open(options.primers, 'rU').readlines()
	
	
	if options.header:
		header=lines[0]
		lines=lines[1:]
		primer1name=header.split()[1]
		primer2name=header.split()[2]
	else:
		primer1name="primer1"
		primer2name="primer2"
	
	
	patterns={}
	primeroutput={}
	
	print "Reading primer file"
	
	for x, line in enumerate(lines):
		words=line.strip().split()
		if len(words)==0:
			continue
		elif len(words)<3:
			print "Found invalid line in primer file", line.strip()
			continue
		primername=words[0]
		primer1seq=replace_ambiguity_codes(words[1].upper())
		primer2seq=replace_ambiguity_codes(words[2].upper())
		if options.rcreverse:
			primer2seq=revcomp(primer2seq)
	
		
		if patterns.has_key(primername):
			print "Error: Two of your regions have the same name:", primername
			print "Skipping..."
			continue
		
		patterns[primername]=[[],[]]
		patterns[primername][0].append(tre.compile(primer1seq, tre.EXTENDED))
		patterns[primername][0].append(tre.compile(revcomp(primer1seq), tre.EXTENDED))
		patterns[primername][1].append(tre.compile(primer2seq, tre.EXTENDED))
		patterns[primername][1].append(tre.compile(revcomp(primer2seq), tre.EXTENDED))
		
		
		if len(words)>=5:
			try:
				min_product_len=int(words[4])
			except ValueError:
				print "Invalid minimum product length for region", primername
				min_product_len=0
		if len(words)>=4:
			try:
				max_product_len=int(words[3])
			except ValueError:
				print "Invalid product length for region", primername
				max_product_len=10000
			min_product_len=0
		else:
			max_product_len=10000
			min_product_len=0
		
		primeroutput[primername]=open(options.prefix+primername+".fasta", "w")
		patterns[primername].append(max_product_len)
		patterns[primername].append(min_product_len)
	
	print "Found", len(patterns), "primer pairs"
	
	fz = tre.Fuzzyness(maxerr = options.maxcost)
	
	
	
	for file in seqorder:
	
		alltaboutput=open(options.prefix+file[0]+"_"+"pcr_regions.tab",'w')
		print >> alltaboutput, 'ID   pcr_regions'
		
		print "Searching file: "+file[0]
		
		#foundprimers=[]
		
		totallength=0
		for x, seq in enumerate(file[1]):
#			taboutput=open(options.prefix+seq+"_pcr_regions.tab",'w')
#			print >> taboutput, 'ID   pcr_regions'
			#print x, seq, patterns.keys()
			for primername in patterns.keys():
				
				#print "Searching for", primername
				
				#search for first primer forward
				#print "search for first primer forward"
				matches = patterns[primername][0][0].search(seqs[seq], fz)
				mincost1=4
				primer1matches=[]
				if matches:
					mincost1=matches.cost
					for match in matches.groups():
						primer1matches.append([match,"f"])
				
				#search for first primer reverse
				#print "search for first primer reverse"
				matches = patterns[primername][0][1].search(seqs[seq], fz)
				if matches and matches.cost<=mincost1:
					if matches.cost<mincost1:
						primer1matches=[]
						mincost1=matches.cost
					for match in matches.groups():
						primer1matches.append([match,"r"])
	
				
#				if len(primer1matches)==0:
#					print "  No matches found for primer1. Skipping search for primer2"
#					continue
				
				
				#search for second primer forward
				#print "search for second primer forward"
				matches = patterns[primername][1][0].search(seqs[seq], fz)
				mincost2=4
				primer2matches=[]
				if matches:
					mincost2=matches.cost
					for match in matches.groups():
						primer2matches.append([match,"f"])
				
				#search for second primer reverse
				#print "search for second primer reverse"
				matches = patterns[primername][1][1].search(seqs[seq], fz)
				if matches and matches.cost<=mincost2:
					if matches.cost<mincost2:
						primer2matches=[]
						mincost2=matches.cost
					for match in matches.groups():
						primer2matches.append([match,"r"])
				
	
				#if len(primer2matches)==0:
				#	print "  No matches found for primer2."
#					continue
				
				
				max_product_len=patterns[primername][2]
				
				count=0
				found=False
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
						
						if p1min>=p2max and (p1min-p2max)<=patterns[primername][2] and (p1min-p2max)>=patterns[primername][3] and p1[1]!=p2[1]:
							found=True
							if count==0:
								print >> primeroutput[primername], ">"+seq
							else:
								print >> primeroutput[primername], ">"+seq+"_"+str(count+1)
							count+=1
							if options.incprimer:
								print >> primeroutput[primername], revcomp(seqs[seq][p2min:p1max])
							else:
								print >> primeroutput[primername], revcomp(seqs[seq][p2max:p1min])
								
							print "  Found matching region for "+primername+" on reverse strand of contig "+seq+" from", p2max, "to", p1min, "with total cost of", mincost1+mincost2
							
							
							print >> alltaboutput, 'FT   misc_feature    complement('+str(p1min+totallength+1)+".."+str(p1max+totallength)+")"
							print >> alltaboutput, 'FT                   /sequence='+seq
							print >> alltaboutput, 'FT                   /primer="'+primername+' '+primer1name+'"'
							print >> alltaboutput, 'FT                   /match_cost='+str(mincost1)
							if mincost1==0:
								print >> alltaboutput, 'FT                   /colour=3'
							else:
								print >> alltaboutput, 'FT                   /colour=1'
							print >> alltaboutput, 'FT   misc_feature    complement('+str(p2min+totallength+1)+".."+str(p2max+totallength)+")"
							print >> alltaboutput, 'FT                   /sequence='+seq
							print >> alltaboutput, 'FT                   /primer="'+primername+' '+primer2name+'"'
							print >> alltaboutput, 'FT                   /match_cost='+str(mincost2)
							if mincost2==0:
								print >> alltaboutput, 'FT                   /colour=3'
							else:
								print >> alltaboutput, 'FT                   /colour=1'
							print >> alltaboutput, 'FT   misc_feature    complement('+str(p2max+totallength+1)+".."+str(p1min+totallength)+")"
							print >> alltaboutput, 'FT                   /sequence='+seq
							print >> alltaboutput, 'FT                   /product="'+primername+'"'
							print >> alltaboutput, 'FT                   /colour=2'
						
						elif p2min>=p1max and (p2min-p1max)<=patterns[primername][2] and p1[1]!=p2[1]:
							found=True
							if count==0:
								print >> primeroutput[primername], ">"+seq
							else:
								print >> primeroutput[primername], ">"+seq+"_"+str(count+1)
							count+=1
							if options.incprimer:
								print >> primeroutput[primername], seqs[seq][p1min:p2max]
							else:
								print >> primeroutput[primername], seqs[seq][p1max:p2min]
							print "  Found matching region for "+primername+" on forward strand of contig "+seq+" from", p1max, "to", p2min, "with total cost of", mincost1+mincost2
							
							
							
							print >> alltaboutput, 'FT   misc_feature    '+str(p1min+totallength+1)+".."+str(p1max+totallength)
							print >> alltaboutput, 'FT                   /sequence='+seq
							print >> alltaboutput, 'FT                   /primer="'+primername+' '+primer1name+'"'
							print >> alltaboutput, 'FT                   /match_cost='+str(mincost1)
							if mincost1==0:
								print >> alltaboutput, 'FT                   /colour=3'
							else:
								print >> alltaboutput, 'FT                   /colour=1'
							print >> alltaboutput, 'FT   misc_feature    '+str(p2min+totallength+1)+".."+str(p2max+totallength)
							print >> alltaboutput, 'FT                   /sequence='+seq
							print >> alltaboutput, 'FT                   /primer="'+primername+' '+primer2name+'"'
							print >> alltaboutput, 'FT                   /colour=3'
							print >> alltaboutput, 'FT                   /match_cost='+str(mincost2)
							if mincost2==0:
								print >> alltaboutput, 'FT                   /colour=3'
							else:
								print >> alltaboutput, 'FT                   /colour=1'
							print >> alltaboutput, 'FT   misc_feature    '+str(p1max+totallength+1)+".."+str(p2min+totallength)
							print >> alltaboutput, 'FT                   /sequence='+seq
							print >> alltaboutput, 'FT                   /product="'+primername+'"'
							print >> alltaboutput, 'FT                   /colour=2'
				
				
				
						elif options.circular:
							if (p2min+len(seqs[seq]))>=p1max and ((p2min+len(seqs[seq]))-p1max)<=patterns[primername][2] and p1[1]!=p2[1]:
								found=True
								if count==0:
									print >> primeroutput[primername], ">"+seq
								else:
									print >> primeroutput[primername], ">"+seq+"_"+str(count+1)
								count+=1
								if options.incprimer:
									print >> primeroutput[primername], seqs[seq][p1min:]+seqs[seq][:p2max]
								else:
									print >> primeroutput[primername], seqs[seq][p1max:]+seqs[seq][:p2min]
								print "  Found matching circularised region for "+primername+" on forward strand of contig "+seq+" from", p1max, "to", p2min, "with total cost of", mincost1+mincost2
								
								
								
								print >> alltaboutput, 'FT   misc_feature    '+str(p1min+totallength+1)+".."+str(p1max+totallength)
								print >> alltaboutput, 'FT                   /sequence='+seq
								print >> alltaboutput, 'FT                   /primer="'+primername+' '+primer1name+'"'
								print >> alltaboutput, 'FT                   /match_cost='+str(mincost1)
								if mincost1==0:
									print >> alltaboutput, 'FT                   /colour=3'
								else:
									print >> alltaboutput, 'FT                   /colour=1'
								print >> alltaboutput, 'FT   misc_feature    '+str(p2min+totallength+1)+".."+str(p2max+totallength)
								print >> alltaboutput, 'FT                   /sequence='+seq
								print >> alltaboutput, 'FT                   /primer="'+primername+' '+primer2name+'"'
								print >> alltaboutput, 'FT                   /colour=3'
								print >> alltaboutput, 'FT                   /match_cost='+str(mincost2)
								if mincost2==0:
									print >> alltaboutput, 'FT                   /colour=3'
								else:
									print >> alltaboutput, 'FT                   /colour=1'
								print >> alltaboutput, 'FT   misc_feature    join('+str(p1max+totallength+1)+".."+str(len(seqs[seq])+totallength)+","+str(1+totallength+1)+".."+str(p2min+totallength)+")"
								print >> alltaboutput, 'FT                   /sequence='+seq
								print >> alltaboutput, 'FT                   /product="'+primername+'"'
								print >> alltaboutput, 'FT                   /colour=2'
							
							
							elif (p1min+len(seqs[seq]))>=p2max and ((p1min+len(seqs[seq]))-p2max)<=patterns[primername][2] and ((p1min+len(seqs[seq]))-p2max)>=patterns[primername][3] and p1[1]!=p2[1]:
									found=True
									if count==0:
										print >> primeroutput[primername], ">"+seq
									else:
										print >> primeroutput[primername], ">"+seq+"_"+str(count+1)
									count+=1
									if options.incprimer:
										print >> primeroutput[primername], revcomp(seqs[seq][:p1max])+revcomp(seqs[seq][p2min:])
									else:
										print >> primeroutput[primername], revcomp(seqs[seq][:p1min])+revcomp(seqs[seq][p2max:])
										
									print "  Found matching circularised region for "+primername+" on reverse strand of contig "+seq+" from", p2max, "to", p1min, "with total cost of", mincost1+mincost2
									
									
									print >> alltaboutput, 'FT   misc_feature    complement('+str(p1min+totallength+1)+".."+str(p1max+totallength)+")"
									print >> alltaboutput, 'FT                   /sequence='+seq
									print >> alltaboutput, 'FT                   /primer="'+primername+' '+primer1name+'"'
									print >> alltaboutput, 'FT                   /match_cost='+str(mincost1)
									if mincost1==0:
										print >> alltaboutput, 'FT                   /colour=3'
									else:
										print >> alltaboutput, 'FT                   /colour=1'
									print >> alltaboutput, 'FT   misc_feature    complement('+str(p2min+totallength+1)+".."+str(p2max+totallength)+")"
									print >> alltaboutput, 'FT                   /sequence='+seq
									print >> alltaboutput, 'FT                   /primer="'+primername+' '+primer2name+'"'
									print >> alltaboutput, 'FT                   /match_cost='+str(mincost2)
									if mincost2==0:
										print >> alltaboutput, 'FT                   /colour=3'
									else:
										print >> alltaboutput, 'FT                   /colour=1'
									print >> alltaboutput, 'FT   misc_feature    join(complement('+str(p2min+totallength+1)+".."+str(len(seqs[seq])+totallength)+"),complement("+str(totallength+1)+".."+str(p1max+totallength)+")"
									print >> alltaboutput, 'FT                   /sequence='+seq
									print >> alltaboutput, 'FT                   /product="'+primername+'"'
									print >> alltaboutput, 'FT                   /colour=2'
							
							
						
				
						
				if not found:
					#print "  Found no pairs of matching primers within the specified product size range of "+str(patterns[primername][3])+" to "+str(patterns[primername][2])+" bases"
					if options.potentials:
						#print "  Printing best individual primer matches to tab file"
						for p1 in primer1matches:
							if p1[0][0]>p1[0][1]:
								print >> alltaboutput, 'FT   misc_feature    complement('+str(p1[0][1]+1+totallength)+".."+str(p1[0][0]+totallength)+")"
								print >> alltaboutput, 'FT                   /sequence='+seq
								print >> alltaboutput, 'FT                   /primer="Potential '+primername+' '+primer1name+'"'
								print >> alltaboutput, 'FT                   /match_cost='+str(mincost1)
								print >> alltaboutput, 'FT                   /colour=4'
							else:
								print >> alltaboutput, 'FT   misc_feature    '+str(p1[0][0]+1+totallength)+".."+str(p1[0][1]+totallength)
								print >> alltaboutput, 'FT                   /sequence='+seq
								print >> alltaboutput, 'FT                   /primer="Potential '+primername+' '+primer1name+'"'
								print >> alltaboutput, 'FT                   /match_cost='+str(mincost1)
								print >> alltaboutput, 'FT                   /colour=4'
						for p2 in primer2matches:
							if p2[0][0]>p2[0][1]:
								print >> alltaboutput, 'FT   misc_feature    complement('+str(p2[0][1]+1+totallength)+".."+str(p2[0][0]+totallength)+")"
								print >> alltaboutput, 'FT                   /sequence='+seq
								print >> alltaboutput, 'FT                   /primer="Potential '+primername+' '+primer2name+'"'
								print >> alltaboutput, 'FT                   /match_cost='+str(mincost2)
								print >> alltaboutput, 'FT                   /colour=4'
							else:
								print >> alltaboutput, 'FT   misc_feature    '+str(p2[0][0]+1+totallength)+".."+str(p2[0][1]+totallength)
								print >> alltaboutput, 'FT                   /sequence='+seq
								print >> alltaboutput, 'FT                   /primer="Potential '+primername+' '+primer2name+'"'
								print >> alltaboutput, 'FT                   /match_cost='+str(mincost2)
								print >> alltaboutput, 'FT                   /colour=5'
#				else:
#					foundprimers.append(primername)
				sys.stdout.flush()
#			taboutput.close()
			totallength+=len(seqs[seq])
			alltaboutput.flush()
		alltaboutput.close()
	
	for primername in primeroutput.keys():
		primeroutput[primername].close()
	
	
#	print foundprimers

