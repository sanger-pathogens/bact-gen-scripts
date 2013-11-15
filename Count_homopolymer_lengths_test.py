#!/usr/bin/env python

##################
# Import modules #
##################

import os, sys, string
from optparse import OptionParser
from random import randrange, randint, choice
import pysam

##########################
# Error message function #
##########################

def DoError(errorstring):
	print "\nError:", errorstring
	print "\nFor help use -h or --help\n"
	sys.exit()


##############################
# Get command line arguments #
##############################

def get_user_options():
	usage = "usage: %prog [options] <bam files>"
	version="%prog 1.0. Written by Simon Harris, Wellcome Trust Sanger Institute, 2012"
	parser = OptionParser(usage=usage, version=version)
	
	#do not allow arguments to be interspersed. e.g. -a -b arg1 agr2. MUST be -a arg1 -b arg2.
	parser.disable_interspersed_args()
	
	#No need for reference to be loaded
	parser.add_option("-r", "--reference", action="store", dest="ref", help="Reference fasta file", default="", metavar="FILE")
	parser.add_option("-m", "--min", action="store", dest="min", help="minimum length of homopolymer to analyse", default=6, type="int")
	parser.add_option("-o", "--output", action="store", dest="output", help="output file name prefix", default="", metavar="FILE")
#	parser.add_option("-f", "--features", action="store", dest="features", help="feature keys to use as homopolymer regions [default=%default]", default="misc_feature", metavar="FILE")


	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):

	if options.output=='':
		DoError("No output file selected")

#	elif options.ref=='':
#		DoError("No reference file selected")
#	elif not os.path.isfile(options.ref):	
#		DoError(options.ref+" is not a file")
	
	elif options.min<2:
		DoError("Minimum Homopolymer length must be greater than 1")
	
	if len(args)==0:
		DoError("no input bam files selected")
	
	for arg in args:
		if not os.path.isfile(arg):	
			DoError(arg+" is not a file")


	
	return


########
# Main #
########


if __name__ == "__main__":



	(options, args)=get_user_options()
	#print options, args
	check_input_validity(options, args)
	
	


	refseqs={}
	homopolymers={}
	regions={}
	region_lengths={}
	if sys.argv[1]=="-r":
		doerrors=True
		reflines=open(sys.argv[2]).read().split(">")[1:]
		for ref in reflines:
			seqlines=ref.split("\n")
			refname=seqlines[0].split()[0]
			refseqs[refname]=''.join(seqlines[1:]).upper()
			homopolymers[refname]=[]
			regions[refname]={}
			region_lengths[refname]={}
	
	for ref in refseqs:
		seq=refseqs[ref]
		runbase='N'
		start=0
		count=0
		for x in xrange(len(seq)):
			if seq[x]=='N' or seq[x]!=runbase:
				if x-start>=options.min and runbase!='N':
#					print start+1, x+1, x-start, seq[start:x]
					regions[ref][start]=[start,x]
					count+=1
				start=x
				runbase=seq[x]
		if x-start>=options.min and runbase!='N':
#			print start+1, x+1, x-start, seq[start:x]
			regions[ref][start]=[start,x]
			count+=1
				
	print "Found", count, "homopolymers greater than or equal to ", options.min, "bases long"
#	sys.exit()
	
	
#	if options.features!="":
#		features=options.features.split(",")
#	else:
#		features=["misc_feature"]
#	
#	
#	regions={}
#	region_names={}
#	region_lengths={}
#	
#	#read the tab file
#	print "Reading tab file"
#	sys.stdout.flush()
#	for line in open(options.tab, "rU"):
#		line =line.strip()
#		words=line.split()
#		if len(words)==3 and words[0]=="FT" and words[1] in features:
#			location=words[2]
#			words[2]=words[2].replace("complement(","")
#			words[2]=words[2].replace(")","")
#			try:
#				startpos=int(words[2].split("..")[0])
#				endpos=int(words[2].split("..")[1])
#			except StandardError:
#				DoError("This script currently only works for feature locations of the form xxxx..xxxx")
#				
#			if startpos>endpos:
#				regions[location]=[endpos, startpos]
#			else:
#				regions[location]=[startpos, endpos]
#			
#			region_lengths[location]={}
#		elif len(words)==2 and words[1].split("=")[0]=="/label":
#			region_names[location]=words[1].split("=")[1].replace('"','')

#	for region in regions:
#		try:
#			print regions[region], region_names[region], region_lengths[region]
#		except StandardError:
#			print regions[region], region, region_lengths[region]
	
	print "Finding regions in bam files"
	sys.stdout.flush()
	for filename in args:
		print "\t"+filename+"..."
		sys.stdout.flush()
		if filename.split(".")[-1]=="bam":
			samfile = pysam.Samfile( filename, "rb" )
		elif filename.split(".")[-1]=="sam":
			samfile = pysam.Samfile( filename, "r" )
		else:
			print filename, "not a readable bam file"
			continue
		
		refs=samfile.references
		lengths=samfile.lengths
		
		if len(refs)!=len(refseqs):
			DoError("bam has different number of reference sequences to reference fasta file")
		else:
			for ref in refs:
				if not ref in refseqs:
					DoError("bam and reference fasta file do not match")
		
		ref=refs[0]
		reflen=lengths[0]
		
		for ref in refs:
			
			for region in regions[ref]:
				
				if not region in region_lengths[ref]:
					region_lengths[ref][region]={}
				
				if not filename in region_lengths[ref][region]:
					region_lengths[ref][region][filename]={"A": {}, "G": {},"C": {}, "T": {}}
			
				if regions[ref][region][0]<1 or regions[ref][region][1]>reflen:
					DoError("At least one of your regions is out of the range of your reference")
					
				iter = samfile.fetch( ref, regions[ref][region][0]-2, regions[ref][region][0]-1)
				
				read_dict={}
				
				start_reads=set([])
				for read in iter:
					if read.is_read1:
						fr="f"
					else:
						fr="r"
					start_reads.add(read.qname+"fr")
					read_dict[read.qname+"fr"]=read
				
				
				iter = samfile.fetch( ref, regions[ref][region][1]+1, regions[ref][region][1]+2)
				
				end_reads=set([])
				for read in iter:
					if read.is_read1:
						fr="f"
					else:
						fr="r"
					end_reads.add(read.qname+"fr")
					read_dict[read.qname+"fr"]=read
				
				
				intersection_reads=start_reads.intersection(end_reads)
				
				for readname in intersection_reads:
					read=read_dict[readname]
					
					refpos=read.pos
					readpos=0
					readseq=read.seq
					
					repbit=""
					repstart=-1
					repend=-1
					
					for cig in read.cigar:
						if cig[0]==0:
							for x in range(0,cig[1]):
								if refpos==regions[ref][region][1]:
									repend=readpos
									break
								if refpos==regions[ref][region][0]-1:
									repstart=readpos
								readpos+=1
								refpos+=1
							
						elif cig[0]==1:
							for x in range(0,cig[1]):
								readpos+=1
						elif cig[0]==2:
							for x in range(0,cig[1]):
								if refpos==regions[ref][region][1]:
									repend=readpos
									break
								refpos+=1
						elif cig[0]==4:
							for x in range(0,cig[1]):
								readpos+=1
						elif cig[0]==5:
							continue
						else:
							print cig
						if refpos==regions[ref][region][0]-1:
							repstart=readpos
						if refpos==regions[ref][region][1]:
							repend=readpos
							break
						
					
					if repstart==-1 or repend==-1:
						continue
					
					
					basecount={"A":0, "C":0, "G":0, "T":0, 'N':0}
					
					for base in read.seq[repstart:repend]:
						basecount[base]+=1
					
					
					maxbase=0
					foundbases=[]
					for base in basecount:
						if basecount[base]>maxbase:
							maxbase=basecount[base]
							repbase=base
						if basecount[base]>0 and base in foundbases:
							foundbases.append(base)
						if len(foundbases)>1:
							break
					
					
					if repbase not in ["A", "T", "C", "G"]:
						continue
					
					if len(foundbases)>1:
						continue
									
	
					while read.seq[repstart]==repbase:
						repstart-=1
					repstart+=1
					
					
					while read.seq[repend]==repbase:
						repend+=1
						if repend>=read.rlen:
							break
					if repend>=read.rlen:
						break
					repend-=1
					
					
					runlen=0
					maxrunlen=0
					for x in range(repstart,repend+1):
						if read.seq[x]==repbase:
							runlen+=1
						else:
							if runlen>maxrunlen:
								maxrunlen=runlen
							runlen=0
							
					if runlen<3:
						continue
					
					
					
					#print refpos, readpos, read.seq[readpos], read.seq[readpos-replen-1:readpos+1], read.is_reverse, replen
					repbase=repbase.upper()
					if not repbase in region_lengths[ref][region][filename]:
						DoError("Found illegal homopolymer base: "+repbase)
					
					if not runlen in region_lengths[ref][region][filename][repbase]:
						region_lengths[ref][region][filename][repbase][runlen]=0
					
					region_lengths[ref][region][filename][repbase][runlen]+=1
					
				#print read.pos, read.pos+read.alen, regions[region], read.seq[regions[region][0]-(read.pos+2):(regions[region][0]-read.pos)+10], read.is_reverse, regions[region][0]-read.pos
			
#		print region_lengths[ref]
#	sys.exit()
	
	
	print "Printing output files"
	sys.stdout.flush()
	
	results={}
	
	for arg in args:
		results[arg]=[]
	headings=[]
	for ref in refseqs:
		for region in region_lengths[ref]:
			
			variable=False
			
			#print region_lengths[ref][region]
			bases=[]
			basecounts={}
			for filename in region_lengths[ref][region]:
				for base in region_lengths[ref][region][filename]:
					if len(region_lengths[ref][region][filename][base])>0 and base not in bases:
						bases.append(base)
						basecounts[base]=[]
#					if len(region_lengths[ref][region][filename][base])>1:
#						variable=True
					maxbase=0
					maxbasesup=0
					for basecount in region_lengths[ref][region][filename][base]:
#						print "here", region_lengths[ref][region][filename][base][basecount], basecount
						if region_lengths[ref][region][filename][base][basecount]>maxbasesup and region_lengths[ref][region][filename][base][basecount]>10:#NOTE: this 10 means it will only count a maximum value that's greater than 10 to decide if things are different than the reference
							maxbase=basecount
							maxbasesup=region_lengths[ref][region][filename][base][basecount]
						if not basecount in basecounts[base]:
							basecounts[base].append(basecount)
					if maxbase!=0 and maxbase!=regions[ref][region][1]-regions[ref][region][0]:
#						print "Found difference", filename, region_lengths[ref][region][filename][base], maxbase, regions[ref][region][1]-regions[ref][region][0]
						variable=True
						
#			print basecounts	
			if not variable:
				print "Skipping", region, "as the most common length does not vary from the reference"
				continue
			
			
			for base in bases:
				basecounts[base].sort()
#				if ref in region_names and region in region_names[ref]:
#					output=open(options.output+"_"+region_names[ref][region]+"_"+base+".csv", "w")
#				else:
#					output=open(options.output+"_"+region+"_"+base+".csv", "w")
				output=open(options.output+"_"+str(region)+"_"+base+".csv", "w")
				headings.append(options.output+"_"+str(region)+"_"+base)
				print >> output, ",".join([str(region)]+map(str,basecounts[base])+["Max_length"])
				
				for filename in args:
				
					filebasecount=[]
					filemax=0
					filemaxlen=0
					for basecount in basecounts[base]:
						if basecount in region_lengths[ref][region][filename][base]:
							filebasecount.append(region_lengths[ref][region][filename][base][basecount])
							if region_lengths[ref][region][filename][base][basecount]>filemax:
								filemax=region_lengths[ref][region][filename][base][basecount]
								filemaxlen=basecount
						else:
							filebasecount.append(0)
				
					print >> output, ",".join([filename.split("/")[-1]]+map(str,filebasecount)+[str(filemaxlen)])
					results[filename].append(str(filemaxlen))
				output.close()
	output=open("test.csv","w")
	
	print >> output, ','.join(["File"]+headings)
	for filename in results:
		print >> output, ','.join([filename.split("/")[-1].rstrip(".bam")]+results[filename])
	output.close()
	print "Done"
	sys.stdout.flush()	
		
		
