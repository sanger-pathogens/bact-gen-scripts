#!/usr/bin/env python
# !/usr/bin/env python
import string
import os, sys
from Bio import Seq
from Bio import AlignIO
from Bio.Align import Generic
from Bio.Alphabet import IUPAC, Gapped
from optparse import OptionParser, OptionGroup
import glob
import shlex, subprocess

sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_general import *
from Si_SeqIO import *



####################
# Set some globals #
####################


MUMMER_DIR=""


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	
	parser.add_option("-r", "--reference", action="store", dest="ref", help="Reference fasta file name", default="", metavar="FILE")
	parser.add_option("-q", "--query", action="store", dest="query", help="Query (contigs) fasta file name", default="", metavar="FILE")
	parser.add_option("-o", "--output", action="store", dest="output", help="Prefix for output files", default="", metavar="STRING")
	parser.add_option("-p", "--promer", action="store_true", dest="promer", help="Use promer instead of nucmer", default=False)
	parser.add_option("-t", "--tab", action="store_true", dest="tab", help="Create tab file showing aligned blocks", default=False)
	
	
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


#	if options.processors>cpu_count():
#		DoError('You do not have '+str(options.processors)+' processors')
#
	if options.ref=='':
		DoError('No reference fasta file selected (-r)')
	elif not os.path.isfile(options.ref):
		DoError('Cannot find file '+options.ref)
	
	if options.query=='':
		DoError('No query fasta file selected (-q)')
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
	
	
	try:
		refseqs=read_seq_file(options.ref)
	except StandardError:
		DoError("Cannot open reference file")
	
	reflist=[]
	reflens=[]
	lens={}
	seqs={}
	totlen=0
	for seq in refseqs:
		if seq.id in reflist:
			DoError("You have two identically named sequences in your reference fasta file")
		reflist.append(seq.id)
		reflens.append(totlen)
		
	
	
	try:
		queryseqs=read_seq_file(options.query)
	except StandardError:
		DoError("Cannot open query file")
	
	querylist=[]
	for seq in queryseqs:
		if seq.id in querylist:
			DoError("You have two identically named sequences in your query fasta file")
		querylist.append(seq.id)
		seqs[seq.id]=str(seq.seq).upper()
		totlen+=len(seq.seq)
		lens[seq.id]=len(seq.seq)
	
	
	if options.promer:
		promerargs=shlex.split(MUMMER_DIR+"promer  --maxmatch --nosimplify -p "+options.output+" "+options.ref+" "+options.query)
		returnval = subprocess.call(promerargs)
		if returnval!=0:
			DoError("Promer failed with return value "+str(returnval))
	else:
		#Run nucmer
		nucmerargs=shlex.split(MUMMER_DIR+"nucmer  --maxmatch --nosimplify -p "+options.output+" "+options.ref+" "+options.query)
		returnval = subprocess.call(nucmerargs)
		if returnval!=0:
			fastahandle=open(options.output+"_unmatched.fasta", "w")
			fastahandle.close()
			DoError("Nucmer failed with return value "+str(returnval))
	
	#Run delta-filter to filter repetitive regions
	deltafilterargs=shlex.split(MUMMER_DIR+"delta-filter -l 1000 -i 98 "+options.output+".delta")
	handle=open(options.output+".filter", "w")
	returnval = subprocess.call(deltafilterargs, stdout=handle)
	handle.close()
	if returnval!=0:
			DoError("delta-filter failed with return value "+str(returnval))
	
	blocks=[]
	for x, refname in enumerate(reflist):
		print "Aligning query sequences against", refname
		for queryname in querylist:
			if refname==queryname:
				continue
			showalignsarg = shlex.split(MUMMER_DIR+"show-aligns -r "+options.output+".filter "+refname+" "+queryname)
			returnval = subprocess.Popen(showalignsarg, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			align = returnval.stdout.readlines()
			inblock=False
			
			for line in align:
				
				words=line.replace("^","").strip().split()
				
				
				if inblock and len(words)>1:
					if words[1]=="END":
						inblock=False
						continue
					blockline+=1
					if blockline%2==0:
						blocks[-1]["queryseq"]=blocks[-1]["queryseq"]+words[1].replace(".","-")
					else:
						blocks[-1]["refseq"]=blocks[-1]["refseq"]+words[1].replace(".","-")
						
				elif len(words)>1 and words[1]=="BEGIN":
					inblock=True
					blockline=0
					blocks.append({"start":int(words[5])+reflens[x], "ref_start": words[5], "ref_end": words[7],"end":int(words[7])+reflens[x], "refseq":"", "queryseq":"", "query":queryname, "ref":refname})
			
	
	
	blockorder=[]
	for x, block in enumerate(blocks):
		blockorder.append([block["start"],x])
	
	blockorder.sort()
	blockorder.reverse
	filtered={}
	for block in blockorder:
		filtered[blocks[block[1]]["ref"]]={}
	
	
	
	if options.tab:
	
		tabhandle=open(options.output+".tab", "w")
		print >> tabhandle, "ID   mummer_alignment"
		
		
	for block in blockorder:
		percent_ID=blocks[block[1]]["end"]-blocks[block[1]]["start"]
	
		currbase=blocks[block[1]]["start"]
		in_insertion=False
		insertion_start=-1
		
		for x, base in enumerate(blocks[block[1]]["refseq"]):
			if base!="-":
				currbase+=1
				if in_insertion:
					in_insertion=False
					if options.tab:
						print >> tabhandle, "FT   insertion       "+str(insertion_start-1)+".."+str(insertion_start)
						print >> tabhandle, "FT                   /query="+blocks[block[1]]["query"]
						print >> tabhandle, "FT                   /reference="+blocks[block[1]]["ref"]
						print >> tabhandle, "FT                   /reference_start="+blocks[block[1]]["ref_start"]
						print >> tabhandle, "FT                   /reference_end="+blocks[block[1]]["ref_end"]
						#print >> tabhandle, "FT                   /refseq="+blocks[block[1]]["refseq"][start_x-5:x+5]
						print >> tabhandle, "FT                   /insertion="+blocks[block[1]]["queryseq"][start_x:x]
						print >> tabhandle, "FT                   /colour=2"
			elif not in_insertion:
				in_insertion=True
				insertion_start=currbase
				start_x=x
		
		currbase=blocks[block[1]]["start"]
		in_deletion=False
		deletion_start=-1
		for x, base in enumerate(blocks[block[1]]["queryseq"]):
			if blocks[block[1]]["refseq"][x]!="-":
				currbase+=1
			if base!="-":
				if in_deletion:
					in_deletion=False
					if options.tab:
						print >> tabhandle, "FT   deletion        "+str(deletion_start)+".."+str(currbase-1)
						print >> tabhandle, "FT                   /query="+blocks[block[1]]["query"]
						print >> tabhandle, "FT                   /reference="+blocks[block[1]]["ref"]
						print >> tabhandle, "FT                   /deletion="+blocks[block[1]]["refseq"][start_x:x]
						print >> tabhandle, "FT                   /reference_start="+blocks[block[1]]["ref_start"]
						print >> tabhandle, "FT                   /reference_end="+blocks[block[1]]["ref_end"]
						print >> tabhandle, "FT                   /queryseq="+blocks[block[1]]["queryseq"][start_x-5:x+5]
						print >> tabhandle, "FT                   /colour=5"
					percent_ID-=((currbase-1)-deletion_start)
				if base!=blocks[block[1]]["refseq"][x] and blocks[block[1]]["refseq"][x]!="-":
					if options.tab:
						print >> tabhandle, "FT   SNP             "+str(currbase-1)
						print >> tabhandle, "FT                   /query="+blocks[block[1]]["query"]
						print >> tabhandle, "FT                   /reference="+blocks[block[1]]["ref"]
						print >> tabhandle, "FT                   /reference_start="+blocks[block[1]]["ref_start"]
						print >> tabhandle, "FT                   /reference_end="+blocks[block[1]]["ref_end"]
						print >> tabhandle, "FT                   /querybase="+base
						print >> tabhandle, "FT                   /refseq="+blocks[block[1]]["refseq"][x]
						print >> tabhandle, "FT                   /colour=4"
					percent_ID-=1
			elif not in_deletion:
				in_deletion=True
				deletion_start=currbase
				start_x=x
		
		final_ID=((float(percent_ID)/(blocks[block[1]]["end"]-blocks[block[1]]["start"]))*100)
		if options.tab:
			print >> tabhandle, "FT   misc_feature    "+str(blocks[block[1]]["start"])+".."+str(blocks[block[1]]["end"])
			print >> tabhandle, "FT                   /query="+blocks[block[1]]["query"]
			print >> tabhandle, "FT                   /reference="+blocks[block[1]]["ref"]
			print >> tabhandle, "FT                   /reference_start="+blocks[block[1]]["ref_start"]
			print >> tabhandle, "FT                   /reference_end="+blocks[block[1]]["ref_end"]
			print >> tabhandle, "FT                   /ID="+str(final_ID)
		if final_ID>98 and (blocks[block[1]]["end"]-blocks[block[1]]["start"])>1000:
			if options.tab:
				print >> tabhandle, "FT                   /colour=2"
			if not blocks[block[1]]["query"] in filtered[blocks[block[1]]["ref"]]:
				filtered[blocks[block[1]]["ref"]][blocks[block[1]]["query"]]=[]
			filtered[blocks[block[1]]["ref"]][blocks[block[1]]["query"]].append([int(blocks[block[1]]["ref_start"]), int(blocks[block[1]]["ref_end"])])

		
		
	fastahandle=open(options.output+"_unmatched.fasta", "w")
	for subject in filtered:
		ssequence=seqs[subject]
		slength_orig=len(ssequence)
		if slength_orig<1000:
			continue
		matched_regions=[0]*slength_orig
		prev=0
		count=1
		for query in filtered[subject]:
			qsequence=seqs[query]
			qlength=len(qsequence)
			if qlength<1000:
				continue
			for region in filtered[subject][query]:
				slength=matched_regions.count(0)
				print region, subject, slength, query, qlength
				if (region[0]<200 or (slength-region[1])<200) and (qlength>slength or (qlength==slength and query>subject)):
					for x in xrange(region[0],region[1]):
						matched_regions[x]=1
				
			print subject, query, matched_regions
				
		
		
		
		
		
		
#			for block in filtered[name]:
#				if block[0]-prev>1000:
#					#unmatched_regions.append([prev,block[0], sequence[prev:block[0]]])
#					print >> fastahandle, ">"+name+"#"+str(count)
#					print >> fastahandle, sequence[prev:block[0]]
#					count+=1
#				prev=block[1]
#			if length-prev>1000:
#				#unmatched_regions.append([prev,length, sequence[prev:]])
#				print >> fastahandle, ">"+name+"#"+str(count)
#				print >> fastahandle, sequence[prev:]
#		else:
#			print >> fastahandle, ">"+name+"#"+str(count)
#			print >> fastahandle, sequence
	
	
	if options.tab:
		tabhandle.close()
	fastahandle.close()
			
				
				
			
		
		
	
	
	
	
	
