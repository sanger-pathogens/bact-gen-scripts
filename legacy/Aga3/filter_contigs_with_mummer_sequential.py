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

sys.path.extend(map(os.path.abspath, ['/nfs/pathogen/sh16_scripts/modules/']))
from Si_general import *
from Si_SeqIO import *



####################
# Set some globals #
####################


MUMMER_DIR=""
MIN_LENGTH=1000


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options] <directory containing novel sequence files>"
	parser = OptionParser(usage=usage)
	
	parser.add_option("-o", "--output", action="store", dest="output", help="Prefix for output files", default="", metavar="STRING")
	parser.add_option("-p", "--promer", action="store_true", dest="promer", help="Use promer instead of nucmer", default=False)
	parser.add_option("-t", "--tab", action="store_true", dest="tab", help="Create tab file showing aligned blocks", default=False)
	
	
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


	
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
	
	novelties=[]
	
	for arg in args:
		if os.path.isfile(arg):
			novelties.append(arg)
		if os.path.isdir(arg):
			for file in os.listdir(arg):
			    if file.endswith("novel.fasta"):
			        novelties.append(arg+"/"+file)
	
	
	chars = string.ascii_letters + string.digits
	tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))
	
	os.system("cp "+novelties[0]+" "+tmpname+".fa")
	
	print novelties
	
	for novnum, novelty in enumerate(novelties[1:]):
		
		try:
			refseqs=read_seq_file(tmpname+".fa")
		except StandardError:
			DoError("Cannot open reference file")
		
		reflist=[]
		reflens=[]
		lens={}
		seqs={}
		totlen=0
		count=0
		for seq in refseqs:
			newname=seq.id
	#		if seq.id in querylist:
				#DoError("You have two identically named sequences in your query fasta file")
			count=1
			while newname in reflist:
				newname=str(count)+seq.id
				count+=1
			seqs[newname]=str(seq.seq).upper()
			reflist.append(newname)
			reflens.append(totlen)
			
		if count==0:
			os.system("cp "+novelties[novnum]+" "+tmpname+".fa")
			continue
		
		try:
			queryseqs=read_seq_file(novelty)
		except StandardError:
			DoError("Cannot open query file")
		
		querylist=[]
		count=0
		for seq in queryseqs:
			newname=seq.id
	#		if seq.id in querylist:
				#DoError("You have two identically named sequences in your query fasta file")
			count=1
			while newname in querylist:
				newname=str(count)+seq.id
				print newname
				count+=1
			querylist.append(newname)
			seqs[newname]=str(seq.seq).upper()
			totlen+=len(seq.seq)
			lens[newname]=len(seq.seq)
		if count==0:
			continue
		
		if options.promer:
			promerargs=shlex.split(MUMMER_DIR+"promer  --maxmatch --nosimplify -p "+options.output+" "+tmpname+".fa "+novelty)
			returnval = subprocess.call(promerargs)
			if returnval!=0:
				DoError("Promer failed with return value "+str(returnval))
		else:
			#Run nucmer
			nucmerargs=shlex.split(MUMMER_DIR+"nucmer  --maxmatch --nosimplify -p "+options.output+" "+tmpname+".fa "+novelty)
			returnval = subprocess.call(nucmerargs)
			if returnval!=0:
				fastahandle=open(options.output+"_unmatched.fasta", "w")
				fastahandle.close()
				DoError("Nucmer failed with return value "+str(returnval))
		
		#Run delta-filter to filter repetitive regions
	#	deltafilterargs=shlex.split(MUMMER_DIR+"delta-filter -l 1000 -i 80 "+options.output+".delta")
	#	handle=open(options.output+".filter", "w")
	#	returnval = subprocess.call(deltafilterargs, stdout=handle)
	#	handle.close()
	#	if returnval!=0:
	#			DoError("delta-filter failed with return value "+str(returnval))
		
		blocks={}
		for refname in reflist:
			if not refname in blocks:
				blocks[refname]={}
		for queryname in querylist:
			if not queryname in blocks:
				blocks[queryname]={}
	
		
		for line in open(options.output+".delta"):
			line=line.strip()
			words=line.split()
			if len(words)==4 and line[0]==">":
				refname=words[0][1:]
				queryname=words[1]
			elif len(words)==7:
				if int(words[0])<int(words[1]):
					rstart=words[0]
					rend=words[1]
				else:
					rstart=words[1]
					rend=words[0]
				
				if int(words[2])<int(words[3]):
					qstart=words[2]
					qend=words[3]
				else:
					qstart=words[3]
					qend=words[2]
				if refname!=queryname:
					
					if not queryname in blocks[refname]:
						blocks[refname][queryname]=[]
					blocks[refname][queryname].append([int(rstart), int(rend), int(qstart), int(qend), "s"])
					if not refname in blocks[queryname]:
						blocks[queryname][refname]=[]
					blocks[queryname][refname].append([int(qstart), int(qend), int(rstart), int(rend), "r"])
				
		
		
		blockorder=[]
		for x, block in enumerate(blocks):
			blockorder.append([block[0],x])
		
		blockorder.sort()
		blockorder.reverse
	
	
		
		fastahandle=open(tmpname+".fa", "w")
		for subject in blocks:
			
			ssequence=seqs[subject]
			slength=len(ssequence)
			if slength<MIN_LENGTH:
				continue
			matched_regions=[0]*slength
			prev=0
			count=1
			s_isolate=subject.split("_50x")[0]
			old_q=""
			rtrim=slength
			ltrim=0
			keep=True
			old_rtrim=ltrim
			old_ltrim=rtrim	
			while old_rtrim!=rtrim or old_ltrim!=ltrim:
				old_rtrim=rtrim
				old_ltrim=ltrim	
				for query in blocks[subject]:
					
					if keep==False:
						break
					qsequence=seqs[query]
					qlength=len(qsequence)
		#			if qlength<MIN_LENGTH:
		#				continue
					q_isolate=query.split("_50x")[0]
					
					if old_q!="" and old_q!=query:
						
							
						max_run=0
						run=0
						if rtrim>ltrim:
							for m in matched_regions[ltrim:rtrim]:
								if m==0:
									run+=1
								else:
									if run>max_run:
										max_run=run
									run=0
							if run>max_run:
								max_run=run
		#				print subject, old_q, max_run,ltrim, rtrim
						if max_run<MIN_LENGTH:
							keep=False
					
						matched_regions=[0]*slength
					
					for region in blocks[subject][query]:
		#				slength=matched_regions.count(0)
						#print  subject, s_isolate, slength, region[0], region[1], query, q_isolate, qlength, region[2], region[3]
						if (qlength>slength or (qlength==slength and query>subject)):
							for x in xrange(region[0],region[1]):
								matched_regions[x]=1
							if region[0]<ltrim+MIN_LENGTH and region[1]>ltrim:
								ltrim=region[1]
							if region[1]>rtrim-MIN_LENGTH and region[0]<rtrim:
								rtrim=region[0]
					old_q=query
			
				#print subject, query, matched_regions
				if old_q!="":
					max_run=0
					run=0
					if rtrim>ltrim:
						for m in matched_regions[ltrim:rtrim]:
							if m==0:
								run+=1
							else:
								if run>max_run:
									max_run=run
								run=0
						if run>max_run:
							max_run=run
		#			print subject, old_q, max_run,ltrim, rtrim
					if max_run<MIN_LENGTH:
						keep=False
					
			if keep==False:
				print "Removing sequence", subject
			elif ltrim!=0 or rtrim!=slength:
				print "Trimming", ltrim, "bases from 5' and", slength-rtrim, "bases from 3' of", subject
				print >> fastahandle, ">"+subject
				print >> fastahandle, ssequence[ltrim:rtrim]
			else:
				print "keeping", subject
				print >> fastahandle, ">"+subject
				print >> fastahandle, ssequence
		fastahandle.close()
				
					
			
			
			
			
			
		
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
	
	
#	if options.tab:
#		tabhandle.close()
	os.system("mv "+tmpname+".fa "+options.output+"_filtered.fasta")
			
				
				
			
		
		
	
	
	
	
	
