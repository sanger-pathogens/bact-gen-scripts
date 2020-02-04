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


#############################################
# Function to reverse complement a sequence #
#############################################

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
	parser.add_option("-s", "--supress", action="store_true", dest="supress", help="Supress creation of fasta of novel regions", default=False)
	
	
	
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
	refstarts={}
	refends={}
	lens={}
	seqs={}
	totlen=0
	for seq in refseqs:
		if seq.id in reflist:
			DoError("You have two identically named sequences in your reference fasta file")
		reflist.append(seq.id)
		reflens.append(totlen)
		refstarts[seq.id]=totlen
		refends[seq.id]=totlen+len(seq.seq)
		seqs[seq.id]=str(seq.seq).upper()
		totlen+=len(seq.seq)
		lens[seq.id]=len(seq.seq)

		try:
			refstarts[seq.id]=int(seq.description.split("#")[1].split("..")[0].strip())+1
		except ValueError:
			refstarts[seq.id]=float("Inf")
		except IndexError:
			refstarts[seq.id]=float("Inf")

	
	
	try:
		queryseqs=read_seq_file(options.query)
	except StandardError:
		DoError("Cannot open query file")
	
	querylist=[]
	qstarts={}
	qlengths={}
	qseqs={}
	totlen=0
	name_conversion={}
	qname_prefix=options.query.split("/")[-1].split(".")[0]+"$"
	for seq in queryseqs:
		if seq.id in querylist:
			DoError("You have two identically named sequences in your query fasta file")
		qname=qname_prefix+seq.id
		count=1
		while qname in querylist:
			qname=qname_prefix+seq.id+"_"+str(count)
			count+=1
		name_conversion[seq.id]=qname
		querylist.append(qname)
		qseqs[qname]=str(seq.seq).upper()
		qstarts[qname]=totlen
		totlen+=len(seq.seq)
		qlengths[qname]=len(seq.seq)
	
	if options.promer:
		promerargs=shlex.split(MUMMER_DIR+"promer -p "+options.output+" "+options.ref+" "+options.query)
		returnval = subprocess.call(promerargs)
		if returnval!=0:
			DoError("Promer failed with return value "+str(returnval))
	else:
		#Run nucmer
		nucmerargs=shlex.split(MUMMER_DIR+"nucmer -p "+options.output+" "+options.ref+" "+options.query)
		returnval = subprocess.call(nucmerargs)
		if returnval!=0:
#			fastahandle=open(options.output+"_unmatched.fasta", "w")
#			fastahandle.close()
			DoError("Nucmer failed with return value "+str(returnval))
	
	#Run delta-filter to filter repetitive regions
	deltafilterargs=shlex.split(MUMMER_DIR+"delta-filter -i 98 -l 1000 -r -q "+options.output+".delta")
	handle=open(options.output+".filter", "w")
	returnval = subprocess.call(deltafilterargs, stdout=handle)
	handle.close()
	if returnval!=0:
			DoError("delta-filter failed with return value "+str(returnval))
	
	
	blocks={}
	for qname in querylist:
		if not qname in blocks:
			blocks[qname]=[]

	
	for line in open(options.output+".delta"):
		line=line.strip()
		words=line.split()
		if len(words)==4 and line[0]==">":
			refname=words[0][1:]
			queryname=name_conversion[words[1]]
		elif len(words)==7:
			rstrand="f"
			if int(words[0])<int(words[1]):
				rstart=words[0]
				rend=words[1]
			else:
				rstart=words[1]
				rend=words[0]
				rstrand="r"
			
			if int(words[2])<int(words[3]):
				qstart=words[2]
				qend=words[3]
			else:
				qstart=words[3]
				qend=words[2]
				if rstrand=="r":
					rstrand="f"
				else:
					rstrand="r"
					
			if refname!=queryname:
				
				blocks[queryname].append([int(rstart), int(rend), int(qstart), int(qend), 100*((float(int(rend)-int(rstart)))-int(words[4]))/(float(int(rend)-int(rstart))), rstrand, refname])
	
	
	
	if options.tab:
	
		tabhandle=open(options.output+".tab", "w")
		print >> tabhandle, "ID   mummer_alignment"
		
		
	for query in blocks:
		blockorder=[]
		for x, block in enumerate(blocks[query]):
			blockorder.append([block[2],x])
		blockorder.sort()
		blockorder.reverse()
	
		for x in blockorder:
			block=blocks[query][x[1]]
			#print query, block
			rstart=block[0]
			rend=block[1]
			qstart=block[2]+qstarts[query]
			qend=block[3]+qstarts[query]
			percent_ID=float(block[4])
			strand=block[5]
			ref=block[6]
			
			if options.tab:
				print >> tabhandle, "FT   misc_feature    "+str(qstart)+".."+str(qend)
				print >> tabhandle, "FT                   /query="+query
				print >> tabhandle, "FT                   /reference="+ref
				print >> tabhandle, "FT                   /reference_start="+str(rstart)
				print >> tabhandle, "FT                   /reference_end="+str(rend)
				if block[5]=="f":
					print >> tabhandle, "FT                   /reference_strand=+"
				else:
					print >> tabhandle, "FT                   /reference_strand=-"
					
				print >> tabhandle, "FT                   /ID="+str(percent_ID)
				print >> tabhandle, 'FT                   /colour='+str(int(((percent_ID-50)/50)*255)), 0, str(int(255-(((percent_ID-50)/50)*255)))
	
	if options.tab:
		tabhandle.close()
	if options.supress:
		sys.exit()
	fastahandle=open(options.output+"_novel.fasta", "w")
	for query in qseqs:
		sequence=qseqs[query]
		length=len(sequence)
		if length<1000:
			continue
		#unmatched_regions=[]
		
		ordered_blocks={}
		
		count=1
		if query in blocks:
			prevrefblockend=float("Inf")
			blockorder=[]
			for x, block in enumerate(blocks[query]):
				blockorder.append([block[2],x])
			blockorder.sort()
#			blockorder.reverse()
			prev=0
			ordered_blocks[query]=[]
			for y in blockorder:
				ordered_blocks[query].append(blocks[query][y[1]])
			for x in xrange(len(ordered_blocks[query])):
				block=ordered_blocks[query][x]
				ref=block[6]
#				print query, ref,block
#				print query, ref, block[2],block[3], prev
				if block[3]<prev:
					continue
				if block[2]-prev>1000:
					
					#unmatched_regions.append([prev,block[0], sequence[prev:block[0]]])
					loc_start=""
					if prevrefblockend!=float("Inf"):
						loc_start=str(prevrefblockend)
#						else:
#							loc_string=loc_string+str(refstarts[ref])
					loc_end=""
#					if len(ordered_blocks[query])>x+1:
					if block[5]=="f":
						loc_end=str(refstarts[ref]+block[0])
					else:
						loc_end=str(refstarts[ref]+block[1])
					if block[5]=="f":
						print >> fastahandle, ">"+query.replace("#","_")+"_"+str(count)+"#"+loc_start+".."+loc_end
						print >> fastahandle, sequence[prev:block[2]]
						count+=1
#						print sequence[prev:block[2]]
#						print len(sequence), prev, block[2]
					else:
						if x>0 and ordered_blocks[query][x-1][5]=="f":
							print >> fastahandle, ">"+query.replace("#","_")+"_"+str(count)+"#"+loc_start+".."+loc_end
							print >> fastahandle, sequence[prev:block[2]]
							count+=1
#							print sequence[prev:block[2]]
#							print len(sequence), prev, block[2]
						else:
							print >> fastahandle, ">"+query.replace("#","_")+"_"+str(count)+"#"+loc_end+".."+loc_start
							print >> fastahandle, revcomp(sequence[prev:block[2]])
							count+=1
#							print sequence[prev:block[2]]
#							print len(sequence), prev, block[2]
						
					
				prev=block[3]
				if block[5]=="f":
					prevrefblockend=refstarts[ref]+block[1]
				else:
					prevrefblockend=refstarts[ref]+block[0]
				#print ref, query, prevrefblockend, blocks[query][x]
			#print length, prev
			if length-prev>1000:
				#unmatched_regions.append([prev,length, sequence[prev:]])
				loc_start=""
				if prevrefblockend!=float("Inf"):
					loc_start=str(prevrefblockend)
#						else:
#							loc_string=loc_string+str(refstarts[ref])
				loc_end=""
#					if len(ordered_blocks[query])>x+1:
				
				if len(ordered_blocks[query])==0 or block[5]=="f":
					print >> fastahandle, ">"+query.replace("#","_")+"_"+str(count)+"#"+loc_start+".."+loc_end
					print >> fastahandle, sequence[prev:]
					count+=1
				else:
					print >> fastahandle, ">"+query.replace("#","_")+"_"+str(count)+"#"+loc_end+".."+loc_start
					print >> fastahandle, revcomp(sequence[prev:])
					count+=1
		else:
			print >> fastahandle, ">"+query.replace("#","_")+"_"+str(count)+"#.."
			print >> fastahandle, sequence
			count+=1
	
	
	
	fastahandle.close()
			
				
				
			
		
		
	
	
	
	
	
