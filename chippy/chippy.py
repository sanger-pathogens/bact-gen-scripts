#!/usr/bin/env python

import os, sys
import random
from numpy import log
from math import floor
from mash_bounds import bounds

MASH_DIR="/nfs/users/nfs_s/sh16/mash-Linux64-v1.1.1/"
ref_sketch_size=10000
contig_sketch_size=10000
max_fragments=100
kmer=16
expected_identity=0.95
expected_length_shared=0.95
circular=True
allow_self=True

def sketch_file(filename, s=10000, k=16, i=True):
	if i:
		os.system(MASH_DIR+"mash sketch -i -k "+str(k)+" -s "+str(s)+" "+filename)
	else:
		os.system(MASH_DIR+"mash sketch -k "+str(k)+" -s "+str(s)+" "+filename)


def split_references(reffile):
	
	def calculate_sketch(seqlen):
		if seqlen>10000000:
			s=100000
		elif seqlen>1000000:
			s=10000
		elif seqlen>100000:
			s=1000
		elif seqlen>10000:
			s=100
		else:
			s=10
		return s
	
	ref_sequences={}
	sketches={}
	ref_lengths={}
	seq=[]
	name=""
	description=""
	for line in open(references, "rU"):
	        line=line.strip()
	        if len(line)>0:
	                if line[0]==">":
	                        if name!="":
					sketch_size=calculate_sketch(len(''.join(seq)))
					if not sketch_size in sketches:
						sketches[sketch_size]=[]
					sketches[sketch_size].append([name,description,''.join(seq)])
					ref_lengths[name]=len(''.join(seq))
	                        name=line.split()[0][1:]
	                        seq=[]
	                        if len(line.split())>1:
	                                description=" ".join(line.split()[1:])
	                else:
	                        seq.append(line.replace(" ",""))
	if name!="":
	        sketch_size=calculate_sketch(len(''.join(seq)))
		if not sketch_size in sketches:
			sketches[sketch_size]=[]
		sketches[sketch_size].append([name,description,''.join(seq)])
		ref_lengths[name]=len(''.join(seq))
	
	for sketch in sketches:
		filename=reffile+"."+str(sketch)+".fasta"
		output=open(filename, "w")
		for seq in sketches[sketch]:
			print >> output, ">"+seq[0]+" "+seq[1]
			print >> output, seq[2]
		output.close()
		#sketch_file(filename, s=sketch, k=kmer, i=True)
	
	return ref_lengths, sketches.keys()
	


try:
	references=sys.argv[1]
	assembly=sys.argv[2]
except:
	print "chippy.py <reference fasta> <contigs fasta>"
	sys.exit()
	

ref_lengths,sketches=split_references(references)
		
print len(ref_lengths), "contigs split into", len(sketches), "reference sketches"


for x, sketch in enumerate(sketches):
	sketch_file(assembly, s=sketch, k=kmer, i=True)
	
	if x==0:
		contigs={}
		os.system(MASH_DIR+"mash info -t "+assembly+".msh > "+assembly+".msh.info")
		for line in open(assembly+".msh.info", "rU"):
			line=line.strip()
			words=line.split()
			if len(words)>0:
				try:
					hashes=int(words[0])
					length=int(words[1])
					contigs[words[2]]={"length":length,"hashes":hashes,"distance":float("Inf"),"p":1.0,"hit":"", "hash_matches":"None"}
				except:
					continue
				
		print len(contigs), "contigs in sketch file of contigs"


	os.system(MASH_DIR+"mash dist "+references+"."+str(sketch)+".fasta.msh "+assembly+".msh > "+assembly+"."+str(sketch)+".dist.tab")
	os.system("rm -f "+assembly+".msh")
	for line in open(assembly+"."+str(sketch)+".dist.tab", "rU"):
		line=line.strip()
		words=line.split()
		if float(words[2])<contigs[words[1]]["distance"]:
			distance=float(words[2])
			p=float(words[3])
			hashes=words[4]
			
			contigs[words[1]]["distance"]=distance
			contigs[words[1]]["p"]=p
			contigs[words[1]]["hit"]=words[0]
			contigs[words[1]]["ref_length"]=ref_lengths[words[0]]
			contigs[words[1]]["hash_matches"]=hashes
			if float(ref_lengths[words[0]])>float(contigs[words[1]]["length"]):
				j_exp=(float(contigs[words[1]]["length"])/float(ref_lengths[words[0]]))*expected_identity*expected_length_shared
			else:
				j_exp=(float(ref_lengths[words[0]])/float(contigs[words[1]]["length"]))*expected_identity*expected_length_shared
			contigs[words[1]]["expected"]=(-1/float(kmer))*log((2*j_exp)/(1+j_exp))
			contigs[words[1]]["bounds"]=bounds(contigs[words[1]]["expected"], float(hashes.split("/")[1]), k=16, p=0.99)
			if contigs[words[1]]["expected"]<0:
				contigs[words[1]]["expected"]=0
			
			if contigs[words[1]]["expected"]>contigs[words[1]]["distance"]:
				contigs[words[1]]["match"]="yes"
			elif contigs[words[1]]["bounds"]==float("Inf"):
				contigs[words[1]]["match"]="uncertain"
			elif contigs[words[1]]["expected"]+contigs[words[1]]["bounds"]>contigs[words[1]]["distance"]:
				contigs[words[1]]["match"]="yes"
			else:
				contigs[words[1]]["match"]="no"
	
	

output = open(assembly+".chippy.tsv", "w")
print >> output, "\t".join(["contig", "contig length", "ref length", "top hit", "hash matches", "p-value", "distance", "expected", "cutoff", "chip?"])
	

for contig in contigs:
	print >> output, "\t".join(map(str,[contig, contigs[contig]["length"], contigs[contig]["ref_length"], contigs[contig]["hit"], contigs[contig]["hash_matches"], contigs[contig]["p"], contigs[contig]["distance"], contigs[contig]["expected"], contigs[contig]["distance"]-contigs[contig]["bounds"], contigs[contig]["match"]]))

output.close()



contig_sequences={"yes":[], "no":[], "uncertain":[]}
seq=[]
name=""
description=""
for line in open(assembly, "rU"):
        line=line.strip()
        if len(line)>0:
                if line[0]==">":
                        if name!="":
				contig_sequences[contigs[name]["match"]].append([name,description,''.join(seq)])
                        name=line.split()[0][1:]
                        seq=[]
                        if len(line.split())>1:
                                description=" ".join(line.split()[1:])
                else:
                        seq.append(line.replace(" ",""))
if name!="":
	contig_sequences[contigs[name]["match"]].append([name,description,''.join(seq)])


cout=open(assembly+".chippy.chips.fasta", "w")
fout=open(assembly+".chippy.fish.fasta", "w")
sout=open(assembly+".chippy.small_fries.fasta", "w")
if len(contig_sequences["yes"])>0:
	for seq in contig_sequences["yes"]:
		print >> cout, ">"+seq[0]+" "+seq[1]
		print >> cout, seq[2]
if len(contig_sequences["no"])>0:
	for seq in contig_sequences["no"]:
		print >> fout, ">"+seq[0]+" "+seq[1]
		print >> fout, seq[2]
	fout.close()
if len(contig_sequences["uncertain"])>0:
	for seq in contig_sequences["uncertain"]:
		print >> sout, ">"+seq[0]+" "+seq[1]
		print >> sout, seq[2]
	sout.close()
cout.close()
fout.close()
sout.close()
