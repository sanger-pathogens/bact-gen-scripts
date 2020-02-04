#!/usr/bin/env python

import os, sys
import cPickle
import numpy as np
import scipy.spatial.distance as ssd
import scipy.cluster.hierarchy as cluster
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist


def cluster_sequences(filename, cutoff=0.05, method='average', print_dendrogram=False):
	#cutoff=0.005
	#method='average' #clustering method. Can choose from 'single', 'complete', 'average', 'weighted'. Testing shows average performs best
	#print_dendrogram=True
	distances=[]
	matches=[]
	my_labels=[]
	subject=""
	for line in open(filename, "rU"):
		words=line.strip().split("\t")
		if words[1]!=subject:
			distances.append([])
			matches.append([])
			subject=words[1]
			my_labels.append(words[1])
		matches[-1].append(float(words[4].split("/")[0]))
		distances[-1].append(float(words[2]))
	Num_sequences=len(distances)
	print Num_sequences, "sequences found"
	X=np.array(distances)
	Y = ssd.squareform(X)
	Z=cluster.linkage(Y, method=method)
	c, coph_dists = cophenet(Z, pdist(X))
	num_clusters=Num_sequences
	clusters={}
	for x in xrange(Num_sequences):
		clusters[x]=[x]
	toremove=[]
	for z in Z:
		if z[2]>cutoff:
			break
		seq1=int(z[0])
		seq2=int(z[1])
		clusters[len(clusters)]=clusters[seq1]+clusters[seq2]
		toremove.append(seq1)
		toremove.append(seq2)
		#print z, seq1, seq2,len(clusters)
	for c in toremove:
		del clusters[c]
	cluster_reps={}
	for c in clusters:
		names=[]
		min_mean=float("Inf")
		if len(clusters[c])==1:
			cluster_reps[c]=my_labels[clusters[c][0]]
			names.append(my_labels[clusters[c][0]])
		else:
			for i in clusters[c]:
				names.append(my_labels[i])
				dists=[]
				for j in clusters[c]:
					if i!=j:
						dists.append(distances[i][j])
				mean_dist=np.average(dists)
				if mean_dist<min_mean:
					cluster_reps[c]=my_labels[i]
					min_mean=mean_dist
		clusters[c]=names
	print len(clusters), "clusters identified using", cutoff, "cutoff and", method, "method", cluster_reps
	
	return clusters, cluster_reps


def print_sequence_to_file(sequence, in_handle, out_handle):
	in_handle.seek(sequence["locations"][0])
	while in_handle.tell()<sequence["locations"][1]:
		line = in_handle.readline()
		print >> out_handle, line.strip()


try:
	f=open(sys.argv[1], "rU")
except:
	print "Could not open file "+sys.argv[1]
try:
	idx=open(sys.argv[1]+".chippy.index", "rU")
except:
	print "Could not open Chippy index file "+sys.argv[1]+".chippy.index"
	sys.exit()
index=cPickle.load(idx)
idx.close()
coutput=open(sys.argv[1]+".chromosomes1.chippy.cluster_reps.fasta", "w")
poutput=open(sys.argv[1]+".plasmids1.chippy.cluster_reps.fasta", "w")

if len(sys.argv)>2:
	species_list=[" ".join(sys.argv[2:])]
else:
	species_list=index.keys()
	
species_list.sort()
#print species_list
i=0
ptot=0
ctot=0
for species in species_list:
	if not species in index or not species[0].isupper():
		print species, "not in index. Skipping"
		continue
	print "found", len(index[species]), "sequences for", species
	i+=1
	coutname=species.replace(" ","_")+".chromosomes.chippy.temporary.fasta"
	poutname=species.replace(" ","_")+".plasmids.chippy.temporary.fasta"
	cout=open(coutname, "w")
	pout=open(poutname, "w")
	cs=[]
	ps=[]
	for accession in index[species]:
		if index[species][accession]["plasmid"]==True:
			print_sequence_to_file(index[species][accession], f, pout)
			ps.append(accession)
			if index[species][accession]["locations"][1]-index[species][accession]["locations"][0]>500000:
				print "Warning: Large plasmid", index[species][accession]["locations"][1]-index[species][accession]["locations"][0]
		else:
			print_sequence_to_file(index[species][accession], f, cout)
			cs.append(accession)
			if index[species][accession]["locations"][1]-index[species][accession]["locations"][0]<500000:
				print "Warning: small chromosome", index[species][accession]["locations"][1]-index[species][accession]["locations"][0]
			
		
	cout.close()
	pout.close()
	#os.system("~/scripts/chippy/chippy.py "+foutname+" "+foutname)
	print "found", len(ps), "plasmid sequences for", species
	print "found", len(cs), "chromosome sequences for", species
	files_to_cluster=[]
	if len(ps)>1:
		files_to_cluster.append(poutname)
	else:
		for accession in ps:
			print_sequence_to_file(index[species][accession], f, poutput)
			ptot+=1
	if len(cs)>1:
		files_to_cluster.append(coutname)
	else:
		for accession in cs:
			print_sequence_to_file(index[species][accession], f, coutput)
			ctot+=1
			
	for fname in files_to_cluster:
		os.system("mash sketch -i -k 16 -s 1000 "+fname)
		os.system("mash dist "+fname+".msh "+fname+".msh > "+fname+".dist.tab")
		clusters, representatives=cluster_sequences(fname+".dist.tab")
		for c in representatives:
			accession=representatives[c]
			if index[species][accession]["plasmid"]==True:
				print_sequence_to_file(index[species][accession], f, poutput)
				ptot+=1
			else:
				print_sequence_to_file(index[species][accession], f, coutput)
				ctot+=1
	noclean=False
	if noclean==True:
		os.system("rm "+coutname)
		os.system("rm "+poutname)
	
coutput.close()
putput.close()

	
f.close()
print "Found", ctot, "chromosomes and", ptot, "plasmids"
