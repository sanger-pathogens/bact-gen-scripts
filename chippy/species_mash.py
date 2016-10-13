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
	print len(clusters), "clusters identified using", cutoff, "cutoff and", method, "method"
	
	return clusters, cluster_reps


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
output=open(sys.argv[1]+".chippy.cluster_reps.fasta", "w")
i=0
for species in index:
	print "found", len(index[species]), "sequences for", species
	if len(index[species])>1:
		i+=1
		foutname=species.replace(" ","_")+".chippy.temporary.fasta"
		fout=open(foutname, "w")
		for accession in index[species]:
			f.seek(index[species][accession][0])
			while f.tell()<index[species][accession][1]:
				line = f.readline()
				print >> fout, line.strip()
		fout.close()
		#os.system("~/scripts/chippy/chippy.py "+foutname+" "+foutname)
		os.system("mash sketch -i -k 16 -s 1000 "+foutname)
		os.system("mash dist "+foutname+".msh "+foutname+".msh > "+foutname+".dist.tab")
		clusters, representatives=cluster_sequences(foutname+".dist.tab")
		for c in representatives:
			accession=representatives[c]
			f.seek(index[species][accession][0])
			while f.tell()<index[species][accession][1]:
				line = f.readline()
				print >> output, line.strip()
	else:
		for accession in index[species]:
			f.seek(index[species][accession][0])
			while f.tell()<index[species][accession][1]:
				line = f.readline()
				print >> output, line.strip()
	
output.close()

	
f.close()
