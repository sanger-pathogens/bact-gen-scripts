#!/usr/bin/env python

import os, sys
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_general import *
import string
from random import *

alignment_file=sys.argv[1]

seqs={}
inread=False
clusterorder=[]
for line in open(alignment_file):
	words=line.strip().split()
	if len(words)==0:
		continue
	if words[0][0]==">":
		if inread:
			if not cluster in seqs:
				seqs[cluster]=[]
				clusterorder.append(cluster)
				
			seqs[cluster].append([len(sequence), name, ''.join(sequence)])
		name=words[0][1:]
		cluster=words[0][1:].split("_")[0]
		sequence=[]
		inread=True
	else:
		sequence.append(words[0].upper())
if inread:
	if not cluster in seqs:
		seqs[cluster]=[]
		clusterorder.append(cluster)
	seqs[cluster].append([len(sequence), name, ''.join(sequence)])


output=open(sys.argv[1]+"_clusters", "w")
outputb=open(sys.argv[1]+"_clusters.clstr", "w")
for x in clusterorder:
	cluster=seqs[x]
	cluster.sort()
	cluster.reverse()
	print >> output, ">"+cluster[0][1]
	print >> output, cluster[0][2]
	print >> outputb, ">Cluster", x
	for i in xrange(0,len(cluster)):
		if i == 0:
			print >> outputb, '\t'.join(map(str, [i, str(cluster[i][0])+"nt", ">"+cluster[i][1]+"...", "*"]))
		else:
			print >> outputb, '\t'.join(map(str, [i, str(cluster[i][0])+"nt", ">"+cluster[i][1]+"...", "at +/100%"]))
output.close()
outputb.close()