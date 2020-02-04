#!/usr/bin/env python

import os, sys
import numpy as np
from sklearn.cluster import DBSCAN, AffinityPropagation
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from itertools import cycle
import scipy.spatial.distance as ssd
import scipy.cluster.hierarchy as cluster
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist

cutoff=0.005
method='average' #clustering method. Can choose from 'single', 'complete', 'average', 'weighted'. Testing shows average performs best
print_dendrogram=True

distances=[]
matches=[]
my_labels=[]
subject=""
for line in open(sys.argv[1], "rU"):
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
		print "cutoff"
		break
	seq1=int(z[0])
	seq2=int(z[1])
	clusters[len(clusters)]=clusters[seq1]+clusters[seq2]
	toremove.append(seq1)
	toremove.append(seq2)
	#print z, seq1, seq2,len(clusters)
for c in toremove:
	del clusters[c]
print len(clusters), "clusters identified using", cutoff, "cutoff and", method, "method"
print clusters

if print_dendrogram:
	# calculate full dendrogram
	plt.title('Hierarchical Clustering Dendrogram')
	plt.xlabel('sample index')
	plt.ylabel('distance')
	cluster.dendrogram(
	    Z,
	    leaf_rotation=90.,  # rotates the x axis labels
	    leaf_font_size=8.,  # font size for the x axis labels
	)
	plt.savefig(sys.argv[1]+'.dendrogram.png')
	
sys.exit()
	
X=np.array(matches)

af = AffinityPropagation().fit(X)
cluster_centers_indices = af.cluster_centers_indices_
labels = af.labels_
n_clusters_ = len(cluster_centers_indices)

print('Estimated number of clusters (affinity propagation): %d' % n_clusters_)

for x in xrange(len(labels)):
	print labels[x], my_labels[x]

sys.exit()

plt.close('all')
plt.figure(1)
plt.clf()

colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
for k, col in zip(range(n_clusters_), colors):
    class_members = labels == k
    cluster_center = X[cluster_centers_indices[k]]
    plt.plot(X[class_members, 0], X[class_members, 1], col + '.')
    plt.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=14)
    for x in X[class_members]:
        plt.plot([cluster_center[0], x[0]], [cluster_center[1], x[1]], col)

plt.title('Estimated number of clusters: %d' % n_clusters_)
plt.savefig("AffinityPropagation.png")
#print distances
X=np.array(distances)
db = DBSCAN(metric="precomputed").fit(X)
core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
core_samples_mask[db.core_sample_indices_] = True
labels = db.labels_
print labels
# Number of clusters in labels, ignoring noise if present.
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

print('Estimated number of clusters (DBSCAN): %d' % n_clusters_)
#print("Homogeneity: %0.3f" % metrics.homogeneity_score(labels_true, labels))
#print("Completeness: %0.3f" % metrics.completeness_score(labels_true, labels))
#print("V-measure: %0.3f" % metrics.v_measure_score(labels_true, labels))
#print("Adjusted Rand Index: %0.3f"
#      % metrics.adjusted_rand_score(labels_true, labels))
#print("Adjusted Mutual Information: %0.3f"
#      % metrics.adjusted_mutual_info_score(labels_true, labels))
#print("Silhouette Coefficient: %0.3f"
#      % metrics.silhouette_score(X, labels))



# Black removed and is used for noise instead.
unique_labels = set(labels)
colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
for k, col in zip(unique_labels, colors):
    if k == -1:
        # Black used for noise.
        col = 'k'

    class_member_mask = (labels == k)

    xy = X[class_member_mask & core_samples_mask]
    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=14)

    xy = X[class_member_mask & ~core_samples_mask]
    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=6)

plt.title('Estimated number of clusters: %d' % n_clusters_)
plt.savefig('clusters.png')
