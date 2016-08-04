#!/usr/bin/env python

from scipy import linalg
import numpy
import dendropy
import sys
import os
import time
from math import log
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_general import *
from Si_SeqIO import *

#Calculate Pij from Q matrix and branch length
def calculate_pij(branch_length,rate_matrix):
	return numpy.log(linalg.expm(numpy.multiply(branch_length,rate_matrix)))

#Read the tree file and root
def read_tree(treefile):
	t=dendropy.Tree.get(path=treefile, schema="newick", preserve_underscores=True, rooting="force-rooted")
	t.reroot_at_midpoint()
	t.update_bipartitions()
	return t


#Read the RAxML info file to get rates and frequencies
def read_info(infofile):
	f=[0,0,0,0]
	r=[0,0,0,0,0,0]
	for line in open(infofile, "rU"):
		line=line.strip()
		if "Base frequencies:" in line:
			words=line.split()
			f=map(float, words[2:6])
		elif "ac ag at cg ct gt:" in line:
			
			words=line.split()
			r=map(float, words[9:])
		elif "rate A <-> C:" in line:
			words=line.split()
			r[0]=float(words[4])
		elif "rate A <-> G:" in line:
			words=line.split()
			r[1]=float(words[4])
		elif "rate A <-> T:" in line:
			words=line.split()
			r[2]=float(words[4])
		elif "rate C <-> G:" in line:
			words=line.split()
			r[3]=float(words[4])
		elif "rate C <-> T:" in line:
			words=line.split()
			r[4]=float(words[4])
		elif "rate G <-> T:" in line:
			words=line.split()
			r[5]=float(words[4])
		elif "freq pi(A):" in line:
			words=line.split()
			f[0]=float(words[2])
		elif "freq pi(C):" in line:
			words=line.split()
			f[1]=float(words[2])
		elif "freq pi(G):" in line:
			words=line.split()
			f[2]=float(words[2])
		elif "freq pi(T):" in line:
			words=line.split()
			f[3]=float(words[2])
	return f, r


def create_rate_matrix(f, r):
	#convert f and r to Q matrix
	rm=numpy.array([[0, f[0]*r[1], f[0]*r[2], f[0]*r[3]],[f[1]*r[0], 0, f[1]*r[3],f[1]*r[4]],[f[2]*r[1], f[2]*r[3], 0, f[2]*r[5]],[f[3]*r[2], f[3]*r[4], f[3]*r[5], 0]])
	
	rm[0][0]=numpy.sum(rm[0])*-1
	rm[1][1]=numpy.sum(rm[1])*-1
	rm[2][2]=numpy.sum(rm[2])*-1
	rm[3][3]=numpy.sum(rm[3])*-1
	
	return rm

#Lookup for each base
mb={"A": 0, "C": 1, "G": 2, "T":3 }


# read the alignment
alignment=read_alignment(sys.argv[1])

#Create a new alignment for the output containing all taxa in the input alignment
new_alignment={}
for i, x in enumerate(alignment):
	new_alignment[x.id]=list(str(x.seq))

# read the tree
#tree=read_tree("RAxML_result.ml_sl.small.snps")
tree=read_tree(sys.argv[2])


#read the info file and get frequencids and rates
f, r=read_info(sys.argv[3])
if f=="" or r=="":
	print "Failed to read info file"
	sys.exit()

#create rate matrix from f and r
rm=create_rate_matrix(f,r)

#label internal nodes in tree and add these to the new alignment and calculate pij per non-root branch
mytree={}
nodecounter=0
myleaves=[]
myinternals=[]
preorder=[]
myroot=None
for node in tree.preorder_node_iter():

	if node.taxon==None:
		nodecounter+=1
		nodename='Node_'+str(nodecounter)
		tree.taxon_namespace.add_taxon(dendropy.Taxon(nodename))
		node.taxon=tree.taxon_namespace.get_taxon(nodename)
		if nodename in new_alignment:
			print nodename, "already in alignment. Quitting"
			sys.exit()
		new_alignment[nodename]=["?"]*len(alignment[0])
	if node.parent_node!=None:
		node.pij=calculate_pij(node.edge_length, rm)


#convert tree object to dictionary or nodes to make traversal much faster
		
for node in tree.preorder_node_iter():	
	taxon=str(node.taxon).strip("'")
	mytree[taxon]={'taxon': taxon, "children": [], 'L':{}, "C":{}, "parent": None, "is_leaf": False}
	node.L=nodecounter
	if node.parent_node!=None:
		mytree[taxon]["parent"]=str(node.parent_node.taxon).strip("'")
		mytree[taxon]["pij"]=node.pij
	else:
		myroot=taxon
	if node.is_leaf():
		myleaves.append(taxon)
		mytree[taxon]["is_leaf"]=True
	else:
		myinternals.append(taxon)
		for child in node.child_node_iter():
			mytree[taxon]["children"].append(str(child.taxon).strip("'"))
	preorder.append(taxon)

postorder=preorder[::-1]
postorderlen=len(postorder)
preorderlen=len(preorder)
leaflen=len(myleaves)
		

print "Alignment size:", len(alignment), "taxa and", len(alignment[0]), "sites"

#Find unique base patterns to speed up calculations
print "Finding unique base patterns"
base_patterns={}
t1=time.clock()
for x in xrange(len(alignment[0])):
	try:
		base_patterns[alignment[:,x]].append(x)
	except KeyError:
		base_patterns[alignment[:,x]]=[x]
t2=time.clock()
print "Time taken to find unique base patterns:", t2-t1, "seconds"
print "Unique base patterns:", len(base_patterns)



allbases=set(["A", "C", "G", "T"])

# For each base
t=0
timestart=time.clock()
onetime=0.0
twotime=0.0
threetime=0.0
for x, column in enumerate(base_patterns):
	t=t+1
	if t==1000:
		timeend=time.clock()
		print "Reconstructed", x+1, "of", len(base_patterns), "base patterns at", (timeend-timestart)/(x+1)*1000, "seconds per 1000 patterns"
		t=0
		#break
		
	
	columnbases=set([])
	base={}
	for i, y in enumerate(column):
		#print alignment[i].id, y
		base[str(alignment[i].id).strip("'")]=y
		if y in allbases:
			columnbases.add(y)
	
	#1 For each OTU y perform the following:

	#Visit a nonroot internal node, z, which has not been visited yet, but both of whose sons, nodes x and y, have already been visited, i.e., Lx(j), Cx(j), Ly(j), and Cy(j) have already been defined for each j. Let tz be the length of the branch connecting node z and its father. For each amino acid i, compute Lz(i) and Cz(i) according to the following formulae:
	
	#Denote the three sons of the root by x, y, and z. For each amino acid k, compute the expression Pk x Lx(k) x Ly(k) x Lz(k). Reconstruct r by choosing the amino acid k maximizing this expression. The maximum value found is the likelihood of the best reconstruction.
	
	
	onestart=time.clock()
	for y in xrange(postorderlen):
		nt=postorder[y]
		node=mytree[nt]
		if node["parent"]==None:
			continue
		#calculate the transistion matrix for the branch
		pij=node["pij"]
		
		if node["is_leaf"]:
			taxon=nt
			try:
				if base[taxon] in allbases:
					#1a. Let j be the amino acid at y. Set, for each amino acid i: Cy(i)= j. This implies that no matter w	hat is the amino acid in the father of y, j is assigned to node y.
					node["C"]={"A": base[taxon], "C": base[taxon], "G": base[taxon], "T": base[taxon]}
				
				#1b. Set for each amino acid i: Ly(i) = Pij(ty), where ty is the branch length between y and its father.
					node["L"]={"A": pij[mb["A"]][mb[base[taxon]]], "C": pij[mb["C"]][mb[base[taxon]]], "G": pij[mb["G"]][mb[base[taxon]]], "T": pij[mb["T"]][mb[base[taxon]]]}
				else:
					
					node["C"]={"A": "A", "C": "C", "G": "G", "T": "T"}
					node["L"]={"A": pij[mb["A"]][mb["A"]], "C": pij[mb["C"]][mb["C"]], "G": pij[mb["G"]][mb["G"]], "T": pij[mb["T"]][mb["T"]]}	
				
			except KeyError:
				print "Cannot find", taxon, "in base", x+1
				sys.exit()
		
		
		else:
			node["L"]={}
			node["C"]={}
			
			#2a. Lz(i) = maxj Pij(tz) x Lx(j) x Ly(j)
			#2b. Cz(i) = the value of j attaining the above maximum.
			
			for basenum in columnbases:
				node["L"][basenum]=float("-inf")
				node["C"][basenum]=None
			for end in columnbases:
				c=0.0
				for child in node["children"]:
					c+=mytree[child]["L"][end]
				for start in columnbases:
					j=pij[mb[start],mb[end]]+c
					
					
					if j>node["L"][start]:
						node["L"][start]=j
						node["C"][start]=end
					#print start, z, maxj
				
	node=mytree[myroot]
	node["L"]={}
	node["C"]={}
	for basenum in columnbases:
		node["L"][basenum]=float("-inf")
		node["C"][basenum]=None
	for end in columnbases:
		c=0
		for child in node["children"]:
			c+=mytree[child]["L"][end]
		for start in columnbases:
			j=log(f[mb[end]])+c

			if j>node["L"][start]:
				node["L"][start]=j
				node["C"][start]=end
		
	
	max_root_base=None
	max_root_base_likelihood=float("-inf")
	for root_base in columnbases:
		#print max_root_base, max_root_base_likelihood, root_base, node.L[root_base]
		if node["L"][root_base]>max_root_base_likelihood:
			max_root_base_likelihood=node["L"][root_base]
			max_root_base=node["C"][root_base]
	node["r"]=max_root_base
	#sys.exit()
	#print node.r, "root"
	#new_alignment[node.taxon.label].append(node.r)
	
	
	oneend=time.clock()
	onetime+=oneend-onestart
	
	
	
	#Traverse the tree from the root in the direction of the OTUs, assigning to each node its most likely ancestral character as follows:
	twostart=time.clock()
	for y in xrange(preorderlen):
		nt=preorder[y]
		node=mytree[nt]
		try:
			#5a. Visit an unreconstructed internal node x whose father y has already been reconstructed. Denote by i the reconstructed amino acid at node y.
			i=mytree[node["parent"]]["r"]
		except KeyError:
			continue
		#5b. Reconstruct node x by choosing Cx(i).
		node["r"]=node["C"][i]
		#new_alignment[node.taxon.label].append(node.r)

	twoend=time.clock()
	twotime+=twoend-twostart
	
	
	
	
	threestart=time.clock()
	#Put gaps back in and check that any ancestor with only gaps downstream is made a gap
	for y in xrange(postorderlen):
		nt=postorder[y]
		node=mytree[nt]
		if node["is_leaf"]:
			node["r"]=base[nt]
		else:
			has_child_base=False
			for child in node["children"]:
				if mytree[child]["r"] in allbases:
					has_child_base=True
					break
			if not has_child_base:
				node["r"]="-"
			for bp in base_patterns[column]:
				new_alignment[nt][bp]=node["r"]
	threeend=time.clock()
	threetime+=threeend-threestart
	

print onetime, twotime, threetime

output = open("jar.aln", "w")
for taxon in new_alignment:
	print >> output, ">"+str(taxon)
	print >> output, ''.join(new_alignment[taxon])
output.close()

output=open("jar.tre", "w")
print >> output, tree.as_string("newick")
output.close()


