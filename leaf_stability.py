#!/usr/bin/env python

#################################
# Import some necessary modules #
#################################

import dendropy
import os, sys

def print_help():
	print "leaf_stability.py"
	print
	print "Calculates stability of leaves in bootstrap trees"
	print
	print"Usage:"
	print "\tleaf_stability.py <bootstrap tree file>"
	print
	sys.exit()

def read_tree_file(treefile):

	opened=False
	for treeschema in ["beast-summary-tree", "nexus", "newick"]:
		try:
			t = dendropy.TreeList.get_from_path(treefile, schema=treeschema, as_rooted=True)
			opened=True
			break
		except dendropy.utility.error.DataParseError:
			continue
	if not opened:
		print "Failed to open tree file"
		sys.exit()
	
	return t


def get_split_frequencies(trees, taxa, split_type="quartets"):

	splits={}

	for tree in trees:
		root=tree.seed_node
		root_children = root.child_nodes()
		if len(root_children) != 2:
			tree.reroot_at_edge(root_children[0].edge, update_splits=True)

		leaves=[]
		for leaf in tree.leaf_iter():
			#leaves.append(leaf.taxon.get_bitmask())
			leaves.append(taxa.taxon_bitmask(leaf.taxon))

		leaves.reverse()
		print leaves
		for node in tree.level_order_node_iter():
			
			if not node.is_internal():
				print "l", taxa.taxon_bitmask(node.taxon)
				leaves.pop()
				continue
			print leaves
			for child in node.child_nodes():
				if not child.is_internal():
					print "c", taxa.taxon_bitmask(child.taxon)
			continue
			children=node.child_nodes()
			
			if len(children)>2:
				print "Trees must be bifurcating"
				sys.exit()
			
			for leaf1 in children[0].leaf_iter():
				for leaf2 in children[1].leaf_iter():

					
					
					print taxa.taxon_bitmask(leaf1.taxon), taxa.taxon_bitmask(leaf2.taxon), taxa.get_taxa_bitmask(taxa=dendropy.TaxonSet([leaf1.taxon, leaf2.taxon])), taxa.taxon_bitmask(leaf1.taxon)+taxa.taxon_bitmask(leaf2.taxon)
		

def get_quartet_frequencies(trees, taxa):
	
	quartets={}
	
	for treenum, tree in enumerate(trees):
		print treenum
		
		root = tree.seed_node
		root_children = root.child_nodes()
		if len(root_children) != 2:
			tree.reroot_at_edge(root_children[0].edge, update_splits=True)

		donequartets={}
		
		for node in tree.postorder_internal_node_iter():
		
			children=node.child_nodes()
			
			if len(children)>2:
				print "Trees must be bifurcating"
				sys.exit()
				
			node_taxa=dendropy.TaxonSet()
			for leaf in node.leaf_iter():
				node_taxa.append(leaf.taxon)
				
			
			
			node_taxa_bitmask=taxa.get_taxa_bitmask(taxa=node_taxa)
			complement_taxa_bitmask=taxa.complement_split_bitmask(node_taxa_bitmask)
			complement_taxa=taxa.split_taxa_list(complement_taxa_bitmask)
			
			for leaf1 in children[0].leaf_iter():
				for leaf2 in children[1].leaf_iter():
					inset=dendropy.TaxonSet([leaf1.taxon, leaf2.taxon])
					inmask=taxa.get_taxa_bitmask(taxa=inset)
					
					for x, leaf3 in enumerate(complement_taxa):
						for leaf4 in complement_taxa[x+1:]:
							outset=dendropy.TaxonSet([leaf3, leaf4])
							
							outmask=taxa.get_taxa_bitmask(taxa=outset)
							
							
							if tree.weight is None:
								weight_to_use = 1.0
							else:
								weight_to_use = float(tree.weight)
							
							if outmask>inmask:
								if not outmask in quartets:
									quartets[outmask]={}
								if not inmask in quartets[outmask]:
									quartets[outmask][inmask]=[0.0,0.0]
								if not outmask in donequartets:
									donequartets[outmask]=set([inmask])
									quartets[outmask][inmask][0]+=1.0
									quartets[outmask][inmask][1]+=weight_to_use
								elif not inmask in donequartets[outmask]:
									donequartets[outmask].add(inmask)
									quartets[outmask][inmask][0]+=1.0
									quartets[outmask][inmask][1]+=weight_to_use
									
							else:
								if not inmask in quartets:
									quartets[inmask]={}
								if not outmask in quartets[inmask]:
									quartets[inmask][outmask]=[0.0,0.0]
								if not inmask in donequartets:
									donequartets[inmask]=set([outmask])
									quartets[inmask][outmask][0]+=1.0
									quartets[inmask][outmask][1]+=weight_to_use
								elif not outmask in donequartets[inmask]:
									donequartets[inmask].add(outmask)
									quartets[inmask][outmask][0]+=1.0
									quartets[inmask][outmask][1]+=weight_to_use
									
		
	return quartets


def comb(N,k):
	if (k > N) or (N < 0) or (k < 0):
		return 0L
	N,k = map(long,(N,k))
	top = N
	val = 1L
	while (top > (N-k)):
		val *= top
		top -= 1
	n = 1L
	while (n < k+1L):
		val /= n
		n += 1
	return val
    

def quartets_to_positional_congruence(quartets, taxa, numtrees):
	
	positional_congruence=0.0
	nquartets=0.0
	leaf_positional_congruence={}
	for taxon in taxa:
		leaf_positional_congruence[taxon.label]=[0.0,0.0]
	
	for outside in quartets:
		for inside in quartets[outside]:
			nquartets+=1
			for taxon in map(lambda x: x.label,taxa.split_taxa_list(outside)+taxa.split_taxa_list(inside)):
				leaf_positional_congruence[taxon][0]+=1
			if quartets[outside][inside][0]==numtrees:
				positional_congruence+=1
				for taxon in map(lambda x: x.label,taxa.split_taxa_list(outside)+taxa.split_taxa_list(inside)):
					leaf_positional_congruence[taxon][1]+=1
	
	positional_congruence=positional_congruence/(comb(len(taxa),4))
	for taxon in taxa:
		print taxon.label, leaf_positional_congruence[taxon.label][1]/(comb(len(taxa)-1,3))
	print "Total", positional_congruence
	
			
	print positional_congruence, leaf_positional_congruence


def quartets_to_leaf_stability(quartets, taxa, numtrees):
	
	stability=0.0
	nquartets=0.0
	leaf_stability={}
	for taxon in taxa:
		leaf_stability[taxon.label]=[0.0,0.0]
	
	for outside in quartets:
		for inside in quartets[outside]:
			nquartets+=1
			for taxon in map(lambda x: x.label,taxa.split_taxa_list(outside)+taxa.split_taxa_list(inside)):
				leaf_stability[taxon][0]+=1
			if quartets[outside][inside][0]==numtrees:
				stability+=1
				for taxon in map(lambda x: x.label,taxa.split_taxa_list(outside)+taxa.split_taxa_list(inside)):
					leaf_stability[taxon][1]+=quartets[outside][inside][1]/numtrees
	
	stability=stability/(comb(len(taxa),4))
	for taxon in taxa:
		print taxon.label, (leaf_stability[taxon.label][1]/leaf_stability[taxon.label][0])#/(comb(len(taxa)-1,3))
	print "Total", stability
	
			
	print leaf_stability, stability
	

if len(sys.argv)!=2:
	print "Expecting one argument (bootstrap tree file)"
	print_help()
if sys.argv[1]=="-h":
	print_help()
if not os.path.isfile(sys.argv[1]):
	print "Cannot find file", sys.argv[1]
	print_help()
	
	
trees=read_tree_file(sys.argv[1])
numtrees=len(trees)
taxa=trees[0].taxon_set

print comb(len(taxa),4)
print comb(8,4)*3
#sys.exit()

get_split_frequencies(trees,taxa,"quartets")
sys.exit()
quartets=get_quartet_frequencies(trees, taxa)

quartets_to_positional_congruence(quartets, taxa, numtrees)


leaf_stability=quartets_to_leaf_stability(quartets, taxa, numtrees)
sys.exit()
	
#cons=trees.consensus()
#print(cons.as_ascii_plot(show_internal_node_labels=True))


