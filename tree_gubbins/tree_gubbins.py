#!/software/python-2.7.6/bin/python
##!/usr/bin/env python

#################################
# Import some necessary modules #
#################################

import os, sys
import dendropy
import math
import copy
import random
from optparse import OptionParser

##########################
# Error message function #
##########################

def DoError(errorstring):
	print "\nError:", errorstring
	print "\nFor help use -h or --help\n"
	sys.exit()

##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	
	parser.add_option("-t", "--tree", action="store", dest="tree", help="Input tree file", default="", metavar="FILE")
	parser.add_option("-o", "--output_prefix", action="store", dest="output", help="Output prefix", default="")
	parser.add_option("-s", "--significance", action="store", dest="significance", help="Significance butoff level [default= %default]", default=0.05, type="float", metavar="FLOAT")
	parser.add_option("-p", "--permutations", action="store", dest="permutations", help="Number of permutations to run to test significance [default= %default]", default=100, type="int", metavar="INT")
	


	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


	if options.tree=='':
		DoError('No tree file selected')
	elif not os.path.isfile(options.tree):
		DoError('Cannot find file '+options.tree)
		
	if options.output=='':
		DoError('No output prefix selected')
	
	if options.significance<0 or options.significance>1:
		DoError('Significance must be between 0 and 1')
	
	
	if options.permutations<(1.0/options.significance):
		DoError('Number of permutations must allow significance to be reached')
	
	return




def read_tree(treefile):
	#Try opening the tree using various schemas
	opened=False
	for treeschema in ["nexus", "newick"]:#["beast-summary-tree",  "nexus", "newick"]:
		try: 
			t = dendropy.Tree.get_from_path(treefile, schema=treeschema, as_rooted=True, preserve_underscores=True, case_insensitive_taxon_labels=False, set_node_attributes=True, extract_comment_metadata=True)
			opened=True
			t.schema=treeschema
			break
		except dendropy.utility.error.DataParseError:
			continue
		except ValueError as e:
			print "Encountered ValueError while trying to read tree file as", treeschema+":", e
			continue
			
	if not opened:
		print "Failed to open tree file"
		sys.exit()
	
	t.deroot()
	
	return t


def get_clade_lengths(t, treelength, node_count):
	for node in t.postorder_node_iter():
		
		node.downstream_length=0
		node.downstream_count=0
		node.downstream_set=set([])
		
		if node.is_internal():
			for child in node.child_nodes():
				node.downstream_length+=child.downstream_length
				node.downstream_length+=child.edge_length
				node.downstream_count+=child.downstream_count
				node.downstream_count+=1
				node.downstream_set.update(child.downstream_set)
		else:
			node.downstream_set.add(node.taxon.label)
			
	return t



def getlikelihood(N, C, n, c):
	#Where N = treelength, C = node_count, n=clade_length, c= clade size

#	print c, n, C, N
#	
#	print C-c, N-n,
	
	part1=math.log((c/n),10)*c
	if n-c==0:
		part2=0
	else:
		part2=math.log((((n-c)/n)),10)*(n-c)
	if C-c==0:
		part3=0
	else:
		part3=math.log((((C-c)/(N-n))),10)*(C-c)
	if ((N-n)-(C-c))==0:
		part4=0
	else:
		part4=math.log(((((N-n)-(C-c))/(N-n))),10)*((N-n)-(C-c))
	
	likelihood=(part1+part2+part3+part4)*-1
	
#	print c, n, C, N, likelihood
	
	return likelihood


def get_tree_node_likelihoods(t, treelength, node_count, min_length):
	likelihoods=[]
	for node in t.postorder_internal_node_iter():
		if node.edge_length!=None:
			node.upstream_length=treelength-(node.downstream_length+node.edge_length)
		else:
			continue
		node.upstream_count=(node_count-node.downstream_count)-1
#		print node_count, node.downstream_count, node.upstream_count, min_length
		if node.downstream_count>=(node_count-2) or node.upstream_count>=(node_count-2) or node.downstream_count<2 or node.upstream_count<2:
			continue
		node.downstream_mean=node.downstream_length/node.downstream_count
		node.upstream_mean=node.upstream_length/node.upstream_count
		
		if((1/min_length)*treelength/node_count)>((1/min_length)*node.downstream_length/node.downstream_count):
			node.downstream_likelihood=getlikelihood((1/min_length)*treelength, node_count, (1/min_length)*node.downstream_length, node.downstream_count)
			likelihoods.append([node.downstream_likelihood, node, "d"])
		if ((1/min_length)*treelength/node_count)>((1/min_length)*node.upstream_length/node.upstream_count):
			node.upstream_likelihood=getlikelihood((1/min_length)*treelength, node_count, (1/min_length)*node.upstream_length, node.upstream_count)
			likelihoods.append([node.upstream_likelihood, node, "u"])
	
	return likelihoods


################
# Main program #
################		

if __name__ == "__main__":

	#Get command line arguments

	(options, args) = main()
	
	#Do some checking of the input files
	
	check_input_validity(options, args)

	
	try:
		original_tree=read_tree(options.tree)
	except StandardError:
		DoError("Failed to read tree file")
	
	clusters=[]
	
	tree=copy.deepcopy(original_tree)
	
	
	significant=True
	i=1
	while significant:
		
		print "Iteration", i
		
		node_count=0
		
		taxa=set([])
		
		for node in tree.leaf_iter():
			taxa.add(node.taxon.label)
		
		if len(taxa)<4:
			break
		
		blengths=[]
		min_length=float("Inf")
		for node in tree.postorder_node_iter():
			node_count+=1
			if node.edge_length!=None:
				blengths.append(node.edge_length)
				if node.edge_length==0:
					node.edge_length=0.00000000001
				if node.edge_length<min_length:
					min_length=node.edge_length
		min_length=min_length/10
		treelength=tree.length()
		
		get_clade_lengths(tree, treelength, node_count)
		
		
		likelihoods=get_tree_node_likelihoods(tree, treelength, node_count, min_length)
		
		likelihoods.sort()
		
		
		test_values=[]
		for x in xrange(0,options.permutations-1):
			new_tree=copy.deepcopy(tree)
			rand_blengths=blengths[:]
			random.shuffle(rand_blengths)
			x=xrange(int((1/min_length)*treelength))
			bits=random.sample(x, node_count-1)
			
			
			bits.sort()
			
			lengths=[]
			lastbit=0
			for bit in bits:
				bitlen=float(bit)/(1/min_length)
				lengths.append(bitlen-lastbit)
				lastbit=bitlen
			lengths.append(treelength-lastbit)
			
			
			y=0
			for node in new_tree.postorder_node_iter():
				if node.edge_length!=None:
					node.edge_length=lengths[y]
					y+=1
				
			get_clade_lengths(new_tree, treelength, node_count)
			tlikelihoods=get_tree_node_likelihoods(new_tree, treelength, node_count, min_length)
			tlikelihoods.sort()
			try:
				test_values.append(tlikelihoods[0][0])
			except StandardError:
				print tlikelihoods
		
		
		node=likelihoods[0][1]
		
		bettercount=0.0
		for t in test_values:
			if likelihoods[0][0]<t:
				bettercount+=1
		
		pvalue=(options.permutations-bettercount)/options.permutations
		
		if likelihoods[0][2]=="d":	
			cluster=likelihoods[0][1].downstream_set
		else:
			cluster=taxa.difference(likelihoods[0][1].downstream_set)
		
		
		if pvalue<=options.significance:
			clusters.append([i, cluster, pvalue])
			significant=True
			i+=1
			taxa_to_prune=[]
			tree.prune_taxa_with_labels(list(cluster))
			print "Found significant cluster at", pvalue, "level"
			print "Cluster contains", len(list(cluster)), "taxa:"
			print ", ".join(list(cluster))
		else:
			significant=False
			print "No significant clusters found"
		
		
		
	print "found a total of", len(clusters), "significant clusters"
	if len(clusters)>0:
		print "Printing output csv"
		output=open(options.output+".csv", "w")
		print >> output, ",".join(["Taxon", "Cluster", "p-value:c:0:"+str(options.significance)])
		for cluster in clusters:
			for taxon in list(cluster[1]):
				print >> output, ",".join(map(str,[taxon, cluster[0], cluster[2]]))
		output.close()
			
		print "Drawing tree"
		os.system("~/scripts/iCANDY.py -t "+options.tree+" -M -L right -m "+options.output+".csv -a 2 -C 2,2,3 -r deltran -O portrait -o "+options.output+".pdf")
	else:
		print "No clusters found, so no output to print"
	
	sys.exit()