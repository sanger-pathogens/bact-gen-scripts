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
	
	parser.add_option("-t", "--tree", action="store", dest="tree", help="Input tree file (NOTE: tree must be fully bifurcating)", default="", metavar="FILE")
	parser.add_option("-o", "--output_prefix", action="store", dest="output", help="Output prefix", default="")
	parser.add_option("-s", "--significance", action="store", dest="significance", help="Significance cutoff level [default= %default]", default=0.05, type="float", metavar="FLOAT")
	parser.add_option("-p", "--permutations", action="store", dest="permutations", help="Number of permutations to run to test significance [default= %default]", default=100, type="int", metavar="INT")
	parser.add_option("-m", "--midpoint", action="store_true", dest="midpoint", help="Midpoint root output tree pdf [default= %default]", default=False)
	parser.add_option("-l", "--ladderise", action="store", choices=['right', 'left'], dest="ladderise", help="ladderise tree (choose from right or left) [default= %default]", type="choice", default=None)
	parser.add_option("-i", "--iterative", action="store_true", dest="iterative", help="Run in iterative mode, which allows nested clusters. Default = do not run iteratively (much faster)", default=False)
	

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
			t = dendropy.Tree.get_from_path(treefile, schema=treeschema, rooting='force-rooted', preserve_underscores=True, extract_comment_metadata=True)
			opened=True
			t.schema=treeschema
			break
		except dendropy.utility.error.DataParseError:
			print "DataParseError"
			continue
		except StandardError as e:
			print "Encountered ValueError while trying to read tree file as", treeschema+":", e
			continue
			
	if not opened:
		print "Failed to open tree file"
		sys.exit()
	
	#t.deroot()
	
	
	return t


def get_clade_lengths(t, treelength, node_count):
	for node in t.postorder_node_iter():
		
		node.downstream_length=0.0
		node.downstream_count=0.0
		node.downstream_set=set([])
		
		if node.is_internal():
			
			if len(node.child_nodes())>2:
				max_child_node_value=0
				max_child=0
				for x, child in enumerate(node.child_nodes()):
					if child.downstream_length+child.edge_length>max_child_node_value:
						max_child_node_value=child.downstream_length+child.edge_length
						max_child=x
				for x, child in enumerate(node.child_nodes()):
					if x==max_child:
						continue
					node.downstream_length+=child.downstream_length
					node.downstream_length+=child.edge_length
					node.downstream_count+=child.downstream_count
					node.downstream_count+=1
					node.downstream_set.update(child.downstream_set)
			else:
				for child in node.child_nodes():
					node.downstream_length+=child.downstream_length
					node.downstream_length+=child.edge_length
					node.downstream_count+=child.downstream_count
					node.downstream_count+=1
					node.downstream_set.update(child.downstream_set)
		else:
			node.downstream_set.add(node.taxon.label)
		#print float(node.downstream_count)/2, len(node.downstream_set)
			
	return t



def getlikelihood(N, C, n, c):
	#Where N = treelength, C = node_count, n=clade_length, c= clade size
	
	
#	print c, n, c/n, C, N, C/N, C-c, N-n
	

	try:
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
	except ValueError:
		print "Failed to calculate likelihoods"
		print c, n, C, N
		sys.exit()
	
#	print c, n, C, N, likelihood
	
	return likelihood


def get_tree_node_likelihoods(t, treelength, node_count, multiplier, verbose=False):
	likelihoods=[]
	for node in t.postorder_internal_node_iter():
		if node.edge_length!=None:
			node.upstream_length=treelength-(node.downstream_length+node.edge_length)
		else:
			continue
		node.upstream_count=float(node_count-node.downstream_count)-1
		
		
		if (node.downstream_count)>1 and (node.upstream_count+1)>0 and (multiplier*treelength)/node_count>((multiplier)*node.downstream_length)/node.downstream_count:
			node.downstream_likelihood=getlikelihood(multiplier*treelength, node_count, (multiplier)*node.downstream_length, node.downstream_count)
			likelihoods.append([node.downstream_likelihood, node, "d", multiplier, node.downstream_length, node.downstream_count])
			if verbose:
				print node, len(node.downstream_set), node.downstream_length, node.downstream_count, node.downstream_likelihood
		if (node.upstream_count)>1 and (node.downstream_count+1)>0 and (multiplier*treelength)/node_count>((multiplier)*node.upstream_length)/node.upstream_count:
			node.upstream_likelihood=getlikelihood(multiplier*treelength, node_count, (multiplier)*node.upstream_length, node.upstream_count)
			likelihoods.append([node.upstream_likelihood, node, "u", multiplier, node.upstream_length, node.upstream_count])
			if verbose:
				print node, (float(node_count)/2)-len(node.downstream_set), node.upstream_length, node.upstream_count, node.upstream_likelihood
	
	return likelihoods


################
# Main program #
################		

if __name__ == "__main__":

	#Get command line arguments

	(options, args) = main()
	
	#Do some checking of the input files
	
	check_input_validity(options, args)
	
	if options.iterative:
		print "Running in iterative mode"

	
	try:
		original_tree=read_tree(options.tree)
	except StandardError:
		DoError("Failed to read tree file")
	
	clusters=[]
	singletons=[]
	
	tree=copy.deepcopy(original_tree)
	
	significant=True
	i=1
	while significant:
		significant=False
		if options.iterative:
			print "Iteration", i
		
		node_count=0.0
		
		taxa=set([])
		
		tree.deroot()
		
		for node in tree.leaf_iter():
			taxa.add(node.taxon.label)
		
		if len(taxa)<3:
			print "Too few taxa to cluster"
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
		multiplier=node_count/min_length
		treelength=float(tree.length())
		
		tree=get_clade_lengths(tree, treelength, node_count)
		
		
		likelihoods=get_tree_node_likelihoods(tree, treelength, node_count, multiplier, verbose=False)
		
		likelihoods.sort()
		
		
		if len(likelihoods)==0:
			print "\tNo significant clusters found"
			break
				
		
#		for likelihood in likelihoods:
#			if likelihood[2]=="d":	
#				cluster=list(likelihood[1].downstream_set)
#			else:
#				cluster=list(taxa.difference(likelihood[1].downstream_set))
#			print likelihood[0], cluster
#		sys.exit()
		test_values=[]
		for x in xrange(0,options.permutations-1):
			new_tree=copy.deepcopy(tree)
			rand_blengths=blengths[:]
			random.shuffle(rand_blengths)
			x=xrange(int((1/min_length)*treelength))
			bits=random.sample(x, int(node_count)-1)
			
			
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
				
			new_tree=get_clade_lengths(new_tree, treelength, node_count)
			tlikelihoods=get_tree_node_likelihoods(new_tree, treelength, node_count, multiplier)
			tlikelihoods.sort()
			try:
				test_values.append(tlikelihoods[0][0])
			except StandardError:
				test_values.append(float("Inf"))
		
		test_values.sort()
		test_values.reverse()
		
		removedset=set([])
#		for likelihood in likelihoods[:1]:
		for x, likelihood in enumerate(likelihoods):
			if options.iterative and x>0:
				break
			node=likelihood[1]
			
			bettercount=0.0
			t=0
			while t<len(test_values) and likelihood[0]<test_values[t]:
				bettercount+=1
				t+=1
			
			pvalue=(float(options.permutations)-bettercount)/options.permutations
			
			#correction for multiple testing
			
			if likelihood[2]=="d":	
				cluster=likelihood[1].downstream_set
			else:
				cluster=taxa.difference(likelihood[1].downstream_set)
			
			
			if len(cluster.intersection(removedset))>0:
				continue
			
			if pvalue<=options.significance:
				clusters.append([i+x, cluster, pvalue, likelihood[0], likelihood[3], likelihood[4], likelihood[5]])
				significant=True
				
				taxa_to_prune=[]
				try:
					tree.prune_taxa_with_labels(list(cluster))
				except StandardError:
					print "Could not prune taxa"
				print "\tFound significant cluster at", pvalue, "level"
				print "\t\tCluster contains", len(list(cluster)), "taxa:"
				print "\t\t"+", ".join(list(cluster))
				removedset.update(cluster)
				
				
				
		if not significant:
			print "\tNo significant clusters found"
			if len(likelihoods)>0:
				print "\tTop cluster has p-value of", pvalue
			for j, taxon in enumerate(taxa):
				singletons.append([taxon, i+j])
		i+=1
		if not options.iterative:
			break
		
		
		
	print "Found a total of", len(clusters), "significant clusters"
	if len(clusters)>0:
		print "Printing output csv"
		output=open(options.output+".csv", "w")
		if options.iterative:
			iterative="Iteration"
		else:
			iterative="Cluster"
		print >> output, ",".join(["Taxon", iterative, "p-value:c:"+str(1.0/options.permutations)+":"+str(options.significance), "Likelihood:c", "Multiplier:c", "Length:c", "Count:c"])
		for x, cluster in enumerate(clusters):
			for taxon in list(cluster[1]):
				print >> output, ",".join(map(str,[taxon, cluster[0], cluster[2], cluster[3], cluster[4], cluster[5], cluster[6],  cluster[5]/cluster[6]]))
		output.close()
			
		print "Drawing tree"
		
		if options.midpoint:
			mid="-M"
		else:
			mid=""
		
		if options.ladderise==None:
			lad=""
		else:
			lad="-L "+options.ladderise
		if options.iterative:
			os.system("~sh16/scripts/iCANDY.py -t "+options.tree+" "+mid+" "+lad+" -m "+options.output+".csv -a 2 -C 2,2,3,6,7 -r deltran -O portrait -o "+options.output+".pdf")
		else:
			os.system("~sh16/scripts/iCANDY.py -t "+options.tree+" "+mid+" "+lad+" -m "+options.output+".csv -a 2 -C 2,2,4,3,6,7 -r deltran -O portrait -o "+options.output+".pdf")
		print "Printing PLINK output file"
		output=open(options.output+"_plink.txt", "w")
		
		for x, cluster in enumerate(clusters):
			for taxon in list(cluster[1]):
				print >> output, "\t".join(map(str,[taxon, taxon, "cluster_"+str(cluster[0])]))
		for x in singletons:
			print >> output, "\t".join(map(str,[x[0], x[0], "cluster_"+str(x[1])]))
		output.close()
	else:
		print "No clusters found, so no output to print"
	
	sys.exit()
