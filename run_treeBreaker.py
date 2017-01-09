#!/usr/bin/env python3

#################################
# Import some necessary modules #
#################################

import os, sys
import dendropy
from optparse import OptionParser, OptionGroup
from random import randrange

##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
        
	group = OptionGroup(parser, "Required options")
	group.add_option("-t", "--tree", action="store", dest="tree", help="Input tree file in newick or nexus format", default="", metavar="FILE")
	group.add_option("-p", "--phenotype", action="store", dest="phenotype", help="Phenotype file (e.g. csv or tsv)", default="", metavar="FILE")
	group.add_option("-o", "--output_prefix", action="store", dest="prefix", help="Output prefix", default="treeBreaker")
	group.add_option("-c", "--column", action="store", dest="column", help="Number of the column in the phenotype file that contains the phenotype of interest [default= %default]", default=2, type="int")
	group.add_option("-s", "--separator", action="store", dest="separator", help="Separator for columns in input file [default= %default]", default=",")
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "treeBreaker options")
	group.add_option("-x", "--postburnin", action="store", dest="postburnin", help="Number of post-burnin iterations [default= %default]", default=500000, type="int")
	group.add_option("-y", "--burnin", action="store", dest="burnin", help="Number of burnin iterations [default= %default]", default=500000, type="int")
	group.add_option("-z", "--sample_iters", action="store", dest="sample_iters", help="Number of iterations between samples [default= %default]", default=1000, type="int")
	group.add_option("-S", "--seed", action="store", dest="seed", help="Starting seed for random number generator [default= %default]", default=-1, type="int")
	
	parser.add_option_group(group)
	group = OptionGroup(parser, "General options")
	group.add_option("-v", "--verbose", action="store_true", dest="verbose", help="Verbose [default= %default]", default=False)
	parser.add_option_group(group)
	
	return parser.parse_args()
	

def check_command_line_options(o):
	
	if  not os.path.isfile(o.tree):
		print("Cannot find file", o.tree)
		sys.exit()
		
	if  not os.path.isfile(o.phenotype):
		print("Cannot find file", o.phenotype)
		sys.exit()
	
	if o.column<2:
		print("Column number must be greater than 1")
		sys.exit()
	
	if o.postburnin<=0:
		print("treeBreaker post-brunin iterations should be greater than 0")
		sys.exit()
	if o.burnin<0:
		print("treeBreaker brunin iterations should be 0 or greater")
		sys.exit()
	if o.sample_iters<0:
		print("treeBreaker sample iterations should be 0 or greater")
		sys.exit()
	
	return


def read_tree(tfile):
	try:
		t=dendropy.Tree.get(path=tfile, schema="nexus", preserve_underscores=True)
	except dendropy.dataio.nexusreader.NexusReader.NotNexusFileError:
		t=dendropy.Tree.get(path=tfile, schema="newick", preserve_underscores=True)
	except:
		print("Cannot read tree file. It shoul dbe in newick or nexus format")
	return t

	
################
# Main program #
################                

if __name__ == "__main__":

	(options, args) = main()
	check_command_line_options(options)
	
	treefile=options.tree
		
	phenotypefile=options.phenotype
	
	output_prefix=options.prefix
	
	column=options.column-1
	
	separator=options.separator
	
	postburnin=options.postburnin
	
	burnin=options.burnin
	
	sample_iters=options.sample_iters
	
	if options.seed<0:
		seed=randrange(1,99999)
	else:
		seed=options.seed
	
	verbose=options.verbose
	
	treeBreakerOptions=["-x", postburnin, "-y", burnin, "-z", sample_iters, "-S", seed]
	if verbose:
		treeBreakerOptions.append("-v")
	treeBreakerOptionString=' '.join(map(str,treeBreakerOptions))
	
	tree=read_tree(treefile)
	
	tree_taxa=list(tree.taxon_namespace.labels())
	
	phenotypes=set([])
	taxon_phenotypes={}
	missing_phenotypes=[]
	not_in_tree=[]
	for line in open(phenotypefile, "rU"):
		words=line.strip().split(separator)
		if len(words)<column:
			if options.verbose:
				print("Invalid line in phenotype file (Too few columns):", line.strip())
		if words[column]=="":
			if options.verbose:
				print(words[0], "has missing phenotype")
			missing_phenotypes.append(words[0])
			continue
		if not words[0] in tree_taxa:
			if options.verbose:
				print(words[0], "not in tree. Skipping")
			not_in_tree.append(words[0])
			continue
		taxon_phenotypes[words[0]]=words[column]
		phenotypes.add(words[column])
	
	phenotype={}
	for x, p in enumerate(phenotypes):
		phenotype[p]=x
	
	print("Found", len(missing_phenotypes), "taxa without phenotypes")
	
	pruned_tree=tree.extract_tree_without_taxa_labels(labels=missing_phenotypes)
	
	pruned_tree.write(path=output_prefix+".nwk", schema="newick", unquoted_underscores=True)
	
	output=open(output_prefix+".tab", "w")
	for t in taxon_phenotypes:
		print(t+"\t"+str(phenotype[taxon_phenotypes[t]]), file=output)
	output.close()
	
	print("Phenotype key:")
	for p in phenotype:
		print(str(phenotype[p])+": "+p)
	
	print("Running treeBreaker")
	os.system("treeBreaker "+treeBreakerOptionString+" "+output_prefix+".nwk "+output_prefix+".tab "+output_prefix+".out")
	
	lines = open(output_prefix+".out","rU").readlines()
	treeBreaker_tree = lines[-1]
	
	treeBreaker=dendropy.Tree.get_from_string(treeBreaker_tree.replace("{", "[&").replace("}", "]").replace("|", ","), schema="newick")
	
	

	for node in treeBreaker.postorder_node_iter():
		if node==treeBreaker.seed_node:
			continue
		for a in node.annotations:
			if a.name=="posterior":
				leaves=[]
				for l in node.leaf_iter():
					leaves.append(l.taxon.label)
				if float(a.value)>0.5:
					print(a.value, len(leaves))

	treeBreaker.write(path=output_prefix+".nexus", schema="nexus", unquoted_underscores=True)
	
	
	os.system("~sh16/scripts/iCANDY.py -t "+output_prefix+".nexus -s posterior -J posterior -O portrait -p A0 -z circle -m "+phenotypefile+" -C "+str(column+1)+" -a 2 -o "+output_prefix+".pdf")
	
	
	
	sys.exit()
