#!/usr/bin/env python3

#################################
# Import some necessary modules #
#################################

import os, sys
import dendropy
from optparse import OptionParser

##########################################
# Function to Get command line arguments #
##########################################


def main():

        usage = "usage: %prog [options]"
        parser = OptionParser(usage=usage)
        
        parser.add_option("-t", "--tree", action="store", dest="tree", help="Input tree file (NOTE: tree must be fully bifurcating)", default="", metavar="FILE")
        parser.add_option("-o", "--output_prefix", action="store", dest="prefix", help="Output prefix", default="treeBreaker")
        parser.add_option("-p", "--phenotype", action="store", dest="phenotype", help="Phenotype file cvs or tsv", default="", metavar="FILE")
        parser.add_option("-c", "--column", action="store", dest="column", help="Numer of the column in the phenotype file that contains the phenotype of interest [default= %default]", default=2, type="int")
        parser.add_option("-T", "--tab", action="store_true", dest="tab", help="phenotype file is tab separated (rather than comma separated [default= %default]", default=False)
        
        return parser.parse_args()
	
	
################
# Main program #
################                

if __name__ == "__main__":

	(options, args) = main()
	
	treefile=options.tree
	if  not os.path.isfile(treefile):
		print("Cannot find file", treefile)
		sys.exit()
	phenotypefile=options.phenotype
	if  not os.path.isfile(phenotypefile):
		print("Cannot find file", phenotypefile)
		sys.exit()
	output_prefix=options.prefix
	column=options.column-1
	if column<1:
		print("Column number must be greater than 1")
		sys.exit()
	separator=","
	if options.tab:
		separator="\t"
	
	tree=dendropy.Tree.get(path=treefile, schema="newick", preserve_underscores=True)
	
	tree_taxa=list(tree.taxon_namespace.labels())
	
	phenotypes=set([])
	taxon_phenotypes={}
	missing_phenotypes=[]
	not_in_tree=[]
	for line in open(phenotypefile, "rU"):
		words=line.strip().split(separator)
		if len(words)<column:
			print("Invalid line in phenotype file (Too few columns):", line.strip())
		if words[column]=="":
			print(words[0], "has missing phenotype")
			missing_phenotypes.append(words[0])
			continue
		if not words[0] in tree_taxa:
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
		
	os.system("treeBreaker "+output_prefix+".nwk "+output_prefix+".tab "+output_prefix+".out")
	
	lines = open(output_prefix+".out","rU").readlines()
	treeBreaker_tree = lines[-1]
	
	treeBreaker=dendropy.Tree.get_from_string(treeBreaker_tree.replace("{", "[&").replace("}", "]").replace("|", ","), schema="newick")
	
	treeBreaker.write(path=output_prefix+".tre", schema="nexus", unquoted_underscores=True)
	
	os.system("/nfs/pathogen/sh16_scripts/iCANDY.py -t "+output_prefix+".tre -s posterior -J posterior -O portrait -p A0 -z circle -m "+phenotypefile+" -C "+str(column+1)+" -a 2 -o "+output_prefix+".pdf")
	
	
