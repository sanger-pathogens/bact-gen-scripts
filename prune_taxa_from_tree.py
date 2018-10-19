#!/usr/bin/env python3

#################################
# Import some necessary modules #
#################################

import os, sys
import dendropy
from optparse import OptionParser, OptionGroup
from random import randrange
from scipy.stats import chisquare
import numpy as np



##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
        
	group = OptionGroup(parser, "Required options")
	group.add_option("-t", "--tree", action="store", dest="tree", help="Input tree file in newick or nexus format", default="", metavar="FILE")
	group.add_option("-f", "--sample_file", action="store", dest="samples", help="File containing list of samples to include/exclude (e.g. csv or tsv)", default="", metavar="FILE")
	group.add_option("-o", "--output_file", action="store", dest="output", help="Output file name", default="pruned_tree.nwk")
	group.add_option("-c", "--column", action="store", dest="column", help="Number of the column in the sample file that contains the sample name [default= %default]", default=1, type="int")
	group.add_option("-s", "--separator", action="store", dest="separator", help="Separator for columns in input file [default= %default]", default=",")
	group.add_option("-k", "--keep", action="store_true", dest="keep", help="Keep specified samples [default=remove]", default=False)
	parser.add_option_group(group)
	group = OptionGroup(parser, "Other options")
	group.add_option("-v", "--verbose", action="store_true", dest="verbose", help="Verbose [default= %default]", default=False)
	parser.add_option_group(group)
	
	return parser.parse_args()
	

def check_command_line_options(o):
	
	if  not os.path.isfile(o.tree):
		print("Cannot find file", o.tree)
		sys.exit()
		
	if o.samples!="" and not os.path.isfile(o.samples):
		print("Cannot find file", o.samples)
		sys.exit()
	
	if o.column<1:
		print("Column number must be greater than 0")
		sys.exit()
	return


def read_tree(tfile):
	try:
		t=dendropy.Tree.get(path=tfile, schema="nexus", preserve_underscores=True)
	except dendropy.dataio.nexusreader.NexusReader.NotNexusFileError:
		t=dendropy.Tree.get(path=tfile, schema="newick", preserve_underscores=True)
	except:
		print("Cannot read tree file. It should be in newick or nexus format")
	return t

	
################
# Main program #
################                

if __name__ == "__main__":

	(options, args) = main()
	check_command_line_options(options)
	
	treefile=options.tree
		
	samplefile=options.samples
	
	output_file=options.output
	
	column=options.column-1
	
	separator=options.separator
	
	verbose=options.verbose
	
	keep=options.keep
	
	tree=read_tree(treefile)
	
	tree_taxa=list(tree.taxon_namespace.labels())
	
	samples=set([])
	not_in_tree=[]
	if samplefile!="":
		for line in open(samplefile, "rU", encoding='ISO-8859-1'):
			words=line.strip().split(separator)
			if len(words)<column:
				if options.verbose:
					print("Invalid line in samples file (Too few columns):", line.strip())
			if not words[0] in tree_taxa:
				if options.verbose:
					print(words[0], "not in tree. Skipping")
				not_in_tree.append(words[0])
				continue
			if words[column]=="":
				if options.verbose:
					print(' '.join(words), "has missing sample name")
				continue
			samples.add(words[column])
	
	if len(args)>0:
		for arg in args:
			if len(arg)>0:
				words=arg.split(",")
				for word in words:
					if not word in tree_taxa:
						if options.verbose:
							print(word, "not in tree. Skipping")
						not_in_tree.append(word)
						continue
					samples.add(word)
	
	if keep:
		samples_to_remove=list(samples)
	else:
		samples_to_remove=list(samples)
	
	if len(not_in_tree)>0:
		print(len(not_in_tree), "samples not in tree")
	
	print(len(samples_to_remove), "samples(s) will be removed from tree")
	pruned_tree=tree.extract_tree_without_taxa_labels(labels=samples_to_remove)
	pruned_tree.write(path=output_file, schema="newick", unquoted_underscores=True)
	
