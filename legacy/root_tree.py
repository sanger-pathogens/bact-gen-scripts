#!/usr/bin/env python

import os, sys
sys.path.extend(map(os.path.abspath, ['/nfs/pathogen/sh16_scripts/modules/']))
from Bio.Nexus import Trees, Nodes
import Si_nexus
from Si_general import *
from optparse import OptionParser


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	
	parser.add_option("-r", "--root", action="store", dest="outgroup", help="taxon to root on (if not specified, will midpoint root)", default="")
	parser.add_option("-t", "--tree", action="store", dest="tree", help="tree file", default="")
	parser.add_option("-l", "--ladderise", action="store", dest="ladderise", help="ladderise tree (choose from left and right)", type="choice", choices=["right", "left"], default=None)
	parser.add_option("-o", "--output", action="store", dest="outputfile", help="output file name", default="")
	
	return parser.parse_args()


################
# Main program #
################		

if __name__ == "__main__":


	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	if options.tree=="":
		DoError("No tree file specified")


	print "Reading tree file"
	sys.stdout.flush()
	
	try:
		tree_strings = open(options.tree).readlines()
	except IOError:
		DoError("Cannot open tree file "+options.tree)
	handle = open(options.outputfile, "w")
	for tree_string in tree_strings:
		
		tree = Trees.Tree(tree_string, rooted=True)
		
		
		if options.outgroup!="":
			print "Rooting tree on", options.outgroup
			sys.stdout.flush()
			tree.root_with_outgroup(outgroup=options.outgroup)
		else:
			print "Midpoint rooting tree"
			tree=Si_nexus.midpoint_root(tree)
	
		
		
		#Print tree
		
		treestring=Si_nexus.tree_to_string(tree,plain=False,ladderize=options.ladderise)#tree with branch lengths and support
	
		if treestring[-1]!=";":
			treestring=treestring+";"
				
		
		
		print >> handle, treestring
	handle.close()
