#!/usr/bin/env python
import string, re, gzip
import os, sys
from Bio.SeqRecord import SeqRecord
from Bio.Nexus import Trees, Nodes
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_nexus import draw_ascii_tree, tree_to_string, midpoint_root
from optparse import OptionParser



##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "%prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-t", "--tree", action="store", dest="tree", help="starting tree (optional)", default="")
	parser.add_option("-c", "--cutoff", action="store", dest="cutoff", type="float", help="Cutoff value for collapsing branches [default=%default]", default=0.0)
	parser.add_option("-o", "--outfile", action="store", dest="outfile", help="Output file name", default="")

	
	return parser.parse_args()


################
# Main program #
################		

if __name__ == "__main__":

	#Get command line arguments
	
	(options, args) = main()

	try:
		tree_string = open(options.tree).read()
	except IOError:
		DoError("Cannot open tree file "+options.tree)
	tree = Trees.Tree(tree_string, rooted=True)
	
	treestring=tree_to_string(tree, False, True, True, True,collapse=True, cutoff=options.cutoff)
	#(treeObject, support_as_branchlengths=False,branchlengths_only=False,plain=True,plain_newick=False,ladderize=None)
	
	handle = open(options.outfile, "w")
	print >> handle, treestring+";"
	handle.close()