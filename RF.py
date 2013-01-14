#!/usr/bin/env python

#################################
# Import some necessary modules #
#################################
import string, re
import os, sys
import dendropy
from optparse import OptionParser



##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "%prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-t", "--trees", action="store", dest="trees", help="file containing trees to compare", default="")

	
	return parser.parse_args()

################
# Main program #
################		

if __name__ == "__main__":

	#Get command line arguments
	
	(options, args) = main()
	
	trees = dendropy.TreeList.get_from_path(options.trees, "newick")
	print '\t'.join(['tree1', 'tree2','symmetric distance', 'false positives and negatives', 'euclidean distance', 'robinson foulds distance'])
	for treeonenum in xrange(0,len(trees)-1):
		treeone=trees[treeonenum]
		for treetwonum in xrange(treeonenum+1,len(trees)):
			treetwo=trees[treetwonum]
			print '\t'.join(map(str,[treeonenum+1, treetwonum+1, treeone.symmetric_difference(treetwo), treeone.false_positives_and_negatives(treetwo), treeone.euclidean_distance(treetwo), treeone.robinson_foulds_distance(treetwo)]))
		