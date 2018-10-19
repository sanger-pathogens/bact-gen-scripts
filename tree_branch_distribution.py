#!/usr/bin/env python
import string, re, gzip
import os, sys
from Bio.SeqRecord import SeqRecord
from Bio.Nexus import Trees, Nodes
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_nexus import draw_ascii_tree, tree_to_string, midpoint_root
from optparse import OptionParser

def DoError(Errorstring):
	print "Error:", Errorstring
	print "Use -h for help"
	sys.exit()

##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "%prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-t", "--tree", action="store", dest="tree", help="starting tree (optional)", default="")
	parser.add_option("-c", "--cutoff", action="store", dest="cutoff", type="float", help="Cutoff value for collapsing branches [default=%default]", default=0.0)
	parser.add_option("-m", "--min", action="store_true", dest="min", help="Identify the minimum branch length form the tree and use that as the cutoff. Will overrider -c option.", default=False)
	
	
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
	
	if options.min:
		if options.cutoff>0:
			print "Warning: Cutoff will be ignored as you also chose the minimum branch length cutoff option"
	else:
		if options.cutoff==0:
			print "Warning: Cutoff value of zero chosen. Only negative branch lengths will be collapsed."
		elif options.cutoff<0:
			print "Warning: negative cutoff chosen. Only negative branch lengths will be collapsed."
	
	if options.min:
		root=tree.root
		brlens=[]
		for n in tree._walk(root):
			if tree.node(n).data.branchlength:
				brlens.append(tree.node(n).data.branchlength)
		brlens.sort()
		if len(brlens)==0:
			DoError("Found no branch lengths in tree")
		
		curbrlen=float("Inf")
		curcount=0
		print "\nDistribution of branch lengths (*=one branch):\n"
		for br in brlens:
			if curbrlen!=float("Inf") and br!=curbrlen:
				print curbrlen, "*"*curcount, curcount
				curcount=0
			curbrlen=br
			curcount+=1
		print
		print type(brlens[0]), brlens[0], brlens[:1]
		minbrlen=brlens[0]
		nextbrlen=float("Inf")
		for br in brlens:
			if br>minbrlen:
				nextbrlen=br
				break
		print "Minimum branch length found =", minbrlen
		print "Next smallest branch length found =", nextbrlen
		options.cutoff=brlens[0]+((nextbrlen-brlens[0])/2)
	
	
