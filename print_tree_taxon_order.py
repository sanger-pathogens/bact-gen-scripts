#!/usr/bin/env python
# !/usr/bin/env python
import string
import os, sys
from optparse import OptionParser, OptionGroup

import dendropy


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
	parser.add_option("-t", "--tree", action="store", dest="tree", help="Tree file to allow ordering by clade", default="")
	parser.add_option("-m", "--midpoint", action="store_true", dest="midpoint", help="Midpoint root output tree pdf [default= %default]", default=False, metavar="INT")
	parser.add_option("-l", "--ladderise", action="store", choices=['right', 'left'], dest="ladderise", help="ladderise tree (choose from right or left) [default= %default]", type="choice", default=None)
	
	
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


	if options.tree=='':
		DoError('No tree file selected (-t)')
	elif not os.path.isfile(options.tree):
		DoError('Cannot find file '+options.tree)
	
	
	return



def read_tree(treefile):


	#And to reroot properly if the tree is already rooted on the midpoint
	def reroot_at_midpoint(self, update_splits=False, delete_outdegree_one=True):
		"""
		Reroots the tree at the the mid-point of the longest distance between
		two taxa in a tree.
		Sets the rooted flag on the tree to True.
		If `update_splits` is True, then the edges' `split_bitmask` and the tree's
		`split_edges` attributes will be updated.
		If the *old* root of the tree had an outdegree of 2, then after this
		operation, it will have an outdegree of one. In this case, unless
		`delete_outdegree_one` is False, then it will be
		removed from the tree.
		"""
		from dendropy import treecalc
		pdm = treecalc.PatristicDistanceMatrix(self)
		n1, n2 = pdm.max_dist_nodes
		plen = float(pdm.max_dist) / 2
		mrca_node = pdm.mrca(n1.taxon, n2.taxon)
        #assert mrca_node is self.mrca(taxa=[n1.taxon, n2.taxon])
        #mrca_node = self.mrca(taxa=[n1.taxon, n2.taxon])
		cur_node = n1

		break_on_node = None # populated *iff* midpoint is exactly at an existing node
		target_edge = None
		head_node_edge_len = None

        # going up ...
		while cur_node is not mrca_node:
			if cur_node.edge.length > plen:
				target_edge = cur_node.edge
				head_node_edge_len = plen #cur_node.edge.length - plen
				plen = 0
				break
			elif cur_node.edge.length < plen:
				plen -= cur_node.edge.length
				cur_node = cur_node.parent_node
			else:
				break_on_node = cur_node
				#FIX
				break         #when find the  midpoint, it should break the loop

		assert break_on_node is not None or target_edge is not None

		if break_on_node:
			self.reseed_at(break_on_node, update_splits=False, delete_outdegree_one=delete_outdegree_one)
			new_seed_node = break_on_node
		else:
			tail_node_edge_len = target_edge.length - head_node_edge_len
			old_head_node = target_edge.head_node
			old_tail_node = target_edge.tail_node
			old_tail_node.remove_child(old_head_node)
			new_seed_node = dendropy.dataobject.Node()
			new_seed_node.add_child(old_head_node, edge_length=head_node_edge_len)
			old_tail_node.add_child(new_seed_node, edge_length=tail_node_edge_len)
			self.reseed_at(new_seed_node, update_splits=False, delete_outdegree_one=delete_outdegree_one)
		self.is_rooted = True
		if update_splits:
			self.update_splits(delete_outdegree_one=False)
		return self.seed_node
	
	dendropy.dataobject.Tree.reroot_at_midpoint=reroot_at_midpoint

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
	
	
	#Midpoint root if the option is selected
	if options.midpoint:
		print "Midpoint rooting tree"
		if t.schema=="beast-summary-tree":
			print "Warning: Midpoint rooting a BEAST tree may destroy temporal information represented in node positions"
		t.reroot_at_midpoint(update_splits=True)
		
	#Ladderise the tree if the option is selected
	if options.ladderise=="left":
		print "ladderising tree to the left"
		#ladderise the tree right
		t.ladderize(ascending=False)
	elif options.ladderise=="right":
		print "ladderising tree to the right"
		#ladderise the tree right
		t.ladderize(ascending=True)
	
	for n in t.leaf_iter():
		print n.taxon.label
	
#	print t.as_ascii_plot()
	
	
	return t




################
# Main program #
################		

if __name__ == "__main__":

	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)
	
	tree=read_tree(options.tree)
	
	