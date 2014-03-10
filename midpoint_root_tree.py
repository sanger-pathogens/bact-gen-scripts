#!/usr/bin/env python

#################################
# Import some necessary modules #
#################################

import dendropy
import string, re
import os, sys
from optparse import OptionParser


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options] <list of mfa files>"
	parser = OptionParser(usage=usage)
	parser.add_option("-t", "--tree", action="store", dest="tree", help="Input tree file", default="")
	parser.add_option("-o", "--output", action="store", dest="output", help="Output tree file name", default="")
	parser.add_option("-f", "--format", choices=["nexus", "newick"], type="choice", action="store", dest="out_format", help="Output format "+", ".join(["nexus", "newick"])+" [default= %default]", default="newick")
	parser.add_option("-l", "--ladderise", choices=["left", "right", "none"], type="choice", action="store", dest="ladderise", help="Laddersie tree "+", ".join(["nexus", "newick"])+" [default= %default]", default="none")
	#Could add more options in here so people can specify similarities etc.
	
	
	return parser.parse_args()




##############################################
# Function toread a tree file using dendropy #
##############################################

def read_dendropy_tree(treefile):
		
			
		def _parse_taxlabels_statement(self, taxon_set=None):
			"""
			Processes a TAXLABELS command. Assumes that the file reader is
			positioned right after the "TAXLABELS" token in a TAXLABELS command.
			"""
			if taxon_set is None:
				taxon_set = self._get_taxon_set()
			token = self.stream_tokenizer.read_next_token()
			while token != ';':
				label = token
				if len(taxon_set) >= self.file_specified_ntax and not self.attached_taxon_set:
					raise self.too_many_taxa_error(taxon_set=taxon_set, label=label)
				taxon = taxon_set.require_taxon(label=label)
				self.stream_tokenizer.store_comment_metadata(taxon)
				self.stream_tokenizer.store_comments(taxon)
				token = self.stream_tokenizer.read_next_token()
				
		
		dendropy.dataio.nexusreader_py.NexusReader._parse_taxlabels_statement=_parse_taxlabels_statement
		

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
			if t.schema=="nexus":
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
		
		#Make sure the tree is rooted on an edge rather than a node
		if t.is_unrooted:
			print "Tree is unrooted. Rooting it now."
			t.is_rooted=True
			t.update_splits(delete_outdegree_one=False)
		root = t.seed_node
		root_children = root.child_nodes()
		
		if len(root_children) != 2:
			print "Tree rooted at node. Rerooting on first edge from that node."
			t.reroot_at_edge(root_children[0].edge, update_splits=True)
		#print the tree in ascii as a cladogram
		#print(t.as_ascii_plot())
		#print the tree in ascii including branch lengths
		#print(t.as_ascii_plot(plot_metric='length'))
		
		return t
			

################
# Main program #
################		

if __name__ == "__main__":

	#starttime=time.clock()

	#Get command line arguments

	(options, args) = main()
	options.midpoint=True
	
	if options.output=="":	
		print " No output file specified"
		sys.exit()
	
	if options.tree!="":
		if not os.path.isfile(options.tree):
			print "Cannot find file:", options.tree
			sys.exit()
		else:
			tree=read_dendropy_tree(options.tree)
	else:
		print "No tree file specified"
		sys.exit()
	
	print "Printing rooted tree in", options.out_format, "format"
	treestring=tree.as_string(options.out_format).replace("[&R] ","")
	outfile=open(options.output, "w")
	print >> outfile, treestring
	outfile.close()
	print "Done"
	
