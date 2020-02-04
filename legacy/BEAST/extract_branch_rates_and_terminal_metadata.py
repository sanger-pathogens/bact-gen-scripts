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
	parser.add_option("-m", "--metadata", action="store", dest="m", help="Metadata", default="")
	#Could add more options in here so people can specify similarities etc.
	
	
	return parser.parse_args()




##############################################
# Function toread a tree file using dendropy #
##############################################

def read_dendropy_tree(treefile):
		
		
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
#		for node in t.postorder_node_iter():
#			for x, a in enumerate(node.annotations):
#				if isinstance(a.value, str):
#					print node.annotations[x].value,
#					a.value=a.value.replace('"','')
#					try:
#						node.annotations[x].value=float(a.value)
#						print "here",
#					except:
#						node.annotations[x].value=a.value
#						print "hereb",
#					print node.annotations[x].value
#				elif isinstance(a.value, list):
#					for y in xrange(len(a.value)):
#						if isinstance(a.value[y], str):
#							a.value[y]=a.value[y].replace('"','')
#							node.annotations[x].value[y]=a.value[y]
#					
#					try:
#						node.annotations[x].value=map(float,node.annotations[x].value)
#					except:		
#						if a.name=="!hilight" and len(a.value)==3 and len(a.value[2])>1 and a.value[2][0]=="#":
#							try:
#								rgbint=int(a.value[2][1:])
#							except:
#								break
#							r,g,b=rgbint2rgbtuple(rgbint)
#							node.annotations[x].name="Figtree_hilight"
#							node.annotations[x].value=colors.Color(float(r)/255,float(g)/255,float(b)/255)
#				
#				
#				if isinstance(a.value, str):
#					if a.name=="!color" and len (a.value)>1 and a.value[0]=="#":
#						try:
#							rgbint=int(a.value[1:])
#						except:
#							break
#						r,g,b=rgbint2rgbtuple(rgbint)
#						node.annotations[x].name="Figtree_colour"
#						node.annotations[x].value=colors.Color(float(r)/255,float(g)/255,float(b)/255)
			
		
		return t

def read_metadata(filehandle, column_nums, header=False, split_value="\t"):
	
	metadata={}
	m2h={}
	for x, line in enumerate(filehandle):
		line=line.strip()
		if x==0 and header:
			headings=line.split(split_value)
			if len(headings)<max(column_nums):
				print "Header has too few columns"
				sys.exit()
			for num in column_nums:
				m2h[num]=headings[num-1]
		else:
			columns=line.split(split_value)
			for num in column_nums:
				value=columns[num-1]
				if num in m2h:
					if not columns[0] in metadata:
						metadata[columns[0]]={}
					if not m2h[num] in metadata[columns[0]]:
						metadata[columns[0]][m2h[num]]={}
					metadata[columns[0]][m2h[num]]=value
				else:
					if not columns[0] in metadata:
						metadata[columns[0]]={}
					if not m2h[num] in metadata[columns[0]]:
						metadata[columns[0]][num]={}
					metadata[columns[0]][num]=value
	return metadata

################
# Main program #
################		

if __name__ == "__main__":

	#starttime=time.clock()

	#Get command line arguments

	(options, args) = main()
	
	
	
	mfile=open(options.m)
	metadata=read_metadata(mfile, [4,14], header=True, split_value=",")
	
	
	if options.tree!="":
		if not os.path.isfile(options.tree):
			print "Cannot find file:", options.tree
			sys.exit()
		else:
			tree=read_dendropy_tree(options.tree)
	else:
		print "No tree file specified"
		sys.exit()
	
	mlabels=metadata[metadata.keys()[0]].keys()
	
	print '\t'.join(["Node", "Mean rate", "Median rate", "Length"]+mlabels)
	
	for node in tree.postorder_node_iter():
		meanrate=0.0
		medianrate=0.0
		length=0.0
		
		downstream_meatadata={}
		downstream_taxa=[]
		for leaf in node.leaf_iter():
			if not leaf.taxon.label in metadata:
				print leaf.taxon.label, "not in metadata"
				sys.exit()
			for m in metadata[leaf.taxon.label]:
				if not m in downstream_meatadata:
					downstream_meatadata[m]=set([])
				downstream_meatadata[m].add(metadata[leaf.taxon.label][m])
			downstream_taxa.append(leaf.taxon.label)
		
		do=[]
		for m in mlabels:
			if len(downstream_meatadata[m])==1 and not "" in downstream_meatadata[m]:
				do.append(list(downstream_meatadata[m])[0])
				do.append("unique")
			elif len(downstream_meatadata[m])==1 and "" in downstream_meatadata[m]:
				do.append("unknown")
				do.append("unknown")
			elif len(downstream_meatadata[m])==2 and "" in downstream_meatadata[m]:
				do.append("unique_plus_unknown")
				do.append("unique_plus_unknown")
			elif len(downstream_meatadata[m])>1:
				do.append("multiple")
				do.append("multiple")
		
		for x, a in enumerate(node.annotations):
			if a.name=="rate":
				meanrate=a.value
			if a.name=="rate_median":
				medianrate=a.value
			if a.name=="length":
				length=a.value
		
		
		
		length=node.edge_length
		
		if node.is_leaf():
			nodename=str(node.taxon)
		else:
			nodename=str(node)
		downstream_taxa.sort()
		taxa=','.join(downstream_taxa)
		print '\t'.join([nodename, str(meanrate), str(medianrate), str(length)]+do+[taxa])
		
	
