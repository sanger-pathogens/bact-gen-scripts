#!/usr/bin/env python

import os, sys
import dendropy, re


def read_tree(treefile):
	dendropy.dataio.nexustokenizer.NHX_COMMENT_FIELD_PATTERN = re.compile(r'(.+?)=({.+?}|.+?)(,|$)')
	dendropy.dataio.nexustokenizer.FIGTREE_COMMENT_FIELD_PATTERN = re.compile(r'(.+?)=({.+?}|.+?)(,|$)')
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
	
	return t
	
def get_leaf_set(treefile):
	leafset=set([])
	for leaf in treefile.leaf_iter():
		leafset.add(leaf.taxon.label)
	return leafset

def get_branch_lengths(treefile, exclude=set([])):
	branch_lengths={}
	for node in treefile.postorder_node_iter():
		if node==treefile.seed_node:
			continue
		downsteam_leaves=set([])
		for dnode in node.postorder_iter():
			if dnode.is_leaf():
				downsteam_leaves.add(dnode.taxon.label)
		brlist=list(downsteam_leaves.difference(exclude))
		brlist.sort()
		brliststring='_'.join(brlist)
				
		branch_lengths[brliststring]=node.edge_length
		
	return branch_lengths

def get_branch_rates(treefile, exclude=set([])):
	branch_rates={}
	for node in treefile.postorder_node_iter():
		if node==treefile.seed_node:
			continue
		downsteam_leaves=set([])
		for dnode in node.postorder_iter():
			if dnode.is_leaf():
				downsteam_leaves.add(dnode.taxon.label)
		brlist=list(downsteam_leaves.difference(exclude))
		brlist.sort()
		brliststring='_'.join(brlist)
		if hasattr(node, "annotations"):
			for a in node.annotations:
				if a.name=='rate':
					branch_rates[brliststring]=a.value
			if not brliststring in branch_rates:
				branch_rates[brliststring]=0.0
		
		
	return branch_rates


t1=read_tree(sys.argv[1])
t2=read_tree(sys.argv[2])

t1_leaf_set=get_leaf_set(t1)
t2_leaf_set=get_leaf_set(t2)

diff_set=t1_leaf_set.symmetric_difference(t2_leaf_set)

bl1=get_branch_lengths(t1, exclude=diff_set)
bl2=get_branch_lengths(t2, exclude=diff_set)
br2=get_branch_rates(t2, exclude=diff_set)

for branch in bl1:
	if branch in bl2:
		print bl1[branch], bl2[branch], br2[branch], branch


	
