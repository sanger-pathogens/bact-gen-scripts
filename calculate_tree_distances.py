#!/usr/bin/env python

#################################
# Import some necessary modules #
#################################

import dendropy
import string, re
import os, sys

##############################################
# Function toread a tree file using dendropy #
##############################################

def read_dendropy_tree(treefile):
		
		#Some fixes to dendropy to make it read annotations properly
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
				treelist = dendropy.TreeList.get_from_path(treefile, schema=treeschema, as_rooted=False, preserve_underscores=True, case_insensitive_taxon_labels=False, set_node_attributes=True, extract_comment_metadata=True)
				opened=True
				treelist.schema=treeschema
				break
			except dendropy.utility.error.DataParseError:
				continue
			except ValueError as e:
				print "Encountered ValueError while trying to read tree file as", treeschema+":", e
				continue
				
		if not opened:
			print "Failed to open tree file"
			sys.exit()

		
		return treelist

if len(sys.argv)!=2:
	print "Usage: ~sh16/scripts/calculate_tree_distances.py <tree file>"
	sys.exit()

treelist=read_dendropy_tree(sys.argv[1])

print '\t'.join(["Tree 1", "Tree 2", "Robinson Foulds Distance", "Weighted Robinson Foulds Distance", "Euclidean Distance"])
for x, t1 in enumerate(treelist):
	
	for y in xrange(x+1, len(treelist)):
		
		t2=treelist[y]
		print '\t'.join(map(str,[x+1, y+1, dendropy.treecalc.symmetric_difference(t1,t2), dendropy.treecalc.robinson_foulds_distance(t1,t2), dendropy.treecalc.euclidean_distance(t1,t2)]))
	
