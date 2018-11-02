#!/usr/bin/env python
# !/usr/bin/env python
import string
import os, sys
from Bio import Seq
from Bio import AlignIO
from Bio.Align import Generic
from Bio.Alphabet import IUPAC, Gapped
from optparse import OptionParser, OptionGroup
import glob
import dendropy
import shlex, subprocess

sys.path.extend(map(os.path.abspath, ['/nfs/pathogen/sh16_scripts/modules/']))
from Si_general import *
from Si_SeqIO import *



####################
# Set some globals #
####################


MUMMER_DIR=""


#############################################
# Function to reverse complement a sequence #
#############################################

def revcomp(sequence):
        rev=sequence[::-1]
        revcomp=''
        d={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'g':'c', 'c':'g'}
        for i in rev:
                if d.has_key(i):
                        revcomp=revcomp+d[i]
                else:
                        revcomp=revcomp+i
        
        return revcomp


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	
	parser.add_option("-c", "--contigs", action="store", dest="contigs", help="Contig fasta file name", default="", metavar="FILE")
	parser.add_option("-o", "--output", action="store", dest="output", help="Output file name", default="", metavar="STRING")
	parser.add_option("-t", "--tree", action="store", dest="tree", help="Tree file to allow ordering by clade", default="")
	parser.add_option("-m", "--midpoint", action="store_true", dest="midpoint", help="Midpoint root output tree pdf [default= %default]", default=False, metavar="INT")
	parser.add_option("-l", "--ladderise", action="store", choices=['right', 'left'], dest="ladderise", help="ladderise tree (choose from right or left) [default= %default]", type="choice", default=None)
	parser.add_option("-s", "--sortby", action="store", choices=['tree', 'location'], dest="sortby", help="Sort by which method first (choose from tree or location) [default= %default]", type="choice", default='location')
	
	
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


	if options.contigs=='':
		DoError('No reference fasta file selected (-r)')
	elif not os.path.isfile(options.contigs):
		DoError('Cannot find file '+options.contigs)
	
	
	if options.output=="":
		DoError('No output prefix selected (-o)')
	
	
	return



def read_tree(treefile):
	#Try opening the tree using various schemas

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
#
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
	
#	dendropy.dataobject.Tree.reroot_at_midpoint=reroot_at_midpoint

	#Try opening the tree using various schemas
	opened=False
	for treeschema in ["nexus", "newick"]:#["beast-summary-tree",  "nexus", "newick"]:
		try: 
			t = dendropy.Tree.get_from_path(treefile, schema=treeschema, rooting='force-rooted', preserve_underscores=True, case_sensitive_taxon_labels=True, extract_comment_metadata=True)
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
	
	leaf_order=[]
	
	for n in t.leaf_node_iter():
		leaf_order.append(n.taxon.label)
	
#	print t.as_ascii_plot()
	
	
	return leaf_order





################
# Main program #
################		

if __name__ == "__main__":

	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)
	
	if options.tree!="":
		tree_order=read_tree(options.tree)
	else:
		tree_order=[]
	tree_positions={}
	for x, t in enumerate(tree_order):
		tree_positions[t]=x+1
#	
#	sys.exit()
	
	
	try:
		refseqs=read_seq_file(options.contigs)
	except StandardError:
		DoError("Cannot open reference file")
	
	reflist=[]
	reflens=[]
	refstarts={}
	refends={}
	lens={}
	seqs={}
	totlen=0
	for seq in refseqs:
		if seq.id in reflist:
			DoError("You have two identically named sequences in your reference fasta file")
		reflist.append(seq.id)
		reflens.append(totlen)
		refstarts[seq.id]=totlen
		refends[seq.id]=totlen+len(seq.seq)
		seqs[seq.id]=str(seq.seq).upper()
		totlen+=len(seq.seq)
		lens[seq.id]=len(seq.seq)

		try:
			refstarts[seq.id]=int(seq.description.split("#")[1].split("..")[0].strip())+1
		except ValueError:
			refstarts[seq.id]=float("Inf")
		except IndexError:
				refstarts[seq.id]=float("Inf")
		
		if refstarts[seq.id]==float("Inf") or refstarts[seq.id]=="":
			try:
				refstarts[seq.id]=int(seq.description.split("#")[2].split("..")[0].strip())+1
			except ValueError:
				refstarts[seq.id]=float("Inf")
			except IndexError:
					refstarts[seq.id]=float("Inf")
		
		if refstarts[seq.id]==float("Inf") or refstarts[seq.id]=="":
			try:
				refstarts[seq.id]=int(seq.description.split("#")[1].split("..")[1].strip())+1
			except ValueError:
				refstarts[seq.id]=float("Inf")
			except IndexError:
				refstarts[seq.id]=float("Inf")
		
		if refstarts[seq.id]==float("Inf") or refstarts[seq.id]=="":
			try:
				refstarts[seq.id]=int(seq.description.split("#")[2].split("..")[1].strip())+1
			except ValueError:
				refstarts[seq.id]=float("Inf")
			except IndexError:
				refstarts[seq.id]=float("Inf")
	
	reforders=[]
	for ref in refstarts:
		if ref.split("$")[0] in tree_positions:
			tree_pos=tree_positions[ref.split("$")[0]]
		else:
			tree_pos=float("Inf")
		if options.sortby=="tree":
			reforders.append([tree_pos, refstarts[ref], ref])
		else:
			reforders.append([refstarts[ref], tree_pos, ref])
			
	
	reforders.sort()
	output=open(options.output, "w")
	for ref in reforders:
		print ref[0], ref[1], ref[2]
		print >> output, ">"+ref[2]
		print >> output, seqs[ref[2]]
	
	
	sys.exit()
	
