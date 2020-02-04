#!/usr/bin/env python

import os, sys, dendropy
from optparse import OptionParser

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
	
	parser.add_option("-t", "--tree", action="store", dest="tree", help="Input tree file", default="", metavar="FILE")
	parser.add_option("-l", "--list", action="store", dest="list", help="File containing list of taxa to remove or keep", default="")
	parser.add_option("-k", "--keep", action="store_true", dest="keep", help="Keep only taxa in list (default is to remove them)", default=False)
	parser.add_option("-o", "--output", action="store", dest="output", help="Output file name", default="")
	
	#parser.add_option("-M", "--midpoint", action="store_true", dest="midpoint", help="Midpoint root output tree pdf [default= %default]", default=False, metavar="INT")
	#parser.add_option("-L", "--ladderise", action="store", choices=['right', 'left'], dest="ladderise", help="ladderise tree (choose from right or left) [default= %default]", type="choice", default=None)
	


	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


	if options.tree=='':
		DoError('No tree file selected')
	elif not os.path.isfile(options.tree):
		DoError('Cannot find file '+options.tree)
	if options.list=='':
		DoError('No list file selected')
	elif not os.path.isfile(options.list):
		DoError('Cannot find file '+options.list)	
	if options.output=='':
		DoError('No output file name selected')
	
	return

##################
# Read tree file #
##################

def read_tree(treefile):
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
	
	#t.deroot()
	
	
	return t

################
# Main program #
################		

if __name__ == "__main__":

	#Get command line arguments

	(options, args) = main()
	
	#Do some checking of the input files
	
	check_input_validity(options, args)
	
	warncount=0
	
	try:
		tree=read_tree(options.tree)
	except StandardError:
		DoError("Failed to read tree file")
	
	tree_taxa=set([])
	for taxon in tree.taxon_set:
		tree_taxa.add(taxon.label)
	
	print "Reading taxon list"
	try:
		taxa_to_remove=set([])
		for line in open(options.list, "rU"):
			if len(line)==0:
				continue
			taxa_to_remove.add(line.strip().split()[0])
	except StandardError:
		DoError("Failed to read list file")
	
	if options.keep:
		print "Found", len(taxa_to_remove), "taxa to keep"
	else:
		print "Found", len(taxa_to_remove), "taxa to remove"
	
	
	difference=taxa_to_remove.difference(tree_taxa)
	if len(difference)>0:
		if len(difference)==1:
			print "Warning:", len(difference), "taxon in list is not in tree:", ", ".join(list(difference))
		else:
			print "Warning:", len(difference), "taxa in list are not in tree", ", ".join(list(difference))
		warncount+=len(difference)
		taxa_to_remove.intersection_update(set(tree_taxa))
		if options.keep:
			print "Found", len(taxa_to_remove), "taxa matching tree to keep"
		else:
			print "Found", len(taxa_to_remove), "taxa matching tree to remove"
		
		
	if options.keep:
		if len(taxa_to_remove)<3:
			print "Error: Too few taxa ("+str(len(taxa_to_remove))+") remaining in tree (must be at least 3)"
			sys.exit()
		tree.retain_taxa_with_labels(taxa_to_remove)
	else:
		if len(tree_taxa)-len(taxa_to_remove)<3:
			print "Error: Too few taxa ("+str(len(tree_taxa)-len(taxa_to_remove))+") remaining in tree (must be at least 3)"
			sys.exit()
		tree.prune_taxa_with_labels(taxa_to_remove)
	
	
	
	
	output=open(options.output, "w")
	print>> output, tree.as_string('newick')
	output.close()
	
	if warncount>0:
		print "Completed with", warncount, "warning(s)"
	else:
		print "Successfully completed"

	
