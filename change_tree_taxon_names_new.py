#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math
from optparse import OptionParser, OptionGroup
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_general import *
import dendropy
	
##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-t", "--tree", action="store", dest="tree", help="Input tree file", default="", metavar="FILE")
	parser.add_option("-i", "--info", action="store", dest="info", help="Input information (csv) file", default="", metavar="FILE")
	parser.add_option("-c", "--columns", action="store", dest="columnsstring", help="Column numbers to include in the new names", default="", metavar="LIST")
	parser.add_option("-o", "--output", action="store", dest="outfile", help="Name for output tree file", default="")
	parser.add_option("-T", "--tab", action="store_true", dest="tab", help="csv file is tab separated (default is comma separated)", default=False)
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="More verbose messages", default=False)

	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):
	
	if len(options.columnsstring)==0:
		DoError('You need to choose at least one column (-c option)')
	
	try: options.columns=map(int,options.columnsstring.split(","))
	except ValueError:
		DoError('column numbers must be integers separated by commas')
	
	if options.tree=='':
		DoError('No tree file selected')
	elif not os.path.isfile(options.tree):
		DoError('Cannot find file '+options.tree)
	
	if options.info=='':
		DoError('No tree file selected')
	elif not os.path.isfile(options.info):
		DoError('Cannot find file '+options.info)

	
	
	if len(options.outfile)==0:
		DoError('You need to give an output tree file name (-o option)')



################
# Main program #
################		

if __name__ == "__main__":

	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)

	#treedata=open(options.tree,'rU').read()
	
	
	namedict={}
	for line in open(options.info,'rU'):
		newname=""
		if options.tab:
			words=line.strip().split('\t')
		else:
			words=line.strip().split(',')
			
		if len(words)==0:
			replace=False
			continue
		count=1
		for column in options.columns:
			try:
				x=int(column)
			except Standarderror:
				continue
			if len(words)<int(column):
				continue
			elif len(words[int(column)-1])==0:
				continue
			if count==1:
				newname=words[int(column)-1].replace(' ','_').replace(':','_').replace('(', '_').replace(')','_')
				count=count+1
			else:
				newname=newname+'_'+words[int(column)-1].replace(' ','_').replace(':','_').replace('(', '_').replace(')','_')
		
		if len(newname)>0:
			if newname[-1]=='_':
				newname=newname[:-1]
			
			namedict[words[0].replace("_"," ")]=newname
	
	
	if not os.path.isfile(options.tree):
		print "Cannot find file:", options.tree
		options.tree=""
	else:
		opened=False
		for treeschema in ["beast-summary-tree", "nexus", "newick"]:
			try:
				t = dendropy.Tree.get_from_path(options.tree, schema=treeschema)
				opened=True
				break
			except dendropy.utility.error.DataParseError:
				continue
		if not opened:
			print "Failed to open tree file"
			sys.exit()
		
		
		for leaf in t.leaf_iter():
			if leaf.taxon.label in namedict:
				leaf.taxon.label=namedict[leaf.taxon.label]
			elif options.verbose:
				print "Cannot find "+leaf.taxon.label+" in metadata file"
			
		
		if treeschema=="beast-summary-tree":
			treeschema="nexus"
		
		
		t.write_to_path(options.outfile, treeschema)
				

        
        
        
        
	
