#!/usr/bin/env python

#################################
# Import some necessary modules #
#################################

import dendropy
import pylab
import string, re
import os, sys
#import random
#from math import sqrt, pow, log
#from numpy import repeat, convolve, mean, median
from numpy import mean, median
from optparse import OptionParser, OptionGroup
#import shlex, subprocess
#on my laptop
#sys.path.extend(map(os.path.abspath, ['/Users/sh16/Documents/scripts/modules/']))
#on pcs4
#sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
#from Si_general import *
#from Si_SeqIO import *
from math import floor

################################
# Get the command line options #
################################


def main():
	usage = "usage: %prog [options] args"
	parser = OptionParser(usage=usage)
	
	
	
	parser.add_option("-t", "--tree", action="store", dest="tree", help="tree file to align tab files to", default="")
	parser.add_option("-T", "--trait", action="store", dest="trait", help="name of trait", default="")
	parser.add_option("-v", "--value", action="store", dest="value", help="name of value to store for metadata", default="")
	parser.add_option("-o", "--output", action="store", dest="output", help="prefix for output files", default="")
	parser.add_option("-f", "--format", action="store", dest="format", help="format for output plot files. Choose from pdf, png or screen (show on screen rather than save to file) [default=%default]", default="pdf", type="choice", choices=["screen", "pdf", "png"])
	parser.add_option("-b", "--burnin", action="store", dest="burnin", help="number of trees to remove from the start of the trees file(s) [default=%default]", default=0, type='int')
	parser.add_option("-s", "--subsample", action="store", dest="subsample", help="subsample trees every x trees [default=%default]", default=1, type='int')
	parser.add_option("-c", "--count_changes", action="store_true", dest="count_changes", help="count changes between states", default=False)
	
	
	return parser.parse_args()



###############################################
# Check the command line options are sensible #
###############################################





def read_dendropy_tree(treefile):

#		opened=False
#		for treeschema in ["beast-summary-tree", "nexus", "newick"]:
#			try:
#				trees = dendropy.TreeList.get_from_path(treefile, schema=treeschema, as_rooted=True)
#				opened=True
#				break
#			except dendropy.utility.error.DataParseError:
#				continue
#		if not opened:
#			print "Failed to open tree file"
#			sys.exit()
#		
#		for tree in trees:
#			print(tree.as_ascii_plot(plot_metric='length'))
#			sys.exit()
		if options.count_changes:
			state_changes={}
		trait_data={}
		count=0
		subsamplecount=1
		
		if options.subsample<1:
			print "subsample must be one or greater"
			sys.exit()
		
#		print dir(dendropy.tree_source_iter), dendropy.tree_source_iter.func_globals
#		sys.exit()
		
		for tree in dendropy.tree_source_iter(stream=open(treefile, "rU"), schema="nexus",  extract_comment_metadata=True, tree_offset=options.burnin):#ignore_missing_node_info=True,
			
			subsamplecount-=1
			if subsamplecount>0:
				continue
			else:
				subsamplecount=options.subsample
			count+=1
			print count
			sys.stdout.flush()
			for node in tree.preorder_node_iter():
				if node.parent_node:
					if options.trait in node.comment_metadata:
						if not node.comment_metadata[options.trait] in trait_data:
							trait_data[node.comment_metadata[options.trait]]=[]
						if options.value!="":
							if options.value in node.comment_metadata:
								try:
									trait_data[node.comment_metadata[options.trait]].append(float(node.comment_metadata[options.value]))
								except StandardError:
									print options.value, "is not numeric... aborting"
									sys.exit()
							else:
								print "No metadata named", options.value, " found on node... aborting"
								sys.exit()
						if options.count_changes:
							parent_state=node.parent_node.comment_metadata[options.trait]
							daughter_state=node.comment_metadata[options.trait]
							if not parent_state in state_changes:
								state_changes[parent_state]={}
							if not daughter_state in state_changes[parent_state]:
								state_changes[parent_state][daughter_state]=[]
							while len(state_changes[parent_state][daughter_state])<(count):
								state_changes[parent_state][daughter_state].append(0)
							
							state_changes[parent_state][daughter_state][-1]+=1

					

			#print(tree.as_ascii_plot(plot_metric='length'))
			#print dir(tree)
#			if count>999:
#				break
		
		if count==0:
			print "No trees visited. Perhaps you subsampled too much?"
			sys.exit()
		
		if options.value!="":
			data=[]
			labels=[]
			for trait in trait_data:
				labels.append(trait.replace('"',''))
				print trait, mean(trait_data[trait])
				data.append(trait_data[trait])
			
			pylab.figure()
			pylab.boxplot(data,0,'')
			pylab.xticks([1,2],labels, rotation=45)
			pylab.xlabel(options.trait)
			pylab.ylabel(options.value)
	#		title('Chickenpox cases in NYC 1931-1971')
			
			if options.output=="" or options.format=="screen":
				pylab.show()
			else:
				pylab.savefig(options.output+"_values_boxplot."+options.format)
		
		if options.count_changes:
			#print state_changes
			labels=[]
			data=[]
			
			all_states=[]
			
			for from_state in state_changes:
				if not from_state in all_states:
					all_states.append(from_state)
				for to_state in state_changes[from_state]:
					if not to_state in all_states:
						all_states.append(to_state)
			
			for from_state in all_states:
				for to_state in all_states:
					data.append(state_changes[from_state][to_state])
					labels.append(from_state.replace('"','')+'->'+to_state.replace('"',''))
			pylab.figure()
			n, bins, patches = pylab.hist(data, bins=100, histtype = 'step')
			pylab.xlabel(options.trait)
			pylab.ylabel(options.value)
			pylab.legend(labels)
			if options.output=="" or options.format=="screen":
				pylab.show()
			else:
				pylab.savefig(options.output+"_state_changes_plot."+options.format)
		
		

			
				
		sys.exit()




################
# Main program #
################		

if __name__ == "__main__":

	#starttime=time.clock()

	#Get command line arguments

	(options, args) = main()
	
	if options.tree!="":
		if not os.path.isfile(options.tree):
			print "Cannot find file:", options.tree
			sys.exit()
		else:
			read_dendropy_tree(options.tree)
			
			
			
				
				
				
				
				
	