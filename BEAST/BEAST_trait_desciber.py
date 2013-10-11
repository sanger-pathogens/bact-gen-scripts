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
#			print count, "trees read"
			sys.stdout.flush()
			
			internal_node_count={}
			
			for node in tree.preorder_node_iter():
#				if not hasattr(node, 'annotations'):
#					continue

#					print dir(node)
#					print node.comments
#					print node.annotations
#					sys.exit()
					
				for x, a in enumerate(node.annotations):
					if isinstance(a.value, str):
						a.value=a.value.replace('"','')
						try:
							node.annotations[x].value=float(a.value)
						except:
							node.annotations[x].value=a.value
					elif isinstance(a.value, list):
						for y in xrange(len(a.value)):
							if isinstance(a.value[y], str):
								a.value[y]=a.value[y].replace('"','')
								node.annotations[x].value[y]=a.value[y]
						
						try:
							node.annotations[x].value=map(float,node.annotations[x].value)
						except:	
							node.annotations[x].value=a.value
				
				#print a.label, a.value, a.name
				annotations={}
				for a in node.annotations:
					annotations[a.name]=a.value
#					print annotations

				if node.parent_node:
				
					for x, a in enumerate(node.parent_node.annotations):
						if isinstance(a.value, str):
							a.value=a.value.replace('"','')
							try:
								node.parent_node.annotations[x].value=float(a.value)
							except:
								node.parent_node.annotations[x].value=a.value
						elif isinstance(a.value, list):
							for y in xrange(len(a.value)):
								if isinstance(a.value[y], str):
									a.value[y]=a.value[y].replace('"','')
									node.parent_node.annotations[x].value[y]=a.value[y]
							
							try:
								node.parent_node.annotations[x].value=map(float,node.parent_node.annotations[x].value)
							except:	
								node.parent_node.annotations[x].value=a.value
					
					#print a.label, a.value, a.name
					parent_annotations={}
					for a in node.parent_node.annotations:
						parent_annotations[a.name]=a.value

						
					
					if options.trait in annotations:
						
						if not annotations[options.trait] in trait_data:
							trait_data[annotations[options.trait]]=[]
						if options.value!="":
							if options.value in annotations:
								try:
									trait_data[annotations[options.trait]].append(float(annotations[options.value]))
								except StandardError:
									print options.value, "is not numeric... aborting"
									sys.exit()
							else:
								print "No metadata named", options.value, " found on node... aborting"
								sys.exit()
						if options.count_changes:
							parent_state=parent_annotations[options.trait]
							
							if not parent_state in internal_node_count:
								internal_node_count[parent_state]=0
							internal_node_count[parent_state]+=1
							
							daughter_state=annotations[options.trait]
							if not parent_state in state_changes:
								state_changes[parent_state]={}
#							if not daughter_state in state_changes:
#								state_changes[daughter_state]={}
							if not daughter_state in state_changes[parent_state]:
								state_changes[parent_state][daughter_state]=[]
#							if not parent_state in state_changes[daughter_state]:
#								state_changes[daughter_state][parent_state]=[]	
							while len(state_changes[parent_state][daughter_state])<(count):
								state_changes[parent_state][daughter_state].append(0)
							
							state_changes[parent_state][daughter_state][-1]+=1
				
				elif options.trait in annotations:
					print "Root Node:"
					print annotations[options.trait]
				
			print "\nInternal Nodes:"
			for internal_state in internal_node_count:
				print internal_state, internal_node_count[internal_state]
					
			#print(tree.as_ascii_plot(plot_metric='length'))
			#print dir(tree)
#			if count>999:
#				break
		
	#	print state_changes, trait_data
		
	
		
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
			print "\nChanges:"
			for from_state in all_states:
				if not from_state in state_changes:
					continue
				for to_state in all_states:
					if not to_state in state_changes[from_state]:
						continue
					print from_state, "->", to_state, state_changes[from_state][to_state][0]
					data.append(state_changes[from_state][to_state])
					labels.append(from_state.replace('"','')+'->'+to_state.replace('"',''))
			
			
			return
			
			pylab.figure()
			n, bins, patches = pylab.hist(data, bins=100, histtype = 'step')
			pylab.xlabel(options.trait)
			pylab.ylabel(options.value)
			pylab.legend(labels)
			if options.output=="" or options.format=="screen":
				pylab.show()
			else:
				pylab.savefig(options.output+"_state_changes_plot."+options.format)
		
			output=open(options.output+'_data.csv', "w")
			outlist=["Tree number"]
			for from_state in all_states:
				for to_state in all_states:
					outlist.append(from_state.replace('"','')+'->'+to_state.replace('"',''))
			
			print >> output, ','.join(outlist)
			for x in xrange(len(data[0])):
				outlist=[str(x+1)]
				for y in xrange(len(data)):
					outlist.append(str(data[y][x]))
				print >> output, ','.join(outlist)
			output.close()
			
				
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
			
			
			
				
				
				
				
				
	