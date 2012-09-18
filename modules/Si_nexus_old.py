from Bio.Nexus import Trees, Nodes
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, SeqFeature
import sys, string, os, glob, copy
from Si_general import *
from random import *
import Si_SNPs_temp
import math
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram._Colors import ColorTranslator

RAXML_DIR="/software/pathogen/external/applications/RAxML/RAxML-7.0.4/"


ambiguity_to_bases={"A":["A"], "C":["C"], "T":["T"], "G":["G"], "M":["A", "C"], "K":["G", "T"], "R":["A", "G"], "Y":["C", "T"], "S":["C", "G"], "W":["A", "T"], "B":["C", "G", "T"], "V":["A", "C", "G"], "H":["A", "C", "T"], "D":["A", "G", "T"], "?":["-", "A", "C", "G", "T"], "N":["-", "A", "C", "G", "T"], "X":["-", "A", "C", "G", "T"], "-":["-"], "a":["-", "A"], "c":["-","C"], "t":["-","T"], "g":["-","G"], "m":["-","A","C"], "k":["-", "G", "T"], "r":["-", "A","G"], "y":["-","C","T"], "s":["-","C","G"], "w":["-","A","T"], "b":["-","C","G","T",], "v":["-","A","C","G"], "h":["-","A","C","T"], "d":["-","A","G","T",], "n":["-","A","C","G","T"]}

bases_to_ambiguity={"A":"A", "C":"C", "T":"T", "G":"G", "AC":"M", "GT":"K", "AG":"R", "CT":"Y", "CG":"S", "AT":"W", "CGT":"B", "ACG":"V", "ACT":"H", "AGT":"D", "ACGT":"N", "-":"-", "?":"N", "X":"N", "-A":"a", "-C":"c", "-T":"t", "-G":"g", "-AC":"m", "-GT":"k", "-AG":"r", "-CT":"y", "-CG":"s", "-AT":"w", "-CGT":"b", "-ACG":"v", "-ACT":"h", "-AGT":"d", "-ACGT":"N"}


transition_matrix={"A":{"A":0, "C":1, "G":1, "T":1, "-":1}, "C":{"A":1, "C":0, "G":1, "T":1, "-":1}, "G":{"A":1, "C":1, "G":0, "T":1, "-":1}, "T":{"A":1, "C":1, "G":1, "T":0, "-":1}, "-":{"A":1, "C":1, "G":1, "T":1, "-":0}}







def RGBToHTMLColor(rgb_tuple):
    """ convert an (R, G, B) tuple to #RRGGBB """
    hexcolor = '#%02x%02x%02x' % rgb_tuple
    # that's it! '%02x' means zero-padded, 2-digit hex values
    return hexcolor

def HTMLColorToRGB(colorstring):
    """ convert #RRGGBB to an (R, G, B) tuple """
    colorstring = colorstring.strip()
    if colorstring[0] == '#': colorstring = colorstring[1:]
    if len(colorstring) != 6:
        raise ValueError, "input #%s is not in #RRGGBB format" % colorstring
    r, g, b = colorstring[:2], colorstring[2:4], colorstring[4:]
    r, g, b = [int(n, 16) for n in (r, g, b)]
    return (r, g, b)

def HTMLColorToPILColor(colorstring):
    """ converts #RRGGBB to PIL-compatible integers"""
    colorstring = colorstring.strip()
    while colorstring[0] == '#': colorstring = colorstring[1:]
    # get bytes in reverse order to deal with PIL quirk
    colorstring = colorstring[-2:] + colorstring[2:4] + colorstring[:2]
    # finally, make it numeric
    color = int(colorstring, 16)
    return color

def PILColorToRGB(pil_color):
    """ convert a PIL-compatible integer into an (r, g, b) tuple """
    hexstr = '%06x' % pil_color
    # reverse byte order
    r, g, b = hexstr[4:], hexstr[2:4], hexstr[:2]
    r, g, b = [int(n, 16) for n in (r, g, b)]
    return (r, g, b)

def PILColorToHTMLColor(pil_integer):
    return RGBToHTMLColor(PILColorToRGB(pil_integer))

def RGBToPILColor(rgb_tuple):
    return HTMLColorToPILColor(RGBToHTMLColor(rgb_tuple))












def add_object_to_node(treeObject, Object, node, Objecttype="sequence"):

	data=treeObject.node(node).get_data()
	data.comment={Objecttype:copy.deepcopy(Object)}
	#data.comment[Objecttype].seq=data.comment[Objecttype].seq.tomutable()
	data.comment[Objecttype]=data.comment[Objecttype].tomutable()
	treeObject.node(node).set_data(data)

	return treeObject

def add_object_to_all_nodes(treeObject, Object, node, Objecttype="sequence"):
	
	for daughter in treeObject.node(node).get_succ():
		treeObject=add_object_to_all_nodes(treeObject, Object, daughter, Objecttype=Objecttype)
	
	treeObject=add_object_to_node(treeObject, Object, node, Objecttype=Objecttype)
	
	return treeObject





def get_downstream_nodes(treeObject, node):

	node_list=[]

	def add_successors(treeObject, node, node_list):
		daughters=treeObject.node(node).get_succ()
		
		for daughter in daughters:
			node_list=add_successors(treeObject, daughter, node_list)
		
		node_list.append(node)
		return node_list

	node_list=add_successors(treeObject, node, node_list)
	
	return node_list


def get_upstream_nodes(treeObject, node):

	node_list=[]

	def add_predecessors(treeObject, node, node_list):
		if node!=treeObject.root:
			parent=treeObject.node(node).get_prev()
			node_list=add_predecessors(treeObject, parent, node_list)
		node_list.append(node)
		return node_list

	node_list=add_predecessors(treeObject, node, node_list)
	
	return node_list



###############################################################################################
# Function to print a tree as a string. Fixes float problem in Bio.Nexus.Trees.Tree.to_string #
###############################################################################################


def tree_to_string(treeObject, support_as_branchlengths=False,branchlengths_only=False,plain=True,plain_newick=False,ladderize=None):

	"""Return a paup compatible tree line.
	
	tree_to_string(treeObject,support_as_branchlengths=False,branchlengths_only=False,plain=True)
	"""
	# if there's a conflict in the arguments, we override plain=True
	if support_as_branchlengths or branchlengths_only:
		plain=False
	treeObject.support_as_branchlengths=support_as_branchlengths
	treeObject.branchlengths_only=branchlengths_only
	treeObject.plain=plain

	def make_info_string(data,terminal=False):
		"""Creates nicely formatted support/branchlengths."""
		# CHECK FORMATTING
		if treeObject.plain: # plain tree only. That's easy.
			return ''
		elif treeObject.support_as_branchlengths: # support as branchlengths (eg. PAUP), ignore actual branchlengths
			if terminal: # terminal branches have 100% support
				return ':%s' % str(treeObject.max_support)
			else:
				return ':%s' % str(data.support)
		elif treeObject.branchlengths_only: # write only branchlengths, ignore support
			return ':%s' % str(data.branchlength)
		else: # write suport and branchlengths (e.g. .con tree of mrbayes)
			if terminal:
				return ':%s' % str(data.branchlength)
			else:
				if data.branchlength is not None and data.support is not None: # we have blen and suppport
					return '%s:%s' % (str(data.support),str(data.branchlength))
				elif data.branchlength is not None: # we have only blen
					return '0:%s' % str(data.branchlength)
				elif data.support is not None: # we have only support
					return ':0'
					#return '%s:0' % str(data.support)
				else:
					return ':0'
					#return '0:0'

	def ladderize_nodes(nodes,ladderize=None):
		"""Sorts node numbers according to the number of terminal nodes."""
		if ladderize in ['left','LEFT','right','RIGHT']:
			succnode_terminals=[(treeObject.count_terminals(node=n),n) for n in nodes]
			succnode_terminals.sort()
			if (ladderize=='right' or ladderize=='RIGHT'):
				succnode_terminals.reverse()
			if succnode_terminals:
				succnodes=zip(*succnode_terminals)[1]
			else:
				succnodes=[]
		else:
			succnodes=nodes
			succnodes.sort()
			return succnodes
   
	def newickize(node,ladderize=None): 
		"""Convert a node tree to a newick tree recursively.""" 
   
		if not treeObject.node(node).succ:	#terminal 
			return treeObject.node(node).data.taxon+make_info_string(treeObject.node(node).data,terminal=True) 
		else: 
			succnodes=ladderize_nodes(treeObject.node(node).succ,ladderize=ladderize) 
			subtrees=[newickize(sn,ladderize=ladderize) for sn in succnodes] 
			return '(%s)%s' % (','.join(subtrees),make_info_string(treeObject.node(node).data)) 
					   
	treeline=['tree'] 
	if treeObject.name: 
		treeline.append(treeObject.name) 
	else: 
		treeline.append('a_tree') 
	treeline.append('=') 
	if treeObject.weight != 1: 
		treeline.append('[&W%s]' % str(round(float(treeObject.weight),3))) 
	if treeObject.rooted: 
		treeline.append('[&R]') 
	succnodes=ladderize_nodes(treeObject.node(treeObject.root).succ) 
	subtrees=[newickize(sn,ladderize=ladderize) for sn in succnodes] 
	treeline.append('(%s)' % ','.join(subtrees)) 
	if plain_newick: 
		return treeline[-1] 
	else: 
		return ' '.join(treeline)+';'
		





###############################################################################################
# Function to print a tree as a string. Fixes float problem in Bio.Nexus.Trees.Tree.to_string #
###############################################################################################


def tree_to_figtree_string(treeObject, support_as_branchlengths=False,branchlengths_only=False,plain=True,plain_newick=False,ladderize=None, node_colours={}):

	"""Return a figtree compatible tree line.
	
	tree_to_string(treeObject,support_as_branchlengths=False,branchlengths_only=False,plain=True)
	"""
	# if there's a conflict in the arguments, we override plain=True
	if support_as_branchlengths or branchlengths_only:
		plain=False
	treeObject.support_as_branchlengths=support_as_branchlengths
	treeObject.branchlengths_only=branchlengths_only
	treeObject.plain=plain

	def make_info_string(node, data,terminal=False, node_colours={}):
		"""Creates nicely formatted support/branchlengths."""
		if node_colours.has_key(node):
			colourbit="[&!color=%s,label=0.0]" % (RGBToHTMLColor(node_colours[node]))
		else:
			colourbit=""
		
		#print node, colourbit, node_colours
		
		# CHECK FORMATTING
		if treeObject.plain: # plain tree only. That's easy.
			return ''
		elif treeObject.support_as_branchlengths: # support as branchlengths (eg. PAUP), ignore actual branchlengths
			if terminal: # terminal branches have 100% support
				return '%s:%s' % (colourbit, str(treeObject.max_support))
			else:
				return '%s:%s' % (colourbit, str(data.support))
		elif treeObject.branchlengths_only: # write only branchlengths, ignore support
			return '%s:%s' % (colourbit, str(data.branchlength))
		else: # write suport and branchlengths (e.g. .con tree of mrbayes)
			if terminal:
				return '%s:%s' % (colourbit, str(data.branchlength))
			else:
				if data.branchlength is not None and data.support is not None: # we have blen and suppport
					return '%s%s:%s' % (str(data.support),colourbit, str(data.branchlength))
				elif data.branchlength is not None: # we have only blen
					return '0%s:%s' % (colourbit, str(data.branchlength))
				elif data.support is not None: # we have only support
					return '%s:0' % (colourbit)
					#return '%s:0' % str(data.support)
				else:
					return '%s:0' % (colourbit)
					#return '0:0'

	def ladderize_nodes(nodes,ladderize=None):
		"""Sorts node numbers according to the number of terminal nodes."""
		if ladderize in ['left','LEFT','right','RIGHT']:
			succnode_terminals=[(treeObject.count_terminals(node=n),n) for n in nodes]
			succnode_terminals.sort()
			if (ladderize=='right' or ladderize=='RIGHT'):
				succnode_terminals.reverse()
			if succnode_terminals:
				succnodes=zip(*succnode_terminals)[1]
			else:
				succnodes=[]
		else:
			succnodes=nodes
			succnodes.sort()
			return succnodes
   
	def newickize(node,ladderize=None, node_colours={}): 
		"""Convert a node tree to a newick tree recursively.""" 
   
		if not treeObject.node(node).succ:	#terminal 
			return treeObject.node(node).data.taxon+make_info_string(node, treeObject.node(node).data,terminal=True, node_colours=node_colours) 
		else: 
			succnodes=ladderize_nodes(treeObject.node(node).succ,ladderize=ladderize) 
			subtrees=[newickize(sn,ladderize=ladderize, node_colours=node_colours) for sn in succnodes]
			return '(%s)%s' % (','.join(subtrees),make_info_string(node, treeObject.node(node).data, node_colours=node_colours))
	
	
	treeline=['#NEXUS\n']
	treeline.append('begin taxa;\n')
	
	treeline.append('\tdimensions ntax='+str(len(treeObject.get_terminals()))+';\n')
	treeline.append('\ttaxlabels\n')
	for taxon in treeObject.get_terminals():
		if node_colours.has_key(taxon):
			treeline.append('\t'+treeObject.node(taxon).get_data().taxon+' [&!color='+RGBToHTMLColor(node_colours[taxon])+']\n')
		else:
			treeline.append('\t'+treeObject.node(taxon).get_data().taxon+'\n')
	treeline.append(';\nend;\n')		
	treeline.append('\nbegin trees;\n\n')		   
	treeline.append('tree')
	if treeObject.name: 
		treeline.append(treeObject.name) 
	else: 
		treeline.append('a_tree') 
	treeline.append('=') 
	if treeObject.weight != 1: 
		treeline.append('[&W%s]' % str(round(float(treeObject.weight),3))) 
	if treeObject.rooted: 
		treeline.append('[&R]') 
	succnodes=ladderize_nodes(treeObject.node(treeObject.root).succ) 
	subtrees=[newickize(sn,ladderize=ladderize, node_colours=node_colours) for sn in succnodes] 
	treeline.append('(%s)' % ','.join(subtrees)) 
	treeline[-1]=treeline[-1]+';'
	treeline.append('\nend;\n')
	treeline.append('\nbegin figtree;\n')
	treeline.append('set appearance.branchLineWidth=5.0;')
	treeline.append('\nend;\n')
	
	if plain_newick: 
		return treeline[-1] 
	else: 
		return ' '.join(treeline)
	
	
	






#######################################################
# Function to check if the tree has any branchlengths #
#######################################################

def has_branchlengths(treeObject,node=None):
	"""Returns True if any of the nodes has data.support != None."""
	for n in treeObject._walk(node):
		if treeObject.node(n).data.branchlength>0:
			return True
	else:
		return False 






#############################################
# Function to root the tree on the midpoint #
#############################################


def midpoint_root(treeObject):

	
	if not has_branchlengths(treeObject):
		print "Tree has no branch lengths, so cannot be midpoint rooted"
		return treeObject



	rootnode=treeObject.root
	
	#function to connect subtrees
	def connect_subtree(parent,child):
		"""Hook subtree starting with node child to parent."""
		for i,branch in enumerate(treeObject.unrooted):
			if parent in branch[:2] and child in branch[:2]:
				branch=treeObject.unrooted.pop(i)
				break 
		else:
			DoError('Unable to connect nodes for rooting: nodes %d and %d are not connected' % (parent,child))
		treeObject.link(parent,child)
		treeObject.node(child).data.branchlength=branch[2]
		treeObject.node(child).data.support=branch[3]
		#now check if there are more branches connected to the child, and if so, connect them
		child_branches=[b for b in treeObject.unrooted if child in b[:2]]
		for b in child_branches:
			if child==b[0]:
				succ=b[1]
			else:
				succ=b[0]
			connect_subtree(child,succ) 
	
	
	#identify all terminal nodes
	
	terminal_nodes=treeObject.get_terminals()
	
	#for each pair of terminal nodes, calculate the distance between them and find the maximum distance
	
	distances=[]
	max_distance=[0, 0, 0]
	
	for x, terminal in enumerate(terminal_nodes):
		for secondterminal in terminal_nodes[x+1:]:
			if float(treeObject.distance(secondterminal, terminal))>max_distance[0]:
				max_distance=[float(treeObject.distance(secondterminal, terminal)), terminal, secondterminal]
			#distances.append([float(treeObject.distance(secondterminal, terminal)), terminal, secondterminal])
	

	
	#create a set for the list of nodes from each of the most distant terminal nodes to the root
	nodes_to_terminal = set(treeObject.trace(rootnode, max_distance[1]))
	nodes_to_secondterminal = set(treeObject.trace(rootnode, max_distance[2]))
	#identify nodes in these two sets that are different (i.e. part of the path between the nodes)
	nodes_joining_terminals = nodes_to_terminal.symmetric_difference(nodes_to_secondterminal)
	#identify nodes common to both paths (one of these will be the join between the two paths)
	intersection_nodes=nodes_to_terminal.intersection(nodes_to_secondterminal)
	
	
	#add the node joining the paths to get the complete path between the two most distant terminal nodes
	for node in intersection_nodes:
		daughters=treeObject.node(node).get_succ()
		x=0
		for daughter in daughters:
			if daughter in nodes_joining_terminals:
				x=x+1
		if x>1:
			nodes_joining_terminals.add(node)
			break
	
	
	
#	if the root is part of the path between the most distant nodes, remove it from the list
	#nodes_joining_terminals.delete(rootnode)
	
	#Calculate which branch is the halfway point between the two most distant nodes
	outgroup_brlen=max_distance[0]
	ingroup_brlen=max_distance[0]
	
	for node in nodes_joining_terminals:
		if (max_distance[0]/2)-float(treeObject.distance(node, max_distance[1]))>0:
			if ((max_distance[0]/2)-float(treeObject.distance(node, max_distance[1])))<outgroup_brlen:
				outgroup_node=node
				outgroup_brlen=((max_distance[0]/2)-float(treeObject.distance(node, max_distance[1])))
				try:
					outgroup_support=float(treeObject.node(node).data.support)
				except:
					 outgroup_support=0.0
					 
		else:
			if (float(treeObject.distance(node, max_distance[1]))-(max_distance[0]/2))<ingroup_brlen:
				ingroup_node=node
				ingroup_brlen=(float(treeObject.distance(node, max_distance[1]))-(max_distance[0]/2))
				try:
					ingroup_support=float(treeObject.node(node).data.support)
				except:
					 ingroup_support=0.0
					 
	#is this correct???
	support=ingroup_support


#	#if the root is one of the nodes identified, exit
#	if rootnode in b[:2]:#need to fix what this does
#		print "Root branch has not changed"
#		return treeObject

	
	#unroot the tree

	treeObject.unroot()

	
	# now we find the branch that connects outgroup and ingroup 

	for i,b in enumerate(treeObject.unrooted): 
		if ingroup_node in b[:2] and outgroup_node in b[:2]: 
			root_branch=treeObject.unrooted.pop(i) 
			break
	else: 
		DoError('Unrooted and rooted Tree do not match') #sometimes crashes here... if midoint is on root branch???


	
	# now we destroy the old tree structure, but keep node data. Nodes will be reconnected according to new outgroup
	for n in treeObject.all_ids():
		treeObject.node(n).prev=None
		treeObject.node(n).succ=[]
	
	# now we just add both subtrees (outgroup and ingroup) branch for branch
	root=Nodes.Node(data=Trees.NodeData())			# new root	
	treeObject.add(root)
	#if rootnode!=ingroup_node and rootnode!=outgroup_node:	# add to tree description
	#	root.id=rootnode								# replace new root id with original one
	treeObject.root=root.id						   # set as root
	
	
	treeObject.unrooted.append([root.id,ingroup_node,ingroup_brlen,support])  # add branch to ingroup to unrooted tree
	treeObject.unrooted.append([root.id,outgroup_node,outgroup_brlen,support])   # add branch to outgroup to unrooted tree
	connect_subtree(root.id,ingroup_node)	  # add ingroup
	connect_subtree(root.id,outgroup_node)	 # add outgroup

	
	# if theres still a lonely node in treeObject.chain, then it's the old root, and we delete it
	oldroot=[i for i in treeObject.all_ids() if treeObject.node(i).prev is None and i!=treeObject.root]
	if len(oldroot)>1:
		raise TreeError('Isolated nodes in tree description: %s' % ','.join(oldroot))
	elif len(oldroot)==1:
		treeObject.kill(oldroot[0])
		

	
	return treeObject.root
	



#######################################################################
# Function to change gaps ("-") to unknowns ("N") in alignment object #
#######################################################################

def gap_to_unknown(alignmentObject, unknown="N"):

	for record in alignmentObject:
		seqstring=str(record.seq.tomutable()).replace("-", unknown)
		record.seq=Seq(seqstring)
	
	return alignmentObject



##################################################
# Function to sum the branch lengths on the tree #
##################################################

def get_total_tree_length(treeObject, node=-1, length=0):
		
	if node==-1:
		node=treeObject.root

	daughters=treeObject.node(node).get_succ()
	
	for daughter in daughters:
		length=get_total_tree_length(treeObject, node=daughter, length=length)
	
	
	length=length+treeObject.node(node).get_data().branchlength

	return length




##############################################################
# Function to parsimoneously reconstruct all sites on a tree #
##############################################################
	
	
def parsimonious_sequence_reconstruction(treeObject, alignmentObject, transformation="acctran", sequence_Objecttype="sequence", locations=[]):


	def fitch(treeObject, sitenumber, node=-1, sequence_Objecttype="sequence"):
		
		if node==-1:
			node=treeObject.root
		
		daughters=treeObject.node(node).get_succ()

		if treeObject.node(node).get_data().comment[sequence_Objecttype][sitenumber]=="-":
			return treeObject
		
		if treeObject.is_internal(node):
		
		
			for daughter in daughters:
				
				treeObject=fitch_algorithm(treeObject, sitenumber, daughter)

				if treeObject.node(daughter).get_data().comment[sequence_Objecttype][sitenumber]=="-":
					continue
					
			first_daughter_states=set(ambiguity_to_bases[treeObject.node(daughters[0]).get_data().comment[sequence_Objecttype][sitenumber]])
			second_daughter_states=set(ambiguity_to_bases[treeObject.node(daughters[1]).get_data().comment[sequence_Objecttype][sitenumber]])
			
			if "-" in first_daughter_states:
				first_daughter_states.remove("-")
			if "-" in second_daughter_states:
				second_daughter_states.remove("-")
			
		   	node_states=first_daughter_states.intersection(second_daughter_states)
		   	if len(node_states)==0:
		   		node_states=first_daughter_states.union(second_daughter_states)

			node_states=list(node_states)
			
			node_states.sort()
			
			data=treeObject.node(node).get_data()
			data.comment[sequence_Objecttype][sitenumber]=bases_to_ambiguity[''.join(node_states)]
			treeObject.node(node).set_data(data)
		
		return treeObject
	
		
	
	def checknodes(treeObject, node, sitenumber):
		daughters=treeObject.node(node).get_succ()
		print node, treeObject.node(node).get_data().comment["sequence"][sitenumber]
		for daughter in daughters:
			checknodes(treeObject, daughter, sitenumber)
	

		
	def reduce_indel_locations(treeObject, node=-1, sequence_Objecttype="sequence", transformation="acctran"):
	
		if node==-1:
			node=treeObject.root
		
		daughters=treeObject.node(node).get_succ()
		
		node_data=treeObject.node(node).get_data()
		
		
		for daughter in daughters:
			
			treeObject=reduce_indel_locations(treeObject, node=daughter, sequence_Objecttype=sequence_Objecttype, transformation=transformation)
			
			daughter_data=treeObject.node(daughter).get_data()
			
			for indel in ["deletion_locations", "insertion_locations"]:
				
				locations=[]
				
				if daughter_data.comment.has_key(indel):
					#print indel, daughter_data.comment[indel]
					location=0
					while location <len( daughter_data.comment[indel]):
						start=daughter_data.comment[indel][location]
						end=start
						location+=1
						while location<len(daughter_data.comment[indel]) and daughter_data.comment[indel][location]==end+1:
							end=end+1
							location+=1
						locations.append([start, end])
					#print daughter, indel, locations
					daughter_data.comment[indel]=locations
					treeObject.node(daughter).set_data(daughter_data)
	
		return treeObject			
	
	
	
	
	def sankoff(treeObject, sitenumber, node=-1, sequence_Objecttype="sequence"):
		
		if node==-1:
			node=treeObject.root
		
		daughters=treeObject.node(node).get_succ()
	
		if treeObject.is_terminal(node):
			
			node_data=treeObject.node(node).get_data()
			node_data.comment["sankoff"]={}
			for state in transition_matrix.keys():
				node_data.comment["sankoff"][state]=100
			for state in ambiguity_to_bases[node_data.comment["sequence"][sitenumber]]:
				node_data.comment["sankoff"][state]=0
			treeObject.node(node).set_data(node_data)
			
		else:
			node_data=treeObject.node(node).get_data()
			node_data.comment["sankoff"]={}
			for state in transition_matrix.keys():
				node_data.comment["sankoff"][state]=0
			for daughter in daughters:
				
				sankoff(treeObject, sitenumber, daughter)
			 	
			 	
			 	daughter_sankoff=treeObject.node(daughter).get_data().comment["sankoff"]
			 	
			 	for state in transition_matrix.keys():
			 		min_cost=100
					for comparison_state in daughter_sankoff.keys():
						cost=transition_matrix[state][comparison_state]+daughter_sankoff[comparison_state]
						if cost<min_cost:
							min_cost=cost
						
					node_data.comment["sankoff"][state]+=min_cost
				
			treeObject.node(node).set_data(node_data)
					
		return treeObject	
		
		
		
	def sankoff_second_traversal(treeObject, sitenumber, node=-1, transformation="acctran", sequence_Objecttype="sequence"):
		
		if node==-1:
			node=treeObject.root
			
		daughters=treeObject.node(node).get_succ()	
		
		if node==treeObject.root:
			min_cost_states=set()
		 	min_cost=100
		 	
		 	
			for state in transition_matrix.keys():
				cost=treeObject.node(node).get_data().comment["sankoff"][state]
		 		if cost<min_cost:
					min_cost=cost
					min_cost_states=set(state)
				elif cost==min_cost:
					min_cost_states.add(state)
			
			root_states=set()
			
			if len(min_cost_states)>1 and "-" in min_cost_states:
				min_cost_states.remove("-")
			
			if len(min_cost_states)>1:# and transformation=="deltran":
				
				parent_states_min_cost={}
			 	
			 	for state in min_cost_states:
			 		parent_states_min_cost[state]=[set()]
			 		daughter_min_cost_states=set()
			 		daughter_min_cost=1000
					for comparison_state in transition_matrix.keys():
						cost=0
						for daughter in daughters:
							cost+=transition_matrix[state][comparison_state]+treeObject.node(daughter).get_data().comment["sankoff"][comparison_state]
						if cost<daughter_min_cost:
							daughter_min_cost=cost
							daughter_min_cost_states=set(comparison_state)
							parent_states_min_cost[state]=[set(comparison_state),cost]
						elif cost==daughter_min_cost:
							daughter_min_cost_states.add(comparison_state)
							parent_states_min_cost[state][0].add(comparison_state)
				
				min_cost=100	
				for state in parent_states_min_cost.keys():
					if parent_states_min_cost[state][1]<min_cost:
						min_cost=parent_states_min_cost[state][1]
						
				root_states=set()
				for state in parent_states_min_cost.keys():
					if parent_states_min_cost[state][1]==min_cost:
						for base in parent_states_min_cost[state][0]:
							root_states.add(base)
				
					
				if len(root_states)>0:
					min_cost_states=root_states
			
			node_data=treeObject.node(node).get_data()
			
			possible_states=list(min_cost_states)
			#shuffle(possible_states)
			if len(possible_states)>1:
				#print sitenumber, possible_states
				min_cost_state=possible_states[0]
				min_cost=100
				for state in possible_states:
					for comparison_state in transition_matrix.keys():
						cost=transition_matrix[state][comparison_state]+treeObject.node(daughters[0]).get_data().comment["sankoff"][comparison_state]
						if cost<min_cost:
							min_cost=cost
							min_cost_state=state
				possible_states=[min_cost_state]
				#print sitenumber, possible_states
				
			node_data.comment["sequence"][sitenumber]=possible_states[0]
				
			treeObject.node(node).set_data(node_data)
			
			
			
		
	
		node_data=treeObject.node(node).get_data()
		node_state=node_data.comment["sequence"][sitenumber]
		node_state_set=set(node_state)
		
		for daughter in daughters:
			
		 	daughter_data=treeObject.node(daughter).get_data()
		 	daughter_sankoff=daughter_data.comment["sankoff"]
		 	
		 	
	 		min_cost_states=set()
	 		min_cost=100
			for comparison_state in daughter_sankoff.keys():
				cost=transition_matrix[node_state][comparison_state]+daughter_sankoff[comparison_state]
				if cost<min_cost:
					min_cost=cost
					min_cost_states=set(comparison_state)
				elif cost==min_cost:
					min_cost_states.add(comparison_state)
			
			if len(min_cost_states)>1 and "-" in min_cost_states:
				min_cost_states.remove("-")
			
			if len(min_cost_states)>1:
			
				if transformation=="acctran":
					min_cost_states=min_cost_states.difference(node_state_set)
				elif transformation=="deltran":
					min_cost_states_new=min_cost_states.intersection(node_state_set)
					if len(min_cost_states_new)==0:
						min_cost_states=min_cost_states
					else:
						min_cost_states=min_cost_states_new
			
			possible_states=list(min_cost_states)
			shuffle(possible_states)
#			if sitenumber>5512 and sitenumber<6811 and daughter==6 or daughter==56:
#				print sitenumber, node, daughter, possible_states, "HERE"
			daughter_data.comment["sequence"][sitenumber]=possible_states[0]
				
			treeObject.node(daughter).set_data(daughter_data)
				
			
			sankoff_second_traversal(treeObject, sitenumber, daughter, transformation=transformation)

	
		return treeObject	
	
	
	
	
	
	
	
	def call_snps(treeObject, node, sitenumber, sequence_Objecttype="sequence", transformation="acctran"):
		
		
		#get all daughters of the current node
		daughters=treeObject.node(node).get_succ()
		#shuffle(daughters)
		#
		node_data=treeObject.node(node).get_data()
		
		
		
		for daughter in daughters:
			
						
			daughter_state=treeObject.node(daughter).get_data().comment[sequence_Objecttype][sitenumber]
			node_state=treeObject.node(node).get_data().comment[sequence_Objecttype][sitenumber]
			
			
			#print node, node_ambiguity, daughter, daughter_ambiguity, new_daughter_states, daughter_states
			
			if daughter_state!=node_state:
#				if daughter==6 or daughter==5:
#					print sitenumber+1, "from node", node, "to node", daughter, node_state, daughter_state
				if node_state=="-" and daughter_state!="-":
					#print "insertion at base", sitenumber+1, "from node", node, "to node", daughter, node_ambiguity, daughter_ambiguity
					daughter_data=treeObject.node(daughter).get_data()
					if not daughter_data.comment.has_key("insertion_locations"):
						daughter_data.comment["insertion_locations"]=[]
					daughter_data.comment["insertion_locations"].append(sitenumber)
					treeObject.node(daughter).set_data(daughter_data)
				elif node_state!="-" and daughter_state=="-":
					#print "deletion at base", sitenumber+1, "from node", node, "to node", daughter, node_ambiguity, daughter_ambiguity
					daughter_data=treeObject.node(daughter).get_data()
					if not daughter_data.comment.has_key("deletion_locations"):
						daughter_data.comment["deletion_locations"]=[]
					daughter_data.comment["deletion_locations"].append(sitenumber)
					treeObject.node(daughter).set_data(daughter_data)
				else:
					#print "change at base", sitenumber+1, "from node", node, "to node", daughter, node_ambiguity+"->"+daughter_ambiguity
					daughter_data=treeObject.node(daughter).get_data()
					if not daughter_data.comment.has_key("SNP_locations"):
						daughter_data.comment["SNP_locations"]={}
					
					SNP=Si_SNPs_temp.SNP()
					SNP.position=sitenumber
					SNP.parent=node
					SNP.daughter=daughter
					SNP.parent_base=ambiguity_to_bases[node_data.comment[sequence_Objecttype][sitenumber]][:]
					SNP.daughter_base=ambiguity_to_bases[daughter_data.comment[sequence_Objecttype][sitenumber]][:]
					
					
					daughter_data.comment["SNP_locations"][sitenumber]=SNP
					
					treeObject.node(daughter).set_data(daughter_data)

				
			

			if treeObject.is_internal(daughter):
				call_snps(treeObject, daughter, sitenumber, sequence_Objecttype=sequence_Objecttype, transformation=transformation)

		#print ambiguous_snps_at_root_node

		return treeObject
	
	
	
	
		
	rootnode=treeObject.root
	
	
	taxa={}
	
	for count, taxon in enumerate(alignmentObject):
		taxa[taxon.name]=count
	
	
	taxadata={}
	substitutions={}
	
	

	
	for node in treeObject.get_terminals():
		data=treeObject.node(node).get_data()
		#data.comment[sequence_Objecttype].seq=alignmentObject[taxa[data.taxon]].seq
		
		mutable_seq=alignmentObject[taxa[data.taxon]].seq.upper().tomutable()
		data.comment[sequence_Objecttype]=mutable_seq
		treeObject.node(node).set_data(data)
	
	#treeObject.display()
	
	if locations==[]:
		locations=range(alignmentObject.get_alignment_length())
	
	#print alignmentObject
	count=0
	total=0.0
	hundredth=float(len(locations))/100
	for columnnumber in locations:
		count=count+1
		if count>=hundredth:
			total=total+count
			count=0
			print "%.0f%% complete\r" % (100*(total/len(locations))),
			sys.stdout.flush()
		
		#when we get the new version of biopython, this can be changed to slice notation alignmentObject[:,base]
		column=alignmentObject.get_column(columnnumber)
		
		column=column.replace("N","").replace("?","").replace("X","")
		firstbase=column[0]
		for base in column[1:]:
			if base!=firstbase:#can add a speedup here by only reconstructing informative sites?
				
				treeObject=sankoff(treeObject, columnnumber, rootnode)
				treeObject=sankoff_second_traversal(treeObject, columnnumber, rootnode, transformation=transformation)
				treeObject=call_snps(treeObject, rootnode, columnnumber, sequence_Objecttype=sequence_Objecttype, transformation=transformation)
				break
	
	treeObject=reduce_indel_locations(treeObject, transformation=transformation, node=rootnode, sequence_Objecttype=sequence_Objecttype)
	
	return treeObject

		







def get_SNPs_from_tree(treeObject, node=-1, sequence_Objecttype="sequence", transformation="acctran", SNP_locations={}):
		
	if node==-1:
		node=treeObject.root
		
	#get all daughters of the current node
	daughters=treeObject.node(node).get_succ()

	#get the data from the node
	node_data=treeObject.node(node).get_data()
	
	
	#for each daughter, recurse to daughter and work out SNPs on that branch
	
	
	for daughter in daughters:
		
		SNP_locations=get_SNPs_from_tree(treeObject, node=daughter, sequence_Objecttype="sequence", transformation="acctran", SNP_locations=SNP_locations)
		
		daughter_data=treeObject.node(daughter).get_data()
		
		if daughter_data.comment.has_key("SNP_locations"):
			for location in daughter_data.comment["SNP_locations"]:
				
				if not SNP_locations.has_key(location):
					SNP_locations[location]={"SNPs":[]}
				
				SNP=Si_SNPs_temp.SNP()
				SNP.parent=node
				SNP.daughter=daughter
				SNP.parent_base=ambiguity_to_bases[node_data.comment[sequence_Objecttype][location]][:]
				SNP.daughter_base=ambiguity_to_bases[daughter_data.comment[sequence_Objecttype][location]][:]
				SNP_locations[location]["SNPs"].append(SNP)
				
				#({"parent":node, "daughter":daughter, "from":ambiguity_to_bases[node_data.comment[sequence_Objecttype][location]][:], "to":ambiguity_to_bases[node_data.comment[sequence_Objecttype][location]][:]})
			
		
		
	return SNP_locations
		
	








	
################################################################
# Function to make branchlength the number of SNPs on a branch #
################################################################
	
	
def branchlengths_to_SNP_count(treeObject, node=-1, lengthtype="SNP_locations"):
	
	if node==-1:
		node=treeObject.root
	
	daughters=treeObject.node(node).get_succ()
		
	for daughter in daughters:
		treeObject=branchlengths_to_SNP_count(treeObject, daughter, lengthtype)	
	
	data=treeObject.node(node).get_data()
	if data.comment.has_key(lengthtype):
		data.branchlength=len(data.comment[lengthtype])
	else:
		data.branchlength=0
	treeObject.node(node).set_data(data)

	return treeObject
	
	
############################################
# Function to make support the node number #
############################################
	
	
def support_to_node_names(treeObject, node=-1):
	
	if node==-1:
		node=treeObject.root
	
	daughters=treeObject.node(node).get_succ()
	
	for daughter in daughters:
		treeObject=support_to_node_names(treeObject, daughter)
		
	
	data=treeObject.node(node).get_data()
	data.support=node
	treeObject.node(node).set_data(data)

	return treeObject
	
	
	

		
		
def run_RAxML(alignmentObject, bootstrap=100, model="GTRGAMMA", cleanup=True):

	#Make a random name for the RAxML run
	chars = string.ascii_letters + string.digits
	tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))

	#remove any previous RAxML runs with this tempname (this should never happen)
	filelist=glob.glob('RAxML_*.'+tmpname)
	for file in filelist:
		os.remove(file)
	
	#convert names so that phylip format doesn't cut them off at 10 characters
	name_dict={}
	count=1
	for record in alignmentObject:
		name="seq"+str(count)
		name_dict[name]=record.id
		record.id=name	
		count+=1
	
	#print alignment to phylip format
	output_handle=open(tmpname+".phy", "w")
	AlignIO.write([alignmentObject], output_handle, "phylip")
	output_handle.close()
	
	#Run RAxML
	sys.stdout.flush()
	if bootstrap==0:
		print "Running RAxML phylogeny with "+model+" model of evolution..."
		os.system(RAXML_DIR+"raxmlHPC -f d -s "+tmpname+".phy -m "+model+" -n "+tmpname+" > /dev/null 2>&1")
		outputname="RAxML_result."+tmpname
	else:
		print "Running RAxML phylogeny with "+model+" model of evolution and "+str(bootstrap)+" bootstrap replicates..."
		os.system("RAxML -f a -x "+str(randrange(1,99999))+" -p "+str(randrange(1,99999))+" -# "+str(bootstrap)+" -m "+model+" -s "+tmpname+".phy -n "+tmpname+" > /dev/null 2>&1")
		outputname="RAxML_bipartitions."+tmpname
	
	#extract stats from raxml output files
	if  bootstrap==0:
		treestats=os.popen('grep "Inference\[0\]" RAxML_info.'+tmpname).read()
		alpha=float(treestats.split(":")[2].split()[0])
	else:
		treestats=os.popen('grep -e "^alpha\|^Tree-Length\|^rate" RAxML_info.'+tmpname).read()
		alpha=float(treestats.split("\n")[0].strip().split(":")[1].strip())
	
	#open the RAxML output tree and convert it to a tree object
	try:
		tree_string = open(outputname).read()
	except IOError:
		DoError("Cannot open tree file")
		
	#convert names back to originals
	for record in alignmentObject:
		tree_string=tree_string.replace(record.id+":", name_dict[record.id]+":")
		record.id=name_dict[record.id]
	
	#Convert tree string to a tree object
	RAxML_tree = Trees.Tree(tree_string, rooted=True)
	
	#Cleanup
	if cleanup:
		filelist=glob.glob('RAxML_*.'+tmpname)
		for file in filelist:
			os.remove(file)
		os.remove(tmpname+".phy")
		if os.path.isfile(tmpname+".phy.reduced"):
			os.remove(tmpname+".phy.reduced")



	return RAxML_tree



def run_phyML(alignmentObject, bootstrap=100, datatype="DNA", model="GTR", gamma=True, pinvar=False, cleanup=True):

	#Make a random name for the RAxML run
	chars = string.ascii_letters + string.digits
	tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))

	#remove any previous RAxML runs with this tempname (this should never happen)
#	filelist=glob.glob('phyML_*.'+tmpname)
#	for file in filelist:
#		os.remove(file)
	
	#convert names so that phylip format doesn't cut them off at 10 characters
	name_dict={}
	count=1
	for record in alignmentObject:
		name="seq"+str(count)
		name_dict[name]=record.id
		record.id=name	
		count+=1
	
	#print alignment to phylip format
	output_handle=open(tmpname+".phy", "w")
	AlignIO.write([alignmentObject], output_handle, "phylip")
	output_handle.close()
	
	#create phyML input string
	
	
	if datatype=="DNA":
		dtype="0"
		k="e"
	else:
		dtype="1"
		k=""
	
	if pinvar:
		i="e"
	else:
		i="0"	
	
	if gamma:
		g="4 e"
	else:
		g="1 1"
	
	#Run phyML
	inputstring=' '.join(["phyml", tmpname+".phy", dtype, "i 1", str(bootstrap), model, k, i, g, "BIONJ y y > /dev/null 2>&1"])
	os.system(inputstring)
	

	#open the phyML output tree and convert it to a tree object
	try:
		tree_string = open(tmpname+".phy_phyml_tree.txt").read()
	except IOError:
		DoError("Cannot open tree file")
		
	#convert names back to originals
	for record in alignmentObject:
		tree_string=tree_string.replace(record.id+":", name_dict[record.id]+":")
		record.id=name_dict[record.id]
	
	#Convert tree string to a tree object
	phyML_tree = Trees.Tree(tree_string, rooted=True)
	
	#Cleanup
	if cleanup:
		filelist=glob.glob(tmpname+".phy_phyml*")
		for file in filelist:
			os.remove(file)
		os.remove(tmpname+".phy")

	return phyML_tree



def get_ref_to_alignment_translations(reference_name, alignmentObject):
	
	refnum=-1
	for x, taxon in enumerate(alignmentObject):
		if taxon.id==reference_name:
			refnum=x
			break
	
	if refnum==-1:
		raise AlignError("")
	
	ref_to_alignment={}
	alignment_to_ref={}
	
	refpos=0
	
	for site in range(alignmentObject.get_alignment_length()):
		if alignmentObject[refnum][site] not in ["-", "N", "?", "X"]:
			ref_to_alignment[refpos]=site
			alignment_to_ref[site]=refpos
			refpos+=1
	
	return ref_to_alignment, alignment_to_ref
	
	
	

	
	
def change_reference_location_to_alignment_location(seqFeature, translation_dict):

	#newFeature=copy.deepcopy(seqFeature)
	
	
	startoffset=translation_dict[int(seqFeature.location.nofuzzy_start)]-int(seqFeature.location.nofuzzy_start)
	endoffset=(translation_dict[int(seqFeature.location.nofuzzy_end)-1]-int(seqFeature.location.nofuzzy_end)+1)
	
	seqFeature.location=FeatureLocation(start = seqFeature.location._start._shift(startoffset),end = seqFeature.location._end._shift(endoffset))
	#seqFeature.location.start=seqFeature.location.start._shift(shift)
	
	#what about fuzzy positions???
	
	for subFeature in seqFeature.sub_features:
		subFeature=change_reference_location_to_alignment_location(subFeature, translation_dict)
	
	
	return seqFeature
	




def annotate_SNPs(treeObject, annotationObject, node=-1):

	if node==-1:
		node=treeObject.root
		
	daughters=treeObject.node(node).get_succ()
	
	node_data=treeObject.node(node).get_data()
	
	for daughter in daughters:
		treeObject=annotate_SNPs(treeObject, annotationObject, daughter)
		
		daughter_data=treeObject.node(daughter).get_data()
		
		if not treeObject.is_internal(daughter):
			daughter_node=daughter_data.taxon
		else:
			daughter_node=daughter
		
		if node==treeObject.root:
			node_node="root"
		else:
			node_node=node
			
		if daughter_data.comment.has_key("SNP_locations"):
		
			SNPlocations=daughter_data.comment["SNP_locations"].keys()
#			
#			SNPlocations.sort()

			loc_count=0
			

			for SNPlocation in SNPlocations:
				#if feature.location.nofuzzy_end>self.position and feature.location.nofuzzy_start<self.position:
					
				
			
				daughter_data.comment["SNP_locations"][SNPlocation].get_annotation_info(annotationObject, treeObject.node(node).get_data().comment["sequence"], treeObject.node(daughter).get_data().comment["sequence"], node_node, daughter_node)
				treeObject.node(daughter).set_data(daughter_data)
			
		
		
		
	
	return treeObject








def dNdS_per_branch(treeObject, node=-1):

	if node==-1:
		node=treeObject.root
	
	
	
	daughters=treeObject.node(node).get_succ()
	node_data=treeObject.node(node).get_data()
	node_revcompseq=Si_SNPs_temp.revcomp(str(treeObject.node(node).get_data().comment["sequence"][:]))
	for daughter in daughters:
		
		treeObject=dNdS_per_branch(treeObject, node=daughter)
		
		daughter_data=treeObject.node(daughter).get_data()

		daughter_revcompseq=Si_SNPs_temp.revcomp(str(treeObject.node(daughter).get_data().comment["sequence"][:]))
		
		for x, feature in enumerate(node_data.comment["annotation"]):
			
			daughter_feature=daughter_data.comment["annotation"][x]
			
			if feature["frameshift"] or daughter_feature["frameshift"]:
				#print "gene is not original length in node or daughter... skipping"
				continue
			
			if feature["strand"]==1:
				
				node_seq=str(treeObject.node(node).get_data().comment["sequence"][feature["location"][0]:feature["location"][1]])
				daughter_seq=str(treeObject.node(daughter).get_data().comment["sequence"][daughter_feature["location"][0]:daughter_feature["location"][1]])
				
				
				try:
					dndsout=Si_SNPs_temp.dnbyds(node_seq,daughter_seq)
				except ValueError:
					daughter_data.comment["annotation"][x]["frameshift"]=True
					#print daughter_data.comment["annotation"][x]
				except StandardError:
					node_data.comment["annotation"][x]["frameshift"]=True
					#print node_data.comment["annotation"][x]
				else:
					daughter_data.comment["annotation"][x]["dNdS"]=dndsout
				
				#daughter_data.comment["annotation"][x]["dNdS"]=Si_SNPs_temp.dnbyds(node_seq,daughter_seq)
				
		
			else:
				
				node_seqlen=len(str(treeObject.node(node).get_data().comment["sequence"]))
				daughter_seqlen=len(str(treeObject.node(daughter).get_data().comment["sequence"]))
			
				#node_start, node_end=Si_SNPs_temp.find_gene_limits(Si_SNPs_temp.revcomp(str(treeObject.node(node).get_data().comment["sequence"])),  seqlen-feature.location.nofuzzy_end,seqlen-feature.location.nofuzzy_start, 1)#need to work out start in revcomp
				
				node_seq=node_revcompseq[node_seqlen-feature["location"][1]:node_seqlen-feature["location"][0]]
				daughter_seq=daughter_revcompseq[daughter_seqlen-daughter_feature["location"][1]:daughter_seqlen-daughter_feature["location"][0]]
				
				

				try:
					dndsout=Si_SNPs_temp.dnbyds(node_seq,daughter_seq)
				except ValueError:
					daughter_data.comment["annotation"][x]["frameshift"]=True
					#print daughter_data.comment["annotation"][x]
				except StandardError:
					node_data.comment["annotation"][x]["frameshift"]=True
					#print node_data.comment["annotation"][x]
				else:
					daughter_data.comment["annotation"][x]["dNdS"]=dndsout

#				 Si_SNPs_temp.dnbyds(Si_SNPs_temp.revcomp(str(treeObject.node(node).get_data().comment["sequence"][feature.location.nofuzzy_start:feature.location.nofuzzy_end])),Si_SNPs_temp.revcomp(str(treeObject.node(daughter).get_data().comment["sequence"][feature.location.nofuzzy_start:feature.location.nofuzzy_end])))
		treeObject.node(daughter).set_data(daughter_data)	
		treeObject.node(node).set_data(node_data)
	
	return treeObject
	
	


	

	
def apply_annotation_to_root(treeObject, annotationObject):

	node=treeObject.root
		
	
	node_data=treeObject.node(node).get_data()

		
	seqlen=len(str(node_data.comment["sequence"]))
	revcompseq=Si_SNPs_temp.revcomp(str(node_data.comment["sequence"][:]))
		
	if not node_data.comment.has_key("annotation"):
		node_data.comment["annotation"]=[]

					
	for feature in annotationObject.features:
		
		if feature.strand==1:
			
			node_start, node_end=Si_SNPs_temp.find_gene_limits(str(node_data.comment["sequence"][feature.location.nofuzzy_start:feature.location.nofuzzy_end+9999]),0,feature.location.nofuzzy_end-feature.location.nofuzzy_start)
			
			node_start=node_start+feature.location.nofuzzy_start
			node_end=node_end+feature.location.nofuzzy_start
			

		else:
			
			node_start, node_end=Si_SNPs_temp.find_gene_limits(revcompseq[seqlen-feature.location.nofuzzy_end:(seqlen-feature.location.nofuzzy_start)+9999],  0,feature.location.nofuzzy_end-feature.location.nofuzzy_start)
			
			node_start=feature.location.nofuzzy_end-node_start
			node_end=feature.location.nofuzzy_end-node_end
			
			
			node_start_new=node_end
			node_end=node_start
			node_start=node_start_new

					
		new_feature={}
		
		new_feature["location"]=[node_start, node_end]
		

		
		for key in ["primary_name", "gene", "systematic_id", "locus_tag"]:
			if feature.qualifiers.has_key(key):
				new_feature["name"]=feature.qualifiers[key][0]
		if not new_feature.has_key("name"):
			new_feature["name"]="no name"
		
		new_feature["strand"]=feature.strand
		
		
		
	
		if node_start==feature.location.nofuzzy_start and node_end==feature.location.nofuzzy_end:
			new_feature["frameshift"]=False
		else:
			new_feature["frameshift"]=True
		
		node_data.comment["annotation"].append(new_feature)
		
	
	treeObject.node(node).set_data(node_data)
	
			
	
	return treeObject






def apply_annotation_to_branches(treeObject, annotationObject, node=-1, handle=None):

	if node==-1:
		node=treeObject.root
		apply_annotation_to_root(treeObject, annotationObject)
	
	def recurse_root_annotation_across_tree(treeObject, node=-1, handle=None):
		
		daughters=treeObject.node(node).get_succ()
		
		node_data=treeObject.node(node).get_data()
		
			
		if node==treeObject.root:
			parent_name="root"
		else:
			parent_name=node
		
		for daughter in daughters:
			
			
			daughter_data=treeObject.node(daughter).get_data()
			seqlen=len(str(daughter_data.comment["sequence"]))
			revcompseq=Si_SNPs_temp.revcomp(str(daughter_data.comment["sequence"][:]))
			
			if not daughter_data.comment.has_key("annotation"):
				daughter_data.comment["annotation"]=[]
			
			if treeObject.is_internal(daughter):
				daughter_name=daughter
			else:
				daughter_name=daughter_data.taxon
					
			for feature in node_data.comment["annotation"]:

				if feature["strand"]==1:
					
					daughter_start, daughter_end=Si_SNPs_temp.find_gene_limits(str(daughter_data.comment["sequence"][feature["location"][0]:feature["location"][1]+9999]),0,feature["location"][1]-feature["location"][0])
					daughter_start=daughter_start+feature["location"][0]
					daughter_end=daughter_end+feature["location"][0]
		
				else:
					
					daughter_start, daughter_end=Si_SNPs_temp.find_gene_limits(revcompseq[seqlen-feature["location"][1]:(seqlen-feature["location"][0])+9999],  0,feature["location"][1]-feature["location"][0])#need to work out start in revcomp
					
					#daughter_start=daughter_start+feature["location"][0]
					#daughter_end=daughter_end+feature["location"][0]
					
					daughter_start=feature["location"][1]-daughter_start
					daughter_end=feature["location"][1]-daughter_end
					
					
					daughter_start_new=daughter_end
					daughter_end=daughter_start
					daughter_start=daughter_start_new
					

					#print daughter_start, daughter_end, feature["location"][0], feature["location"][1]
					
		
				
				if handle!=None and ((feature["strand"]==1 and daughter_end!=feature["location"][1]) or (feature["strand"]==-1 and daughter_start!=feature["location"][0])):
	#				print feature["strand"], daughter_start, daughter_end
	#				print feature["strand"], feature["location"][0], feature["location"][1]		
					if daughter_end>seqlen:
						daughter_end=seqlen-1
					if feature["strand"]==1:
						print >> handle, "FT   CDS             "+str(daughter_start+1)+".."+str(daughter_end)
					else:
						print >> handle, "FT   CDS             complement("+str(daughter_start+1)+".."+str(daughter_end)+")"
						
					print >>handle, 'FT                   /primary_name="'+feature["name"]+'"'
					print >>handle, 'FT                   /node="'+str(parent_name)+'->'+str(daughter_name)+'"'
					if (feature["strand"]==1 and daughter_end>feature["location"][1]) or (feature["strand"]==-1 and daughter_start<feature["location"][0]):
						print >> handle, 'FT                   /colour=4'
					else:
						print >> handle, 'FT                   /colour=11'
				
				new_feature={}
				
				new_feature["location"]=[daughter_start, daughter_end]
				
				new_feature["name"]=feature["name"]
				
				new_feature["strand"]=feature["strand"]
				
				
				
				#print node_start, feature.location.nofuzzy_start, node_end, feature.location.nofuzzy_end
				
				if daughter_start==feature["location"][0] and daughter_end==feature["location"][1]:
					new_feature["frameshift"]=False
				else:
					new_feature["frameshift"]=True
				
				daughter_data.comment["annotation"].append(new_feature)
				
			
			treeObject.node(daughter).set_data(daughter_data)
			
			recurse_root_annotation_across_tree(treeObject, node=daughter, handle=handle)
		
		return treeObject
	

		
	recurse_root_annotation_across_tree(treeObject, node=node, handle=handle)
	
	return treeObject














	
	
def calculate_branch_dNdS(treeObject, handle, node=-1):

	if node==-1:
		node=treeObject.root
	
	daughters=treeObject.node(node).get_succ()
	
	
	
	for daughter in daughters:
		
		treeObject=calculate_branch_dNdS(treeObject, handle, node=daughter)
		node_data=treeObject.node(daughter).get_data()
	
		N=0.0
		S=0.0
		pN=0.0
		pS=0.0
		Nd=0.0
		Sd=0.0
		varianceS=0.0
		varianceN=0.0
		z=0.0
		CDSlen=0
		
		for x, feature in enumerate(node_data.comment["annotation"]):
		
	
		
			if feature.has_key("dNdS"):
				
				N=N+feature["dNdS"]["N"]
				S=S+feature["dNdS"]["S"]
				pN=pN+feature["dNdS"]["pN"]
				pS=pS+feature["dNdS"]["pS"]
				Nd=Nd+feature["dNdS"]["Nd"]
				Sd=Sd+feature["dNdS"]["Sd"]
				CDSlen=CDSlen+feature["dNdS"]["CDSlen"]
		
		if N!=0:
			pN=Nd/N
		if S!=0:
			pS=Sd/S
		
		
		if pS==0:
			#print "No sites are synonymous."
			
			dN=0.0
			dS=0.0
			
			
		elif pS<0.75 and pN<0.75:
			dS=(-3*(math.log(1 - ((pS*4)/3))))/4
			dN=(-3*(math.log(1 - ((pN*4)/3))))/4
			if dN==-0.0:
				dN=0.0
			varianceS=(9 * pS * (1 - pS))/(((3 - (4 *pS)) **2) * CDSlen);
			varianceN=(9 * pN * (1 - pN))/(((3 - (4 *pN)) **2) * CDSlen);
			z=(dN - dS) / math.sqrt(varianceS + varianceN)
			
			node_data.comment["dNdS"]={'N':N, 'S':S, 'dN':dN, 'dS':dS, 'pN':pN, 'pS':pS, 'varS':varianceS, 'varN':varianceN, 'z':z, 'CDSlen':CDSlen}
		
		else:
			#print "Too divergent for JC! Using pN/pS instead."
			dS=pS
			dN=pN
		
		node_data.comment["dNdS"]={'N':N, 'S':S, 'dN':dN, 'dS':dS, 'pN':pN, 'pS':pS, 'varS':varianceS, 'varN':varianceN, 'z':z, 'CDSlen':CDSlen}
		
		
		if not treeObject.is_internal(daughter):
			daughter_node=node_data.taxon
		else:
			daughter_node=str(daughter)
		
		if node==treeObject.root:
			node_node="root"
		else:
			node_node=str(node)
		
		try:
			dnds=dN/dS
		except ZeroDivisionError:
			dnds=0.0
		
		print >> handle, node_node+"->"+daughter_node, dnds
		
		treeObject.node(daughter).set_data(node_data)
		
	return treeObject
			



def write_tab_output(treeObject, handle, node=-1, colour_snps_by="synonymous"):

	if node==-1:
		node=treeObject.root
	
	daughters=treeObject.node(node).get_succ()
	
	node_data=treeObject.node(node).get_data()
	
	for daughter in daughters:
		
		write_tab_output(treeObject, handle, node=daughter, colour_snps_by=colour_snps_by)
		
		daughter_data=treeObject.node(daughter).get_data()
		
		if not treeObject.is_internal(daughter):
			daughter_node=daughter_data.taxon
		else:
			daughter_node=str(daughter)
		
		if node==treeObject.root:
			node_node="root"
		else:
			node_node=str(node)
	
		if daughter_data.comment.has_key("SNP_locations"):
			for SNP_location in daughter_data.comment["SNP_locations"]:
				daughter_data.comment["SNP_locations"][SNP_location].write_tab_format(handle, treeObject.get_taxa(node_id=daughter) ,colourby=colour_snps_by)
		
		
		
		if daughter_data.comment.has_key("insertion_locations"):
			for insertion_location in daughter_data.comment["insertion_locations"]:
			
				print >> handle, 'FT   insertion       '+str(insertion_location[0]+1)+".."+str(insertion_location[1]+1)
				print >> handle, 'FT                   /node="'+str(node_node)+'->'+str(daughter_node)+'"'
				print >> handle, 'FT                   /colour=7'
					
					
		if daughter_data.comment.has_key("deletion_locations"):
			for deletion_location in daughter_data.comment["deletion_locations"]:
				print >> handle, 'FT   deletion        '+str(deletion_location[0]+1)+".."+str(deletion_location[1]+1)
				print >> handle, 'FT                   /node="'+str(node_node)+'->'+str(daughter_node)+'"'
				print >> handle, 'FT                   /colour=10'


		
		if daughter_data.comment.has_key("annotation") and node_data.comment.has_key("annotation"):
		
			for x, feature in enumerate(daughter_data.comment["annotation"]):
		
				if feature["location"][1]!=node_data.comment["annotation"][x]["location"][1]:
		#				print feature["strand"], daughter_start, daughter_end
		#				print feature["strand"], feature["location"][0], feature["location"][1]		
	#				if daughter_data.comment.["annotation"][x]["location"][1]>daughter_data.comment.["annotation"][x]["location"][1]:
	#					daughter_end=seqlen-1
					if feature["strand"]==1:
						print >> handle, "FT   CDS             "+str(feature["location"][0]+1)+".."+str(feature["location"][1]+1)
					else:
						print >> handle, "FT   CDS             complement("+str(feature["location"][0]+1)+".."+str(feature["location"][1]+1)+")"
						
					print >>handle, 'FT                   /primary_name="'+node_data.comment["annotation"][x]["name"]+'"'
					print >>handle, 'FT                   /node="'+str(node_node)+'->'+str(daughter_node)+'"'
					if feature["location"][1]>node_data.comment["annotation"][x]["location"][1]:
						print >> handle, 'FT                   /colour=8'
					else:
						print >> handle, 'FT                   /colour=11'









def identify_homoplasies(treeObject):
	
	
	def add_SNP_locations_to_list(treeObject, node, SNPlist):
		
		daughters=treeObject.node(node).get_succ()
		
		for daughter in daughters:
		
			data=treeObject.node(daughter).get_data().comment
			
			if data.has_key("SNP_locations"):
				for location in data["SNP_locations"]:
					if not SNPlist.has_key(location):
						SNPlist[location]=[]
					SNPlist[location].append([node, daughter, data["SNP_locations"][location].parent_base, data["SNP_locations"][location].daughter_base])
			
			if treeObject.is_internal(daughter):
				SNPlist=add_SNP_locations_to_list(treeObject, daughter, SNPlist)
	
		
	
		return SNPlist
	
	
	SNPlist={}
	SNPlist=add_SNP_locations_to_list(treeObject, treeObject.root, SNPlist)
	#homoplasies={}

	for key in SNPlist.keys():
		if len(SNPlist[key])>1:
			for x, SNP in enumerate(SNPlist[key]):
				for SNPb in SNPlist[key][x+1:]:
				
					if bases_to_ambiguity[''.join(SNP[3])]==bases_to_ambiguity[''.join(SNPb[3])]:
						#print "Base", key, "convergence", str(SNP[0])+"->"+str(SNP[1]), "and", str(SNPb[0])+"->"+str(SNPb[1]), "to", bases_to_ambiguity[''.join(SNP[3])]
						
						node_data=treeObject.node(SNP[1]).get_data()
						node_data.comment["SNP_locations"][key].homoplasy=True
						node_data.comment["SNP_locations"][key].homoplasies.append(["c", SNPb[1]])
						treeObject.node(SNP[1]).set_data(node_data)
						node_data=treeObject.node(SNPb[1]).get_data()
						node_data.comment["SNP_locations"][key].homoplasy=True
						node_data.comment["SNP_locations"][key].homoplasies.append(["c", SNP[1]])
						treeObject.node(SNPb[1]).set_data(node_data)
						
						#homoplasies[key]=[[SNP[1]]+get_downstream_nodes(treeObject, SNP[1]),[SNPb[1]]+get_downstream_nodes(treeObject, SNPb[1])]
						
						
						
					elif bases_to_ambiguity[''.join(SNP[3])]==bases_to_ambiguity[''.join(SNPb[2])] and  treeObject.common_ancestor(SNP[1], SNPb[1])==SNPb[1]:
						#print "Base", key, "reversal", str(SNP[0])+"->"+str(SNP[1]), "to", bases_to_ambiguity[''.join(SNP[3])], "and", str(SNPb[0])+"->"+str(SNPb[1]), "to", bases_to_ambiguity[''.join(SNPb[3])]


						
						node_data=treeObject.node(SNP[1]).get_data()
						node_data.comment["SNP_locations"][key].homoplasy=True
						node_data.comment["SNP_locations"][key].homoplasies.append(["r", SNPb[1]])
						treeObject.node(SNP[1]).set_data(node_data)
						node_data=treeObject.node(SNPb[1]).get_data()
						node_data.comment["SNP_locations"][key].homoplasy=True
						node_data.comment["SNP_locations"][key].homoplasies.append(["r", SNP[1]])
						treeObject.node(SNPb[1]).set_data(node_data)
						
					elif bases_to_ambiguity[''.join(SNP[2])]==bases_to_ambiguity[''.join(SNPb[3])] and  treeObject.common_ancestor(SNP[1], SNPb[1])==SNP[1]:
						#print "Base", key, "reversal", str(SNP[0])+"->"+str(SNP[1]), "to", bases_to_ambiguity[''.join(SNP[3])], "and", str(SNPb[0])+"->"+str(SNPb[1]), "to", bases_to_ambiguity[''.join(SNPb[3])]

						node_data=treeObject.node(SNP[1]).get_data()
						node_data.comment["SNP_locations"][key].homoplasy=True
						node_data.comment["SNP_locations"][key].homoplasies.append(["r", SNPb[1]])
						treeObject.node(SNP[1]).set_data(node_data)
						node_data=treeObject.node(SNPb[1]).get_data()
						node_data.comment["SNP_locations"][key].homoplasy=True
						node_data.comment["SNP_locations"][key].homoplasies.append(["r", SNP[1]])
						treeObject.node(SNPb[1]).set_data(node_data)
						
			#print homoplasies			
				
	return treeObject
			



def get_homoplasy_blocks(treeObject, handle):

	root=treeObject.root
	#handle=open("test.tab","w")


	def count_diffs(node1, node2, start, end):

		seq1=treeObject.node(node1).get_data().comment["sequence"][start:end]

		seq2=treeObject.node(node2).get_data().comment["sequence"][start:end]
		numdiffs=0
		for x, ambiguity1 in enumerate(seq1):

			bases1=set(ambiguity_to_bases[ambiguity1])
			bases2=set(ambiguity_to_bases[seq2[x]])

			if "-" in bases1:
				bases1.remove("-")
			elif "-" in bases2:
				bases2.remove("-")
			
			#print bases1, bases2, bases1.union(bases2)
			
			if len(bases1.intersection(bases2))==0:
				numdiffs+=1
		
		#print node1, node2, start, end, numdiffs
		return numdiffs
	
	
	
	def count_diffs_for_seqs(seq1, seq2):
		numdiffs=0
		for x, ambiguity1 in enumerate(seq1):

			bases1=set(ambiguity_to_bases[ambiguity1])
			bases2=set(ambiguity_to_bases[seq2[x]])

			if "-" in bases1:
				bases1.remove("-")
			elif "-" in bases2:
				bases2.remove("-")
			
			#print bases1, bases2, bases1.union(bases2)
			
			if len(bases1.intersection(bases2))==0:
				numdiffs+=1
		
		#print numdiffs
		return numdiffs
	
	
		
		
	def get_consensus(seq1, seq2):
		
		consensus_seq=""
		
		numdiffs=0
		for x, ambiguity1 in enumerate(seq1):
		
			bases1=set(ambiguity_to_bases[ambiguity1])
			bases2=set(ambiguity_to_bases[seq2[x]])
			
			union=bases1.union(bases2)
			
			union_list=list(union)
			union_list.sort()
			
			consensus_seq=consensus_seq+bases_to_ambiguity[''.join(union_list)]
		
		return consensus_seq
	
	
	
	
	def find_blocks(node, final_blocks, block_improvement, block_key):
	
		data=treeObject.node(node).get_data().comment
		#root_branches=treeObject.node(root).get_succ()
		#if node==root_branches[1]:
		#	return
		
		if data.has_key("SNP_locations"):

			locations=data["SNP_locations"].keys()
			locations.sort()

			
			blocks={}
			new_blocks={}
			
			for location in locations:
				if data["SNP_locations"][location].homoplasy:

				   	
					for homoplasy in data["SNP_locations"][location].homoplasies:
						
						if homoplasy[0]=="r" and treeObject.common_ancestor(node, homoplasy[1])==homoplasy[1]:
							
							blocks[location]=[set(),"r", node,treeObject.node(node).get_prev()]
							blocks[location][0]=blocks[location][0].union(set([homoplasy[1]]+get_upstream_nodes(treeObject, homoplasy[1])))
							
							#remove the root from the list, as it has no prev
							
							if treeObject.root in blocks[location][0]:
								blocks[location][0].remove(treeObject.root)
						elif homoplasy[0]=="c":
							#print "cf"
							blocks[location]=[set(),"f",node, treeObject.node(node).get_prev()]
							#blocks[location][0]=blocks[location][0].union(set([homoplasy[1]]+get_downstream_nodes(treeObject, homoplasy[1])))
							blocks[location][0]=blocks[location][0].union(set(get_downstream_nodes(treeObject, homoplasy[1])))
							#blocks[location][0]=blocks[location][0].union(set([treeObject.node(homoplasy[1]).get_prev()]+get_downstream_nodes(treeObject, homoplasy[1])))#+get_downstream_nodes(treeObject, treeObject.node(homoplasy[1]).get_prev())))
							#print get_downstream_nodes(treeObject, homoplasy[1])
						else:
							#print "rf"
							blocks[location]=[set(),"f",node, treeObject.node(node).get_prev()]
							#blocks[location][0]=blocks[location][0].union(set([homoplasy[1]]+get_downstream_nodes(treeObject, homoplasy[1])))
							blocks[location][0]=blocks[location][0].union(set(get_downstream_nodes(treeObject, homoplasy[1])))
							#blocks[location][0]=blocks[location][0].union(set([treeObject.node(homoplasy[1]).get_prev()]+get_downstream_nodes(treeObject, treeObject.node(homoplasy[1]).get_prev())))
							
			#print blocks
			keys=blocks.keys()
			keys.sort()
			x=0
			while x+1<len(keys):


				if blocks[keys[x]][1]==blocks[keys[x+1]][1]:
					
					start=keys[x]
					end=keys[x+1]+1
				
				
					seq1=treeObject.node(blocks[keys[x]][2]).get_data().comment["sequence"][start:end]

					
					#print start, end, blocks[keys[x]][3], blocks[keys[x]][2], 
					
					seq2=treeObject.node(blocks[keys[x]][3]).get_data().comment["sequence"][start:end]
					
					origin_len=count_diffs_for_seqs(seq1, seq2)
					
					
					#origin_len=count_diffs(blocks[keys[x]][2], blocks[keys[x]][3], keys[x], keys[x+1]+1)
					
					possible_alternatives=blocks[keys[x]][0].intersection(blocks[keys[x+1]][0])
	
					alt_lengths=set()
	
					min_alt_len=10000000
					
					for alt in possible_alternatives:
						#print start, end, str(treeObject.node(alt).get_prev())+"->"+str(alt), blocks[keys[x]][2], 
						
						if treeObject.node(alt).get_prev()==treeObject.root:
							
							seq2=get_consensus(treeObject.node(treeObject.node(treeObject.root).get_succ()[0]).get_data().comment["sequence"][start:end], treeObject.node(treeObject.node(treeObject.root).get_succ()[1]).get_data().comment["sequence"][start:end])
						else:
							seq2=get_consensus(treeObject.node(treeObject.node(alt).get_prev()).get_data().comment["sequence"][start:end], treeObject.node(alt).get_data().comment["sequence"][start:end])
						length=count_diffs_for_seqs(seq1,seq2)
						#print length, alt, origin_len
						if length<min_alt_len and length<=origin_len:
							alt_lengths={}
							alt_lengths=set([alt])
							min_alt_len=length
						elif length==min_alt_len:
							alt_lengths.add(alt)
					
					
					#print "locations=", keys[x], keys[x+1]+1, "node=", blocks[keys[x]][2], "len=",  origin_len, alt_lengths#, blocks[keys[x]], blocks[keys[x+1]]
					
					if len(alt_lengths)>0:
						new_blocks[x+1]=[[],blocks[keys[x]][1],blocks[keys[x]][2], blocks[keys[x]][3]]
						new_blocks[x+1][0]=[keys[x], keys[x+1]+1,origin_len, alt_lengths]
				x=x+1

			#print blocks, len(blocks)
			if len(blocks)>0:
				
				last_homoplasy=keys[-1]
				first_homoplasy=keys[0]
				start=0
				end=len(treeObject.node(node).get_data().comment["sequence"])
				#end=50000
			
			
				seq1=treeObject.node(blocks[keys[x]][2]).get_data().comment["sequence"][start:first_homoplasy]+treeObject.node(blocks[keys[x]][2]).get_data().comment["sequence"][last_homoplasy:end]
				seq2=treeObject.node(blocks[keys[x]][3]).get_data().comment["sequence"][start:first_homoplasy]+treeObject.node(blocks[keys[x]][3]).get_data().comment["sequence"][last_homoplasy:end]
				
				
				origin_len=count_diffs_for_seqs(seq1, seq2)
				
				possible_alternatives=blocks[keys[0]][0]

				alt_lengths=set()

				min_alt_len=10000000
				
				for alt in possible_alternatives:
					if treeObject.node(alt).get_prev()==treeObject.root:
						seq2=get_consensus(treeObject.node(treeObject.node(treeObject.root).get_succ()[0]).get_data().comment["sequence"][start:first_homoplasy], treeObject.node(treeObject.node(treeObject.root).get_succ()[1]).get_data().comment["sequence"][start:first_homoplasy])+get_consensus(treeObject.node(treeObject.node(treeObject.root).get_succ()[0]).get_data().comment["sequence"][last_homoplasy:end], treeObject.node(treeObject.node(treeObject.root).get_succ()[1]).get_data().comment["sequence"][last_homoplasy:end])
					else:
						seq2=get_consensus(treeObject.node(treeObject.node(alt).get_prev()).get_data().comment["sequence"][start:first_homoplasy], treeObject.node(alt).get_data().comment["sequence"][start:first_homoplasy])+get_consensus(treeObject.node(treeObject.node(alt).get_prev()).get_data().comment["sequence"][last_homoplasy:end], treeObject.node(alt).get_data().comment["sequence"][last_homoplasy:end])
					length=count_diffs_for_seqs(seq1,seq2)
					if length<min_alt_len and length<=origin_len:
						alt_lengths={}
						alt_lengths=set([alt])
						min_alt_len=length
					elif length==min_alt_len:
						alt_lengths.add(alt)
				

				
				#print "locations=", 0, keys[0]+1, "node=", blocks[keys[0]][2], "len=",  origin_len, alt_lengths, min_alt_len

				if len(alt_lengths)>0:
						new_blocks[0]=[[],blocks[keys[0]][1],blocks[keys[0]][2], blocks[keys[0]][3]]
						new_blocks[0][0]=[0, keys[0]+1,origin_len, alt_lengths]
						new_blocks[x+1]=[[],blocks[keys[-1]][1], blocks[keys[-1]][2], blocks[keys[-1]][3]]
						new_blocks[x+1][0]=[keys[-1], end ,origin_len, alt_lengths]
			
			
			
			
			
			keys=new_blocks.keys()
			keys.sort()
			
			for key in keys:
				if new_blocks.has_key(key+1) and new_blocks[key][1]==new_blocks[key+1][1] and new_blocks[key][0][3].intersection(new_blocks[key+1][0][3]):
					new_blocks[key+1][0][3]=new_blocks[key][0][3].intersection(new_blocks[key+1][0][3])
					new_blocks[key+1][0][0]=new_blocks[key][0][0]
					del new_blocks[key]
			
			#print "NB", new_blocks
			
			for key in new_blocks.keys():
				if len(new_blocks[key][0][3])>1 and new_blocks[key][1]=="f":
#					common_ancestor=treeObject.common_ancestor(list(new_blocks[key][0][3])[0], list(new_blocks[key][0][3])[1])
#					for x in range(len(new_blocks[key][0][3])):
#						if treeObject.common_ancestor(common_ancestor, list(new_blocks[key][0][3])[x]) not in list(new_blocks[key][0][3]):
#							common_ancestor="None"
#							break
					common_ancestors=set(new_blocks[key][0][3])
					#print common_ancestors
					for x in range(len(new_blocks[key][0][3])):
						common_ancestors=common_ancestors.intersection(set([list(new_blocks[key][0][3])[x]]+get_upstream_nodes(treeObject, list(new_blocks[key][0][3])[x])))
						#print "HERE f", common_ancestors
				elif len(new_blocks[key][0][3])>1 and new_blocks[key][1]=="r":
					common_ancestors=set(new_blocks[key][0][3])
					#print common_ancestors
					for x in range(len(new_blocks[key][0][3])):
						common_ancestors=common_ancestors.intersection(set([list(new_blocks[key][0][3])[x]]+get_downstream_nodes(treeObject, list(new_blocks[key][0][3])[x])))
						#print "HERE r", common_ancestors
				elif len(new_blocks[key][0][3])==1:
						common_ancestors=new_blocks[key][0][3]
				
				if len(common_ancestors)>0:
					new_blocks[key][0][3]=common_ancestors
				else:
					print "no common ancestor", new_blocks[key][3]
			
			
			
			
			for key in new_blocks.keys():
			
				start=new_blocks[key][0][0]
				end=new_blocks[key][0][1]
		
				seq1=treeObject.node(node).get_data().comment["sequence"][start:end]

				
				#print start, end, blocks[keys[x]][3], blocks[keys[x]][2], 
				
				seq2=treeObject.node(treeObject.node(node).get_prev()).get_data().comment["sequence"][start:end]
				
				min_count=count_diffs_for_seqs(seq1, seq2)
				#min_count=count_diffs(node, treeObject.node(node).get_prev(), new_blocks[key][0][0], new_blocks[key][0][1])
				node_snps=min_count
				#print node_snps
				min_alt=""
				for alt in new_blocks[key][0][3]:
					if treeObject.node(alt).get_prev()==treeObject.root:
						seq2=get_consensus(treeObject.node(treeObject.node(treeObject.root).get_succ()[0]).get_data().comment["sequence"][start:end], treeObject.node(treeObject.node(treeObject.root).get_succ()[1]).get_data().comment["sequence"][start:end])
					else:
						seq2=get_consensus(treeObject.node(treeObject.node(alt).get_prev()).get_data().comment["sequence"][start:end], treeObject.node(alt).get_data().comment["sequence"][start:end])
					alt_count=count_diffs_for_seqs(seq1, seq2)
					#alt_count=count_diffs(new_blocks[key][2], alt, new_blocks[key][0][0], new_blocks[key][0][1])
					if alt_count<min_count:
						min_count=alt_count
						min_alt=alt
				#print node_snps, min_count, min_alt
				if min_alt!="" and (min_count+10)<node_snps:# and min_count<node_snps:#
					
					
				
					final_blocks[block_key]={"start":new_blocks[key][0][0], "end":new_blocks[key][0][1], "from_node":new_blocks[key][3], "to_node": new_blocks[key][2], "to_rec":min_alt, "from_rec":treeObject.node(min_alt).get_prev(), "rec_snps":min_count, "node":node, "node_snps":node_snps}
					
					block_improvement.append([node_snps-min_count,block_key])
					block_key+=1
						
			
		#print node, new_blocks
			
			
		return final_blocks, block_improvement, block_key	
			
			
			

	def iterate_nodes(node, final_blocks, block_improvement, block_key):
		daughters=treeObject.node(node).get_succ()

		for daughter in daughters:
			final_blocks, block_improvement, block_key=iterate_nodes(daughter, final_blocks, block_improvement, block_key)

		final_blocks, block_improvement, block_key=find_blocks(node, final_blocks, block_improvement, block_key)
		
		return final_blocks, block_improvement, block_key

	final_blocks={}
	block_improvement=[]
	block_key=0

	final_blocks, block_improvement, block_key=iterate_nodes(root, final_blocks, block_improvement, block_key)
	
	
	
	block_improvement.sort()
	block_improvement=block_improvement[::-1]
	boiled_down={}
	
	for x, num in enumerate(block_improvement):
			block_key=num[1]
			
			if not final_blocks.has_key(block_key):
				continue
			
			start=final_blocks[block_key]["start"]
			end=final_blocks[block_key]["end"]
			for numb in block_improvement[x+1:]:
				block_keyb=numb[1]
				
				if not final_blocks.has_key(block_keyb):
					continue
				if not final_blocks.has_key(block_key):
					break	
				
				startb=final_blocks[block_keyb]["start"]
				endb=final_blocks[block_keyb]["end"]
				if (startb>=start and startb<end) or (endb>start and endb<=end) or (startb<=start and endb>=end):
					if (final_blocks[block_key]["to_rec"] in [final_blocks[block_keyb]["to_rec"], final_blocks[block_keyb]["to_node"]]) or (final_blocks[block_keyb]["to_rec"] in [final_blocks[block_key]["to_rec"], final_blocks[block_key]["to_node"]]):
						#print final_blocks[block_key]
						#print final_blocks[block_keyb]
						
						if start<startb:
							new_start=start
						else:
							new_start=startb
						if end>endb:
							new_end=end
						else:
							new_end=endb
						
						#print new_start, new_end
						
						seq1=treeObject.node(final_blocks[block_key]["to_node"]).get_data().comment["sequence"][new_start:new_end]
						seq2=treeObject.node(final_blocks[block_key]["from_node"]).get_data().comment["sequence"][new_start:new_end]
						originalcount=count_diffs_for_seqs(seq1, seq2)
						
						
						seq2=get_consensus(treeObject.node(final_blocks[block_key]["from_rec"]).get_data().comment["sequence"][new_start:new_end], treeObject.node(final_blocks[block_key]["to_rec"]).get_data().comment["sequence"][new_start:new_end])
						reccount=count_diffs_for_seqs(seq1, seq2)
						
						seq1=treeObject.node(final_blocks[block_keyb]["to_node"]).get_data().comment["sequence"][new_start:new_end]
						seq2=treeObject.node(final_blocks[block_keyb]["from_node"]).get_data().comment["sequence"][new_start:new_end]
						originalcountb=count_diffs_for_seqs(seq1, seq2)
						
						seq2=get_consensus(treeObject.node(final_blocks[block_keyb]["from_rec"]).get_data().comment["sequence"][new_start:new_end], treeObject.node(final_blocks[block_keyb]["to_rec"]).get_data().comment["sequence"][new_start:new_end])
						reccountb=count_diffs_for_seqs(seq1, seq2)
						
						
						
						if (originalcount-reccount)>(originalcountb-reccountb) and (originalcount-reccount)>(final_blocks[block_key]["node_snps"]-final_blocks[block_key]["rec_snps"]):
							#print "A", originalcount, reccount, originalcount-reccount, originalcountb, reccountb, originalcountb-reccountb
							final_blocks[block_key]["start"]=new_start
							final_blocks[block_key]["end"]=new_end
							
						elif (originalcount-reccount)>(originalcountb-reccountb) and (originalcountb-reccountb)>(final_blocks[block_key]["node_snps"]-final_blocks[block_key]["rec_snps"]):
							print "B", originalcount, reccount, originalcount-reccount, originalcountb, reccountb, originalcountb-reccountb
							final_blocks[block_keyb]["start"]=new_start
							final_blocks[block_keyb]["end"]=new_end
							final_blocks[block_key]=final_blocks[block_keyb]
						
						del final_blocks[block_keyb]
					
					
						
						
						

	for key in final_blocks.keys():
		if treeObject.is_terminal(final_blocks[key]["to_rec"]):
			to_rec=treeObject.node(final_blocks[key]["to_rec"]).get_data().taxon
		else:
			to_rec=str(final_blocks[key]["to_rec"])
		if ["from_rec"]==treeObject.root:
			from_rec="Root"
		else:
			from_rec=str(final_blocks[key]["from_rec"])
		if treeObject.is_terminal(final_blocks[key]["to_node"]):
			to_node=treeObject.node(final_blocks[key]["to_node"]).get_data().taxon
		else:
			to_node=str(final_blocks[key]["to_node"])
		if ["from_node"]==treeObject.root:
			from_node="Root"
		else:
			from_node=str(final_blocks[key]["from_node"])
		
		
		print >> handle, "FT   recombination   "+str(final_blocks[key]["start"]+1)+".."+str(final_blocks[key]["end"])
		print >> handle, 'FT                   /recipient_branch='+from_node+"->"+to_node
		print >> handle, 'FT                   /donor_branch='+from_rec+"->"+to_rec
		print >> handle, "FT                   /original_branch_SNPs="+str(final_blocks[key]["node_snps"])
		print >> handle, "FT                   /recombination_branch_SNPs="+str(final_blocks[key]["rec_snps"])
#		print "FT   recombination   "+str(final_blocks[key]["start"]+1)+".."+str(final_blocks[key]["end"])
#		print 'FT                   /recipient_branch='+from_node+"->"+to_node
#		print 'FT                   /donor_branch='+from_rec+"->"+to_rec
#		print "FT                   /original_branch_SNPs="+str(final_blocks[key]["node_snps"])
#		print "FT                   /recombination_branch_SNPs="+str(final_blocks[key]["rec_snps"])
	
	handle.close()





def genome_diagram_for_tree(treeObject,filename,ladderize=None, printtype="SNPs", colourby="SNPtype", referenceObject=None):

	def ladderize_nodes(nodes,ladderize=None):
		"""Sorts node numbers according to the number of terminal nodes."""
		if ladderize in ['left','LEFT','right','RIGHT']:
			succnode_terminals=[(treeObject.count_terminals(node=n),n) for n in nodes]
			succnode_terminals.sort()
			if (ladderize=='right' or ladderize=='RIGHT'):
				succnode_terminals.reverse()
			if succnode_terminals:
				succnodes=zip(*succnode_terminals)[1]
			else:
				succnodes=[]
		else:
			succnodes=nodes
			succnodes.sort()
			return succnodes
   
   
	def SNPs_to_features(node, features=[], depth=0, colourby="SNPtype", taxonsets=[]):
		
		node_data=treeObject.node(node).get_data()
		
		
		
		if node_data.comment.has_key("SNP_locations"):
		
			for SNP in node_data.comment["SNP_locations"]:
				
				if colourby=="SNPtype":
					feature=SeqFeature(FeatureLocation(node_data.comment["SNP_locations"][SNP].position,node_data.comment["SNP_locations"][SNP].position), strand=1, type="SNP", qualifiers={"colour":node_data.comment["SNP_locations"][SNP].colour})
					features.append([feature,depth])
				
				elif colourby=="homoplasy":
					if node_data.comment["SNP_locations"][SNP].homoplasy:
						feature=SeqFeature(FeatureLocation(node_data.comment["SNP_locations"][SNP].position,node_data.comment["SNP_locations"][SNP].position), strand=1, type="SNP", qualifiers={"colour":"2"})
						features.append([feature,depth])
#					else:
#						feature=SeqFeature(FeatureLocation(node_data.comment["SNP_locations"][SNP].position,node_data.comment["SNP_locations"][SNP].position), strand=1, type="SNP", qualifiers={"colour":"1"})

			
				elif colourby=="taxa":
					if node_data.comment["SNP_locations"][SNP].homoplasy:
						taxaset=set()
						for taxon in treeObject.get_terminals():
							if treeObject.node(taxon).get_data().comment["sequence"][SNP]==treeObject.node(node).get_data().comment["sequence"][SNP]:
								taxaset.add(taxon)
								
						snpcolour=-1
						x=0
						for x, tset in enumerate(taxonsets):
							if len(tset.intersection(taxaset))==len(taxaset) and len(tset.intersection(taxaset))==len(tset):
								snpcolour=x+1
								break
					
						if snpcolour==-1:
							snpcolour=x+2
							taxonsets.append(taxaset)

						while snpcolour>17:
							snpcolour=snpcolour-17

					#else:
						#snpcolour=0
					
					
						feature=SeqFeature(FeatureLocation(node_data.comment["SNP_locations"][SNP].position,node_data.comment["SNP_locations"][SNP].position), strand=1, type="SNP", qualifiers={"colour":str(snpcolour)})
						features.append([feature,depth])
				
				
		
		return features, taxonsets
		
		
   
	def add_SNP_tracks(node,track_number=1,ladderize=None, previous_features=[], depth=0, colourby="SNPtype", taxonsets=[]): 
		"""Convert a node tree to a newick tree recursively."""
   			
   		
		if treeObject.node(node).succ:
   			succnodes=ladderize_nodes(treeObject.node(node).succ,ladderize=ladderize)
   			succnodes.reverse()
			features, taxonsets=SNPs_to_features(node, previous_features[:], depth, colourby=colourby, taxonsets=taxonsets)
	   		for succnode in succnodes:
	   			track_number, taxonsets=add_SNP_tracks(succnode,track_number,ladderize=ladderize,previous_features=features, depth=depth+1, colourby=colourby, taxonsets=taxonsets)
   		else:
			gd_track_for_features = gd_diagram.new_track(track_number, scale=0, height=1)
			gd_feature_set = gd_track_for_features.new_set()
	
	
			features, taxonsets=SNPs_to_features(node, previous_features[:], depth, colourby=colourby, taxonsets=taxonsets)
			for x, feature in enumerate(features):
				translator = ColorTranslator()
				if colourby=="depth":
					color=translator.artemis_color(str(feature[1]))
				elif feature[0].qualifiers.has_key("colour"):
					color=translator.artemis_color(feature[0].qualifiers["colour"])
				elif feature[0].qualifiers.has_key("color"):
					color=translator.artemis_color(feature[0].qualifiers["color"])
				else:
					color = translator.artemis_color("0")
				
				gd_feature_set.add_feature(feature[0], color=color, label=0, arrowhead_length=0.25)

			track_number+=1
   
   		return track_number, taxonsets
   		
   	
   	
   	
   	
   	def changed_genes_to_features(node, features=[], depth=0):
		
		node_data=treeObject.node(node).get_data()
		
		
		
		if node_data.comment.has_key("annotation"):
			parent=treeObject.node(node).prev
			parent_data=treeObject.node(parent).get_data()
			for x, gene in enumerate(node_data.comment["annotation"]):
				if gene["frameshift"]:
					if (gene["strand"]==1 and  gene["location"][1]<parent_data.comment["annotation"][x]["location"][1]) or (gene["strand"]==-1 and  gene["location"][0]>parent_data.comment["annotation"][x]["location"][0]):
						feature=SeqFeature(FeatureLocation(gene["location"][0],gene["location"][1]), strand=1, type="SNP", qualifiers={"colour":"11", "name":gene["name"]})
						features.append([feature,depth])
					elif (gene["strand"]==1 and  gene["location"][1]>parent_data.comment["annotation"][x]["location"][1]) or (gene["strand"]==-1 and  gene["location"][0]<parent_data.comment["annotation"][x]["location"][0]):
						feature=SeqFeature(FeatureLocation(gene["location"][0],gene["location"][1]), strand=1, type="SNP", qualifiers={"colour":"4", "name":gene["name"]})
						features.append([feature,depth])
				
				
		
		return features
   	
   	
   	
   	
   	
   	def add_changed_gene_tracks(node,track_number=1,ladderize=None, previous_features=[], depth=0): 
		"""Convert a node tree to a newick tree recursively."""
   			
   		
		if treeObject.node(node).succ:
   			succnodes=ladderize_nodes(treeObject.node(node).succ,ladderize=ladderize)
   			succnodes.reverse()
			features=changed_genes_to_features(node, previous_features[:], depth)
	   		for succnode in succnodes:
	   			track_number=add_changed_gene_tracks(succnode,track_number,ladderize=ladderize,previous_features=features, depth=depth+1)
   		else:
			gd_track_for_features = gd_diagram.new_track(track_number, scale=0, height=1)
			gd_feature_set = gd_track_for_features.new_set()
	
			features=changed_genes_to_features(node, previous_features[:], depth)
			for x, feature in enumerate(features):
				translator = ColorTranslator()
				if feature[0].qualifiers.has_key("colour"):
					color=translator.artemis_color(feature[0].qualifiers["colour"][0])
				elif feature[0].qualifiers.has_key("color"):
					color=translator.artemis_color(feature[0].qualifiers["color"][0])
				else:
					color = translator.artemis_color("0")
				#color=translator.artemis_color(str(feature[1]))
				gd_feature_set.add_feature(feature[0], color=color, label_angle=0, label_size=8, label=0, arrowhead_length=0.25, name=feature[0].qualifiers["name"])

			track_number+=1
   
   		return track_number
   	
   	
   	
   	

					   
	gd_diagram = GenomeDiagram.Diagram()
	track_number=1
	
	if referenceObject:
		gd_track_for_features = gd_diagram.new_track(track_number, scale=0, height=2)
		gd_feature_set = gd_track_for_features.new_set()
		translator = ColorTranslator()
		for x, feature in enumerate(referenceObject.features):
	
			if feature.type not in  ["CDS"]:#,"tRNA","repeat_unit"] :
				#Exclude this feature
				continue
			if feature.qualifiers.has_key("colour"):
				color=translator.artemis_color(feature.qualifiers["colour"][0])
			elif feature.qualifiers.has_key("color"):
				color=translator.artemis_color(feature.qualifiers["color"][0])
			else:
				color = translator.artemis_color("5")
			gd_feature_set.add_feature(feature, color=color, label=0)
		
		
		track_number+=1
	
	succnodes=ladderize_nodes(treeObject.node(treeObject.root).succ)
	succnodes.reverse()
	for node in succnodes:
		if printtype in ["SNPs", "taxa"]:
			track_number, taxonsets=add_SNP_tracks(node, track_number=track_number,ladderize=ladderize, colourby=colourby)
		elif printtype=="changed_genes":
			track_number=add_changed_gene_tracks(node, track_number=track_number,ladderize=ladderize)
	for tset in taxonsets:		
		print tset
	
	
	gd_diagram.draw(format="linear", orientation="landscape", pagesize='A3', fragments=1)#, start=0, end=len(record))
	
	gd_diagram.write(filename, "PDF")




def draw_ascii_tree(treeObject, show_nodes=False):

	try:
		terminal_rows, terminal_columns = os.popen('stty size', 'r').read().split()
	except ValueError:
		return

	def get_terminal_nodes(tree):
		node=tree.root
		node_list=[]
		
		def get_daughter_nodes(tree, node, node_list):
		
			daughters=tree.node(node).get_succ()
			for daughter in daughters:
				node_list=get_daughter_nodes(tree, daughter, node_list)
			
			if tree.is_terminal(node):
				node_list.append(node)
			
			return node_list
		
		node_list=get_daughter_nodes(tree, node, node_list)
		
		return node_list

	
	terminal_node_list=get_terminal_nodes(treeObject)
	
	terminal_taxa=[]
	for node in terminal_node_list:
	
		terminal_taxa.append(treeObject.node(node).get_data().taxon)
		
	
	max_label_width = max(len(str(taxon)) for taxon in terminal_taxa)
	drawing_width = int(terminal_columns) - max_label_width - 1
	drawing_height = 2 * len(terminal_node_list) - 1 
	


	def get_nodes(tree):
		node=tree.root
		node_list=[]
		
		def get_daughter_nodes(tree, node, node_list):
		
			daughters=tree.node(node).get_succ()
			for daughter in daughters:
				node_list=get_daughter_nodes(tree, daughter, node_list)
			
			node_list.append(node)
			
			return node_list
		
		node_list=get_daughter_nodes(tree, node, node_list)
		
		return node_list


	def get_col_positions(tree):
		"""Create a mapping of each clade to its column position."""
		#depths = tree.depths()
		
		depths={}
		for node in get_nodes(tree):
			depths[node]=tree.sum_branchlength(node=node)
		
		# If there are no branch lengths, assume unit branch lengths 
		if not max(depths.itervalues()):
			depths = tree.depths(unit_branch_lengths=True)
		# Potential drawing overflow due to rounding -- 1 char per tree layer
		fudge_margin = int(math.ceil(math.log(len(terminal_node_list), 2)))
		cols_per_branch_unit = ((drawing_width - fudge_margin)/ float(max(depths.itervalues())))
		return dict((clade, int(round(blen*cols_per_branch_unit + 0.5))) for clade, blen in depths.iteritems())
   
	def get_row_positions(tree): 
		positions = dict((taxon, 2*idx) for idx, taxon in enumerate(terminal_node_list))
		
		def calc_row(node):
			daughters=tree.node(node).get_succ()
			for daughter in daughters:
				if daughter not in positions: 
					calc_row(daughter) 
			positions[node] = (positions[daughters[0]] + positions[daughters[1]]) / 2 
		calc_row(tree.root) 
		return positions 

	column_positions=get_col_positions(treeObject)
	
	row_positions=get_row_positions(treeObject)
	
	char_matrix = [[' ' for x in range(drawing_width+max_label_width+1)] for y in range(drawing_height)] 
	
	
	def draw_clade(tree, node, startcol): 
		thiscol = column_positions[node] 
		thisrow = row_positions[node] 
		# Draw a horizontal line 
		
		
		for col in range(startcol, thiscol): 
			char_matrix[thisrow][col] = '-' 
			
		if tree.is_internal(node):
			daughters=tree.node(node).get_succ()
			# Draw a vertical line 
			toprow = row_positions[daughters[0]] 
			botrow = row_positions[daughters[-1]] 
			for row in range(toprow+1, botrow): 
				char_matrix[row][thiscol] = '|'
				if row==thisrow and show_nodes:
					for x in range(len(str(node))):
						char_matrix[row][thiscol+x+2]=str(node)[x]
			
			#Add a little bit to the top and bottom?
			
#			if char_matrix[toprow][thiscol]==" ":
#				char_matrix[toprow][thiscol]=","
#			if char_matrix[botrow][thiscol]==" ":
#				char_matrix[botrow][thiscol]="'"

			
			
			
			for daughter in daughters: 
				draw_clade(tree, daughter, thiscol+1) 
		else:
			taxon=tree.node(node).get_data().taxon
			for x in range(len(taxon)):
				char_matrix[thisrow][thiscol+x+1]=taxon[x]
	
	
	
	
	draw_clade(treeObject, treeObject.root, 0) 

	    # Print the complete drawing 
	for idx, row in enumerate(char_matrix): 
		line = ''.join(row)#.rstrip() 
		# Add labels for terminal taxa in the right margin 
#		if idx % 2 == 0: 
#			line += ' ' + str(terminal_taxa[idx/2]) 
		print line 
		

	print
	




def moving_window_recombination_detection(treeObject):
	
	seqlen=len(treeObject.node(treeObject.root).get_data().comment["sequence"])
	
	windowsize=int(float(seqlen)/1000)
	
	if windowsize<10:
		windowsize=10
	elif windowsize>1000:
		windowsize=1000
	
	
	max_distance=0
	for terminal in treeObject.get_terminals():
		distance=treeObject.distance(treeObject.root, terminal)
		if distance>max_distance:
			max_distance=distance
	
	
	def ladderize_nodes(nodes,ladderize=None):
		"""Sorts node numbers according to the number of terminal nodes."""
		if ladderize in ['left','LEFT','right','RIGHT']:
			succnode_terminals=[(treeObject.count_terminals(node=n),n) for n in nodes]
			succnode_terminals.sort()
			if (ladderize=='right' or ladderize=='RIGHT'):
				succnode_terminals.reverse()
			if succnode_terminals:
				succnodes=zip(*succnode_terminals)[1]
			else:
				succnodes=[]
		else:
			succnodes=nodes
			succnodes.sort()
			return succnodes
	
	def get_node_colours(treeObject, ladderize=None):
		
		node=treeObject.root
		
		def get_node_colour(node, leftdist, rightdist, max_distance, num_nodes, node_num, nodecolours={}):
			node_num+=1
			#print node, (leftdist/max_distance)*255, (float(node_num)/num_nodes)*255, (rightdist/max_distance)*255, node_num, num_nodes,leftdist, rightdist, max_distance
			nodecolours[node]=((leftdist/max_distance)*255, (float(node_num)/num_nodes)*255, (rightdist/max_distance)*255)
			#nodecolours[node]=((leftdist/max_distance)*255, 0, (rightdist/max_distance)*255)
			if treeObject.node(node).succ:
	   			daughters=ladderize_nodes(treeObject.node(node).succ,ladderize=ladderize)
	   			daughters.reverse()
	   			
	   			for daughter in daughters:
	   			
					if daughter==daughters[0]:
						newdist=leftdist+treeObject.node(daughter).get_data().branchlength
						nodecolours, node_num=get_node_colour(daughter, newdist, rightdist, max_distance, num_nodes, node_num, nodecolours)
					else:
						newdist=rightdist+treeObject.node(daughter).get_data().branchlength
						nodecolours, node_num=get_node_colour(daughter, leftdist, newdist, max_distance, num_nodes, node_num, nodecolours)
						
					
					
			
			#colour stuff here :-S
			
			
			
			
			return nodecolours, node_num
		
		
		
		num_nodes=treeObject.count_terminals(node)*2
		
		#print max_distance, num_nodes
		
		nodecolours={}
		
		nodecolours, node_num=get_node_colour(node, 0, 0, max_distance, num_nodes, 0, nodecolours={})
		
		#print nodecolours
		
		return nodecolours
		
	
	
	
	
	def get_max_depth(treeObject, ladderize=None):
		node=treeObject.root
	
		def recurse(node, depth, right, left, max_depth, max_right, max_left, ladderize=None):
		
			if depth>max_depth:
				max_depth=depth
			if left>max_left:
				max_left=left
			if right>max_right:
				max_right=right
			
			if treeObject.node(node).succ:
	   			daughters=ladderize_nodes(treeObject.node(node).succ,ladderize=ladderize)
	   			daughters.reverse()
	   			for daughter in daughters:
	   				
		   			if daughter==daughters[0]:
		   				newdepth=depth+1
		   				newleft=left+1
						max_depth, max_right, max_left=recurse(daughter, newdepth, right, newleft, max_depth, max_right, max_left, ladderize=ladderize)
					else:
						newdepth=depth+1
						newright=right+1
						max_depth, max_right, max_left=recurse(daughter, newdepth, newright, left, max_depth, max_right, max_left, ladderize=ladderize)
					
			
			return max_depth, max_right, max_left
		
		
		max_depth, max_right, max_left=recurse(node, 0, 0, 0, 0, 0, 0, ladderize=ladderize)
		
		return max_depth, max_right, max_left
	
	
	def get_node_colours(treeObject, ladderize=None):
		
		node=treeObject.root
		
		max_depth, max_right, max_left=get_max_depth(treeObject)
		#print max_depth, max_right, max_left
		def get_node_colour(node, num_nodes, node_num, nodecolours={}):
			
			#print node, (leftdist/max_distance)*255, (float(node_num)/num_nodes)*255, (rightdist/max_distance)*255, node_num, num_nodes,leftdist, rightdist, max_distance
			#nodecolours[node]=((leftdist/max_distance)*255, (float(node_num)/num_nodes)*255, (rightdist/max_distance)*255)
			#nodecolours[node]=((leftdist/max_distance)*255, 255-((float(rightdist)/max_depth)*255), (rightdist/max_distance)*255)
			#nodecolours[node]=((float(leftdist)/max_left)*255, 255-((float(node_num)/num_nodes)*255), (float(rightdist)/max_right)*255)
			#nodecolours[node]=((float(node_num)/num_nodes)*255, 0, 255-((float(node_num)/num_nodes)*255))
			
			
			
			
			
			
			if treeObject.node(node).succ:
	   			daughters=ladderize_nodes(treeObject.node(node).succ,ladderize=ladderize)
	   			daughters.reverse()
	   			
	   			
	   			nodecolours, node_num=get_node_colour(daughters[0], num_nodes, node_num, nodecolours)
	   			nodeproportion=(float(node_num)/num_nodes)*1275
			
				red=510-nodeproportion
				blue=nodeproportion-510
				green=nodeproportion
	   			if red<-510:
					reddiff=510+red
					red=reddiff*-1
				elif red<0:
					red=0
				elif red>255:
					red=255
				if blue<0:
					blue=0
				elif blue>255:
					blue=255
				if green>255 and green<765:
					green=255
				elif green>765:
					greendiff=765-green
					green=255+greendiff
					if green<0:
						green=0
				
				
				#print red, green, blue
	   			
	   			nodecolours[node]=(red, green, blue)
	   			node_num+=1
				nodecolours, node_num=get_node_colour(daughters[1], num_nodes, node_num, nodecolours)
				
			else:
				nodeproportion=(float(node_num)/num_nodes)*1275
			
				red=510-nodeproportion
				blue=nodeproportion-510
				green=nodeproportion
				if red<-510:
					reddiff=510+red
					red=reddiff*-1
				elif red<0:
					red=0
				elif red>255:
					red=255
				if blue<0:
					blue=0
				elif blue>255:
					blue=255
				if green>255 and green<765:
					green=255
				elif green>765:
					greendiff=765-green
					green=255+greendiff
					if green<0:
						green=0
				
				
				#print red, green, blue
	   			
	   			nodecolours[node]=(red, green, blue)
				node_num+=1
			
			return nodecolours, node_num
		
		
		
		num_nodes=(treeObject.count_terminals(node)*2)-1
		
		#print max_distance, num_nodes
		
		nodecolours={}
		
		nodecolours, node_num=get_node_colour(node, num_nodes, 0, nodecolours={})
		
		#print nodecolours
		
		return nodecolours
	
	
	
	
		
		
		
	
	def print_node_diagram(treeObject, node, nodecolours, besthits):
	
		
	
		def add_coloured_track(node,track_number=1,ladderize=None):
	   		
			if treeObject.node(node).succ:
	   			succnodes=ladderize_nodes(treeObject.node(node).succ,ladderize=ladderize)
	   			succnodes.reverse()
		   		#for succnode in succnodes:
		   		track_number=add_coloured_track(succnodes[0],track_number,ladderize=ladderize)
	   			gd_track_for_features = gd_diagram.new_track(track_number, scale=0, height=1)
				gd_feature_set = gd_track_for_features.new_set()
		
				feature=SeqFeature(FeatureLocation(1,100), strand=1, type="colour_test")
				color = translator.int255_color(nodecolours[node])
				#print nodecolours[node]
				gd_feature_set.add_feature(feature, color=color, label=0)
	
				track_number+=1
		   		track_number=add_coloured_track(succnodes[1],track_number,ladderize=ladderize)
	   		else:
				gd_track_for_features = gd_diagram.new_track(track_number, scale=0, height=1)
				gd_feature_set = gd_track_for_features.new_set()
		
				feature=SeqFeature(FeatureLocation(1,100), strand=1, type="colour_test")
				color = translator.int255_color(nodecolours[node])
				#print nodecolours[node]
				#print color
				gd_feature_set.add_feature(feature, color=color, label=0)
	
				track_number+=1
	   
	   		return track_number
	   		
	   	translator = ColorTranslator()
	   	gd_diagram = GenomeDiagram.Diagram()
		track_number=1
	   	
		add_coloured_track(treeObject.root,track_number,ladderize=None)
	  
		gd_diagram.draw(format="linear", orientation="landscape", pagesize='A3', fragments=1)#, start=0, end=len(record))
		
		filename="test.pdf"
		
		gd_diagram.write(filename, "PDF")
	
	
	def count_diffs(node1, node2, start, end):

		seq1=treeObject.node(node1).get_data().comment["sequence"][start:end]

		seq2=treeObject.node(node2).get_data().comment["sequence"][start:end]
		numdiffs=0
		for x, ambiguity1 in enumerate(seq1):

			bases1=set(ambiguity_to_bases[ambiguity1])
			bases2=set(ambiguity_to_bases[seq2[x]])

			if "-" in bases1:
				bases1.remove("-")
			if "-" in bases2:
				bases2.remove("-")
			
			
			
			if len(bases1)>0 and len(bases2)>0 and len(bases1.intersection(bases2))==0:
				numdiffs+=1
				#print bases1, bases2, bases1.union(bases2)
		
		#print node1, node2, start, end, numdiffs
		return numdiffs
	

	
	
	
	def moving_window_recombination(nodea,nodeb, besthits={}):
		linked_nodes=treeObject.node(nodea).get_succ()+[treeObject.node(nodea).get_prev()]
		
		if not besthits.has_key(nodea):
				besthits[nodea]={}
		
		for x in range(0, seqlen, windowsize):
			
			hit=count_diffs(nodea, nodeb, x, x+windowsize)
			
			#print nodea, nodeb, x, hit
			
			if not besthits[nodea].has_key(x):
				besthits[nodea][x]=[nodeb, hit, x, x+windowsize]
			elif hit<besthits[nodea][x][1]:
				#print x, besthits[x], nodeb, hit
				besthits[nodea][x]=[nodeb, hit, x, x+windowsize]
			elif hit==besthits[nodea][x][1] and besthits[nodea][x][0] not in linked_nodes and nodeb in linked_nodes:
			
				#print x, besthits[x], nodeb, hit, linked_nodes
			
				besthits[nodea][x]=[nodeb, hit, x, x+windowsize]
				
#			else:
#				print x, besthits[x], nodeb, hit
			 
	
		return besthits
	
	
	def get_all_node_distances(node, nodeb, node_distances=[]):
	
		daughters=treeObject.node(nodeb).get_succ()
		
		
		for daughter in daughters:
			
			
			node_distances=get_all_node_distances(node, daughter, node_distances)
		
			node_distances.append([treeObject.distance(node, daughter),daughter])
		
		return node_distances
	
	
	
	def node_recombination(node, besthits={}):
		
		daughters=treeObject.node(node).get_succ()
		
			
		node_distances=[]
		
		node_distances=get_all_node_distances(node, treeObject.root, node_distances)
		node_distances.sort()
		#print node_distances
		
#		pairwise_comparison(daughter)
#		
#		besthits=node_recombination(node, besthits=besthits)
#		print node, besthits
		
		
#		daughters=treeObject.node(nodeb).get_succ()
#		
		for node_data in node_distances:
			daughter=node_data[1]
			
			if node==daughter:
				continue
			
			besthits=moving_window_recombination(node,daughter, besthits)
			#node_recombination(node, nodeb=daughter)
			
				
		
	
			#print node, daughter, besthits
		
		return besthits
	
	
	
	
	def pairwise_comparison(node, total_nodes, node_count=0, besthits={}):
		
		daughters=treeObject.node(node).get_succ()
		node_count=node_count+1
		
		print "%.0f%% complete\r" % (100*(float(node_count)/total_nodes)),
		sys.stdout.flush()
		
		for daughter in daughters:
		
			besthits, node_count=pairwise_comparison(daughter, total_nodes, node_count=node_count,besthits=besthits)
#		

		besthits=node_recombination(node, besthits=besthits)
		#print node, besthits
		
		#linked_nodes=daughters+[treeObject.node(node).get_prev()]
		
#		for besthit in besthits.keys():
#			if not besthits[besthit][0] in linked_nodes:
#				print node, besthit, besthits[besthit]
		
	
		return besthits, node_count
	
	
	
	
	def print_recombination_tree_diagram(treeObject, node, nodecolours, besthits):
	
		depth=0
	
		def add_coloured_track(node,track_number=1,ladderize=None, depth=0):
	   		
			if treeObject.node(node).succ:
	   			succnodes=ladderize_nodes(treeObject.node(node).succ,ladderize=ladderize)
	   			succnodes.reverse()
		   		#for succnode in succnodes:
		   		
		   		track_number=add_coloured_track(succnodes[0],track_number,ladderize=ladderize, depth=depth+treeObject.node(node).get_data().branchlength)
	   			gd_track_for_features = gd_diagram.new_track(track_number, scale=0, height=1)
				gd_feature_set = gd_track_for_features.new_set()
				
				features=[]
				for location in besthits[node].keys():
					if besthits[node][location][0]!=node:
						#print node, location, besthits[node][location]
						color = translator.int255_color(nodecolours[besthits[node][location][0]])
					else:
						#color = translator.int255_color(nodecolours[node])
						color = translator.int255_color((230,230,230))
					
					depthratio=((treeObject.node(node).get_data().branchlength/max_distance)*1000000)+1
					startposition=(float(besthits[node][location][2])/seqlen)
					endposition=(float(besthits[node][location][3])/seqlen)
					windowlen=(float(windowsize)/seqlen)
					totaldepthratio=(depth/max_distance)*1000000
					#print totaldepthratio
					feature=SeqFeature(FeatureLocation((startposition*depthratio)+totaldepthratio,(endposition*depthratio)+totaldepthratio), strand=1, type="colour_test")
					gd_feature_set.add_feature(feature, color=color, label=0)
				
				#print nodecolours[node]
	
				track_number+=1
		   		
		   		track_number=add_coloured_track(succnodes[1],track_number,ladderize=ladderize, depth=depth+treeObject.node(node).get_data().branchlength)
	   		else:
				gd_track_for_features = gd_diagram.new_track(track_number, scale=0, height=1)
				gd_feature_set = gd_track_for_features.new_set()
		
				features=[]
				for location in besthits[node].keys():
					if besthits[node][location][0]!=node:
						#print node, location, besthits[node][location]
						color = translator.int255_color(nodecolours[besthits[node][location][0]])
					else:
						#color = translator.int255_color(nodecolours[node])
						color = translator.int255_color((230,230,230))
					
					depthratio=((treeObject.node(node).get_data().branchlength/max_distance)*1000000)+1
					startposition=(float(besthits[node][location][2])/seqlen)
					endposition=(float(besthits[node][location][3])/seqlen)
					windowlen=(float(windowsize)/seqlen)
					totaldepthratio=(depth/max_distance)*1000000
					#print totaldepthratio
					feature=SeqFeature(FeatureLocation((startposition*depthratio)+totaldepthratio,(endposition*depthratio)+totaldepthratio), strand=1, type="colour_test")
					gd_feature_set.add_feature(feature, color=color, label=0)
				
				#print nodecolours[node]
	
				track_number+=1
	   
	   		return track_number
	   		
	   	translator = ColorTranslator()
	   	gd_diagram = GenomeDiagram.Diagram()
		track_number=1
	   	
		add_coloured_track(treeObject.root,track_number,ladderize=None)
	  
		gd_diagram.draw(format="linear", orientation="landscape", pagesize='A3', fragments=1)#, start=0, end=len(record))
		
		filename="recombination_tree.pdf"
		
		gd_diagram.write(filename, "PDF")
	
	
	
	
	def print_recombination_node_diagram(treeObject, node, nodecolours, besthits):
	
		depth=0
	
		def add_coloured_track(node,track_number=1,ladderize=None, depth=0):
	   		
			if treeObject.node(node).succ:
	   			succnodes=ladderize_nodes(treeObject.node(node).succ,ladderize=ladderize)
	   			succnodes.reverse()
		   		#for succnode in succnodes:
		   		
		   		track_number=add_coloured_track(succnodes[0],track_number,ladderize=ladderize, depth=depth+treeObject.node(node).get_data().branchlength)
	   			gd_track_for_features = gd_diagram.new_track(track_number, scale=0, height=1)
				gd_feature_set = gd_track_for_features.new_set()
				
				features=[]
				for location in besthits[node].keys():
					if besthits[node][location][0]!=node:
						#print node, location, besthits[node][location]
						color = translator.int255_color(nodecolours[besthits[node][location][0]])
					else:
						#color = translator.int255_color(nodecolours[node])
						color = translator.int255_color((230,230,230))
					
					startposition=besthits[node][location][2]
					endposition=besthits[node][location][3]
					#print totaldepthratio
					feature=SeqFeature(FeatureLocation(startposition,endposition), strand=1, type="colour_test")
					gd_feature_set.add_feature(feature, color=color, label=0)
				
				#print nodecolours[node]
	
				track_number+=1
		   		
		   		track_number=add_coloured_track(succnodes[1],track_number,ladderize=ladderize, depth=depth+treeObject.node(node).get_data().branchlength)
	   		else:
				gd_track_for_features = gd_diagram.new_track(track_number, scale=0, height=1)
				gd_feature_set = gd_track_for_features.new_set()
		
				features=[]
				for location in besthits[node].keys():
					if besthits[node][location][0]!=node:
						#print node, location, besthits[node][location]
						color = translator.int255_color(nodecolours[besthits[node][location][0]])
					else:
						#color = translator.int255_color(nodecolours[node])
						color = translator.int255_color((230,230,230))
					
					startposition=besthits[node][location][2]
					endposition=besthits[node][location][3]
					#print totaldepthratio
					feature=SeqFeature(FeatureLocation(startposition,endposition), strand=1, type="colour_test")
					gd_feature_set.add_feature(feature, color=color, label=0)
				
				#print nodecolours[node]
	
				track_number+=1
	   
	   		return track_number
	   		
	   	translator = ColorTranslator()
	   	gd_diagram = GenomeDiagram.Diagram()
		track_number=1
	   	
		add_coloured_track(treeObject.root,track_number,ladderize=None)
	  
		gd_diagram.draw(format="linear", orientation="landscape", pagesize='A3', fragments=1)#, start=0, end=len(record))
		
		filename="recombination_linear.pdf"
		
		gd_diagram.write(filename, "PDF")
	
	
	
	
	
	
	node=treeObject.root
	
	
	
	#besthits={}

	
	nodecolours=get_node_colours(treeObject)
	
	#print nodecolours
	
	
	output=open("recombination_key.tre","w")
	
	print >> output, tree_to_figtree_string(treeObject, False, False, False, False, node_colours=nodecolours)
	
	output.close()
	
	#print_node_diagram(treeObject, treeObject, nodecolours,besthits)
	besthits, node_count=pairwise_comparison(node, len(treeObject.get_terminals())*2)
	
	for key in besthits.keys():
		locations=besthits[key].keys()
		locations.sort()
		
		linked_nodes=treeObject.node(key).get_succ()+[treeObject.node(key).get_prev()]
		
		lastlocationkey=-1
		lastbesthit=-1
		
		for locationkey in locations:
		
			if besthits[key][locationkey][0] in linked_nodes:
				besthits[key][locationkey][0]=key
		
			if besthits[key][locationkey][0]==lastbesthit:
				besthits[key][locationkey][2]=besthits[key][lastlocationkey][2]
				del besthits[key][lastlocationkey]
			lastlocationkey=locationkey
			lastbesthit=besthits[key][locationkey][0]
			if besthits[key][locationkey][3]>seqlen:
				besthits[key][locationkey][3]=seqlen
		
	
	print_recombination_node_diagram(treeObject, node, nodecolours,besthits)
	print_recombination_tree_diagram(treeObject, node, nodecolours,besthits)



