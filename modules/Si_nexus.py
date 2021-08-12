from Bio.Nexus import Trees, Nodes
from Bio import AlignIO
from Bio.Align.Generic import Alignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.Alphabet import IUPAC, Gapped
import sys, string, os, glob, copy
from Si_general import *
from random import *
import Si_SNPs_temp
import math
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram._Colors import ColorTranslator
import time
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/lib/python2.7/site-packages/fisher-0.1.4-py2.7-linux-x86_64.egg']))
#from scipy import stats
import fisher
RAXML = "raxmlHPC"


ambiguity_to_bases={"A":["A"], "C":["C"], "T":["T"], "G":["G"], "M":["A", "C"], "K":["G", "T"], "R":["A", "G"], "Y":["C", "T"], "S":["C", "G"], "W":["A", "T"], "B":["C", "G", "T"], "V":["A", "C", "G"], "H":["A", "C", "T"], "D":["A", "G", "T"], "?":["-", "A", "C", "G", "T"], "N":["-", "A", "C", "G", "T"], "X":["-", "A", "C", "G", "T"], "-":["-"], "a":["-", "A"], "c":["-","C"], "t":["-","T"], "g":["-","G"], "m":["-","A","C"], "k":["-", "G", "T"], "r":["-", "A","G"], "y":["-","C","T"], "s":["-","C","G"], "w":["-","A","T"], "b":["-","C","G","T",], "v":["-","A","C","G"], "h":["-","A","C","T"], "d":["-","A","G","T",], "n":["-","A","C","G","T"]}

bases_to_ambiguity={"A":"A", "C":"C", "T":"T", "G":"G", "AC":"M", "GT":"K", "AG":"R", "CT":"Y", "CG":"S", "AT":"W", "CGT":"B", "ACG":"V", "ACT":"H", "AGT":"D", "ACGT":"N", "-":"-", "?":"N", "X":"N", "-A":"a", "-C":"c", "-T":"t", "-G":"g", "-AC":"m", "-GT":"k", "-AG":"r", "-CT":"y", "-CG":"s", "-AT":"w", "-CGT":"b", "-ACG":"v", "-ACT":"h", "-AGT":"d", "-ACGT":"N"}


#transition_matrix={"A":{"A":0, "C":1, "G":1, "T":1, "-":1, "N":0}, "C":{"A":1, "C":0, "G":1, "T":1, "-":1, "N":0}, "G":{"A":1, "C":1, "G":0, "T":1, "-":1, "N":0}, "T":{"A":1, "C":1, "G":1, "T":0, "-":1, "N":0}, "-":{"A":1, "C":1, "G":1, "T":1, "-":0, "N":0}, "N":{"A":0, "C":0, "G":0, "T":0, "-":0, "N":0}}
transition_matrix={"A":{"A":0, "C":1, "G":1, "T":1, "-":1}, "C":{"A":1, "C":0, "G":1, "T":1, "-":1}, "G":{"A":1, "C":1, "G":0, "T":1, "-":1}, "T":{"A":1, "C":1, "G":1, "T":0, "-":1}, "-":{"A":1, "C":1, "G":1, "T":1, "-":0}}

four_bases=["A", "C", "G", "T"]
missing_data_and_gaps=["N", "?", "-"]
missing_data=["N", "?"]
missing_set=set(["N", "?"])
missing_and_gaps_set=set(["N", "?", "-"])




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
		succnodes=nodes[:]
		succnodes.sort()
		return succnodes

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



###########################################
# Function to create paml baseml.ctl file #
###########################################

def create_baseml_control_file(datafile, treefile, alpha):

	output=open("baseml.ctl","w")
	
	if alpha>10:
		alpha=0

	print >> output, "       seqfile = "+datafile
	print >> output, "       treefile = "+treefile

	print >> output, """

      outfile = mlb       * main result file
        noisy = 9   * 0,1,2,3: how much rubbish on the screen
      verbose = 0   * 1: detailed output, 0: concise output
      runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
                    * 3: StepwiseAddition; (4,5):PerturbationNNI 

        model = 7   * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
                    * 5:T92, 6:TN93, 7:REV, 8:UNREST, 9:REVu; 10:UNRESTu
        
        Mgene = 0   * 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff

*        ndata = 5
        clock = 0   * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
    fix_kappa = 0   * 0: estimate kappa; 1: fix kappa at value below
        kappa = 5  * initial or fixed kappa

    fix_alpha = 1   * 1: estimate alpha; 1: fix alpha at value below"""
    #fix_alpha = 1   * 0: estimate alpha; 1: fix alpha at value below"""
	print >> output, "        alpha =", alpha, "* initial or fixed alpha, 0:infinity (constant rate)"
	print >> output, """        Malpha = 0   * 1: different alpha's for genes, 0: one alpha
        ncatG = 4   * # of categories in the dG, AdG, or nparK models of rates
        nparK = 0   * rate-class models. 1:rK, 2:rK&fK, 3:rK&MK(1/K), 4:rK&MK 

        nhomo = 0   * 0 & 1: homogeneous, 2: kappa for branches, 3: N1, 4: N2
        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 1   * (0,1,2): rates (alpha>0) or ancestral states

   Small_Diff = 7e-6
    cleandata = 0 * remove sites with ambiguity data (1:yes, 0:no)?
*        icode = 0  * (with RateAncestor=1. try "GC" in data,model=4,Mgene=4)
  fix_blength = 1  * 0: ignore, -1: random, 1: initial, 2: fixed
       method = 1  * Optimization method 0: simultaneous; 1: one branch a time
"""
	output.close()








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


def get_downstream_terminal_nodes(treeObject, node):

	node_list=[]

	def add_successors(treeObject, node, node_list):
		daughters=treeObject.node(node).get_succ()
		
		for daughter in daughters:
			node_list=add_successors(treeObject, daughter, node_list)
		if treeObject.is_terminal(node):
			node_list.append(node)
		return node_list

	node_list=add_successors(treeObject, node, node_list)
	#print node_list
	return node_list


def get_upstream_terminal_nodes(treeObject, node):

	
	downstream_terminals=set(get_downstream_terminal_nodes(treeObject, node))
	all_terminals=set(treeObject.get_terminals())
	
	node_list=list(all_terminals.difference(downstream_terminals))
	
	#print node_list, downstream_terminals, all_terminals
	
	return node_list


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



def get_downstream_total_branch_lengths(treeObject, node):

	brlen=0.0

	def distance_to_daughters(node, brlen):
		daughters=treeObject.node(node).get_succ()
		
		for daughter in daughters:
			brlen+=treeObject.distance(node,daughter)
			brlen=distance_to_daughters(daughter, brlen)
		
		
		return brlen

	brlen=distance_to_daughters(node, brlen)
	
	return brlen





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

def get_all_nodes(treeObject):

	node_list=[]

	def add_successors(treeObject, node, node_list):
		daughters=treeObject.node(node).get_succ()
		
		for daughter in daughters:
			node_list=add_successors(treeObject, daughter, node_list)
		
		node_list.append(node)
		return node_list

	node_list=add_successors(treeObject, treeObject.root, node_list)
	
	return node_list


def get_nodes_not_downstream(treeObject, node):
	downstream_nodes=set(get_downstream_nodes(treeObject, node))

	all_nodes=set(get_downstream_nodes(treeObject, treeObject.root))
	
	return list(all_nodes.difference(downstream_nodes))


def get_common_ancestor_of_taxon_list(treeObject,taxonlist):
	
	if len(taxonlist)==0:
		return -1
	elif len(taxonlist)==1:
		return taxonlist[0]
	
	common_ancestor=taxonlist[0]
	
	for x in taxonlist[1:]:
		common_ancestor=treeObject.common_ancestor(common_ancestor,x)
	
	return common_ancestor



def colour_nodes_by_tree_position(treeObject, ladderize=None):
	node=treeObject.root
	
#	max_depth, max_right, max_left=get_max_depth(treeObject)
	#print max_depth, max_right, max_left
	def get_node_colour(node, num_nodes, node_num):
		
		if treeObject.node(node).get_succ():
			daughtersnew=ladderize_nodes(treeObject.node(node).get_succ(),ladderize=ladderize)
			daughtersnew.reverse()
			
			node_num=get_node_colour(daughtersnew[0], num_nodes, node_num)
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
			data=treeObject.node(node).get_data()
			#data.comment[Objecttype].seq=data.comment[Objecttype].seq.tomutable()
			
			
			if node==treeObject.root:
				data.comment["colour"]=(0, 0, 0)
			else:
				data.comment["colour"]=(red, green, blue)
			treeObject.node(node).set_data(data)
			
			node_num=get_node_colour(daughtersnew[1], num_nodes, node_num+1)
			
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
			data=treeObject.node(node).get_data()
			#data.comment[Objecttype].seq=data.comment[Objecttype].seq.tomutable()
			
			
			if node==treeObject.root:
				data.comment["colour"]=(0, 0, 0)
			else:
				data.comment["colour"]=(red, green, blue)
			treeObject.node(node).set_data(data)
			
			node_num+=1
		
		return node_num
	num_nodes=(treeObject.count_terminals(node)*2)-1
		
	#print max_distance, num_nodes
	
	nodecolours={}
	
	node_num=get_node_colour(node, num_nodes, 0)
	
	#print nodecolours
	
	return treeObject






def colour_nodes_by_tree_distances(treeObject, ladderize=None):
	
	max1, max2, max3, distance12, distance13, distance23=get_max_distances(treeObject)
	#print max_depth, max_right, max_left
	
	numterminals=treeObject.count_terminals(treeObject.root)
	
	def get_node_colour(node):
		
		
		for daughter in treeObject.node(node).get_succ():
			node_num=get_node_colour(daughter)

		
		rbCA=treeObject.common_ancestor(max1,node)
		if rbCA==0:
			rbCA=treeObject.common_ancestor(max2,node)
		if rbCA==0:
			rbCA=treeObject.common_ancestor(max1,max2)
		
		redbluedistance=((distance12)/2)-(treeObject.distance(max1,rbCA))
		
		if redbluedistance<0:
			negative=True
			redbluedistance=redbluedistance*-1
		else:
			negative=False
		
		redblueproportion=redbluedistance/((distance12)/2)
		#positive redbluedistance means closer to red, negative redbluedistance means closer to blue, zero means closer to magenta
		if redbluedistance>0:
			red=((1-(redbluedistance/((distance12)/2)))*125.5)
			blue=125.5
		else:
			blue=125.5+((1-(redbluedistance/((distance12)/2)))*125.5)
			red=125.5
			
		
		rgCA=treeObject.common_ancestor(max1,node)
		if rgCA==0:
			rgCA=treeObject.common_ancestor(max3,node)
		if rgCA==0:
			rgCA=treeObject.common_ancestor(max1,max3)
		
		redgreendistance=((distance13)/2)-treeObject.distance(max1,rgCA)
		#positive redgreendistance means closer to red, negative redgreendistance means closer to green, zero means closer to yellow
		if redgreendistance>0:
			red+=((1-(redgreendistance/((distance13)/2)))*125.5)
			green=125.5
		else:
			green=125.5+((1-(redgreendistance/((distance13)/2)))*125.5)
			red+=125.5
		
		bgCA=treeObject.common_ancestor(max2,node)
		if bgCA==0:
			bgCA=treeObject.common_ancestor(max3,node)
		if bgCA==0:
			bgCA=treeObject.common_ancestor(max2,max3)
		
		bluegreendistance=((distance23)/2)-treeObject.distance(max2,bgCA)
		#positive bluegreendistance means closer to blue, negative bluegreendistance means closer to green, zero means closer to cyan
		if bluegreendistance>0:
			blue+=((1-(bluegreendistance/((distance23)/2)))*125.5)
			green+=125.5
		else:
			green+=125.5+((1-(bluegreendistance/((distance23)/2)))*125.5)
			blue+=125.5
#		blue=(1-(float(treeObject.distance(max1,node))/max_distance1))*255
#		red=(1-(float(treeObject.distance(max2,node))/max_distance2))*255
#		print treeObject.distance(max2,node), treeObject.distance(max3,node), treeObject.distance(max1,node)
#		green=(1-(float(treeObject.distance(max3,node))/max_distance3))*255
##		green=(float(treeObject.count_terminals(node))/numterminals)*255
		
		
		#print (1-(redbluedistance/((distance12)/2))),redblueproportion, distance12, (1-(redgreendistance/((distance13)/2))), distance13, (1-(bluegreendistance/((distance23)/2))), distance23, node, red, green, blue
		
		data=treeObject.node(node).get_data()
		
		data.comment["colour"]=(red, green, blue)
		treeObject.node(node).set_data(data)

	node_num=get_node_colour(treeObject.root)
	return treeObject





def colour_nodes_by_splitting(treeObject):

# WORKS FOR ALL COLOURS (but would like to exclude yellows)
	colourrange=[0.0,1275.0]
	def number_to_colours(nodeproportion):
		
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
		return red, green, blue



	def colour_remaining_daughters(colournum, n):

		ds=treeObject.node(n).succ
		for d in ds:
			colour_remaining_daughters(colournum,d)
		data=treeObject.node(n).get_data()
		#data.comment={}
		if data.comment==None:
			data.comment={}
		data.comment["colour"]=number_to_colours(colournum)
		treeObject.node(n).set_data(data)
	
	
	def split_colours(node, colourrange):
		daughters=treeObject.node(node).get_succ()
		#print node, daughters
		if len(daughters)>0:
			nodebrlens=get_downstream_total_branch_lengths(treeObject, node)
			daughter1brlens=get_downstream_total_branch_lengths(treeObject, daughters[0])
			daughter2brlens=get_downstream_total_branch_lengths(treeObject, daughters[1])
			if len(daughters)==3:
				daughter3brlens=get_downstream_total_branch_lengths(treeObject, daughters[2])
			else:
				daughter3brlens=0
			
			if daughter1brlens==0:
				#print daughter1brlens, treeObject.node(daughters[0]).succ
				colour_remaining_daughters(colourrange[0],daughters[0])
				#print daughters[0], number_to_colours(colourrange[0])
				#data=treeObject.node(daughters[0]).get_data()
				##data.comment={}
				#data.comment["colour"]=number_to_colours(colourrange[0])
				#treeObject.node(daughters[0]).set_data(data)
			if daughter2brlens==0:
				#print daughter2brlens, treeObject.node(daughters[1]).succ
				colour_remaining_daughters(colourrange[1],daughters[1])
				#print daughters[1], number_to_colours(colourrange[1])
				#data=treeObject.node(daughters[1]).get_data()
				#data.comment={}
				#data.comment["colour"]=number_to_colours(colourrange[1])
				#treeObject.node(daughters[1]).set_data(data)
			if len(daughters)==3 and daughter3brlens==0:
				colour_remaining_daughters(colourrange[1],daughters[2])
			
			try:
				d1proportion=daughter1brlens/nodebrlens
			except ZeroDivisionError:
				if len(daughters)==3:
					d1proportion=0.33
				else:
					d1proportion=0.49
			try:
				d2proportion=daughter2brlens/nodebrlens
			except ZeroDivisionError:
				if len(daughters)==3:
					d2proportion=0.33
				else:
					d2proportion=0.49
			if len(daughters)==3:
				try:
					d3proportion=daughter3brlens/nodebrlens
				except ZeroDivisionError:
					d3proportion=0.33	
			
			#print node, daughter1brlens, daughter2brlens, nodebrlens
			#d1proportion=daughter1brlens/nodebrlens
			#d2proportion=daughter2brlens/nodebrlens
			colourrangelength=colourrange[1]-colourrange[0]
			
			
			if len(daughters)==3:
				d1colourrange=[colourrange[0], colourrange[0]+(colourrangelength*d1proportion)]
				d2diff=(1.0-(d1proportion+d3proportion))/2
				d2colourrange=[(colourrange[0]+(colourrangelength*d1proportion))+d2diff, (colourrange[0]+(colourrangelength*d1proportion)+(colourrangelength*d2proportion))+d2diff]
				d3colourrange=[colourrange[1]-(colourrangelength*d3proportion), colourrange[1]]
			else:
				d1colourrange=[colourrange[0], colourrange[0]+(colourrangelength*d1proportion)]
				d2colourrange=[colourrange[1]-(colourrangelength*d2proportion), colourrange[1]]	
			
			if len(daughters)==3:
				nodecolour=colourrange[0]+((colourrangelength*d1proportion)+(colourrangelength*d2proportion))+(3*((colourrangelength-(((colourrangelength*d1proportion)+(colourrangelength*d2proportion))+(colourrangelength*d3proportion)))/2))
			else:	
				nodecolour=colourrange[0]+(colourrangelength*d1proportion)+((colourrangelength-((colourrangelength*d1proportion)+(colourrangelength*d2proportion)))/2)
			data=treeObject.node(node).get_data()
			if data.comment==None:
				data.comment={}
			data.comment["colour"]=number_to_colours(nodecolour)
			treeObject.node(node).set_data(data)
#			if len(daughters)==3:
#				print node, daughters, nodebrlens, daughter1brlens, daughter2brlens, daughter3brlens, d1proportion, d2proportion, d3proportion, colourrange, d1colourrange, d2colourrange, d3colourrange, nodecolour, colourrangelength
#				sys.exit()
			
			split_colours(daughters[0], d1colourrange)
			split_colours(daughters[1], d2colourrange)
			if len(daughters)==3:
				split_colours(daughters[2], d3colourrange)
		
	

	split_colours(treeObject.root, colourrange)
	return treeObject








####################################################
# Function to print the colours of nodes to a file #
####################################################


def print_node_colours(treeObject, handle):
	
	def print_node_colour(node):
		#print node
		colour=treeObject.node(node).data.comment["colour"]
		if treeObject.is_terminal(node):
			printlist=[str(node), ' '.join([treeObject.node(node).data.taxon, str(int(colour[0])), str(int(colour[1])), str(int(colour[2]))]), str(RGBToHTMLColor(colour))]
		else:
			printlist=[str(node), ' '.join(["-", str(int(colour[0])), str(int(colour[1])), str(int(colour[2]))]), str(RGBToHTMLColor(colour))]
		print >> handle, '\t'.join(printlist)
		for daughter in treeObject.node(node).get_succ():
			print_node_colour(daughter)
	
	print >> handle, 'Node\ttaxon\tRGB\tHtml'
	
	print_node_colour(treeObject.root)

###############################################################################################
# Function to print a tree as a string. Fixes float problem in Bio.Nexus.Trees.Tree.to_string #
###############################################################################################


def tree_to_string(treeObject, support_as_branchlengths=False,branchlengths_only=False,plain=True,plain_newick=False,ladderize=None, collapse=False, cutoff=0.0, treename=False, comments=False, node_names=False):
	
	"""Return a paup compatible tree line.
	
	tree_to_string(treeObject,support_as_branchlengths=False,branchlengths_only=False,plain=True)
	"""
	# if there's a conflict in the arguments, we override plain=True
	
	if support_as_branchlengths or branchlengths_only:
		plain=False
	treeObject.support_as_branchlengths=support_as_branchlengths
	treeObject.branchlengths_only=branchlengths_only
	treeObject.plain=plain
	
	if node_names:
		comments=True
	
	def get_comments_line(data, terminal=False):
		commentlist=[]
		if node_names:
			if not terminal and data.taxon!=None:
				return data.taxon
		else:
			if data.comment!=None and type(data.comment) is dict:
				for commenttype in data.comment:
					if type(data.comment[commenttype]) in [list, tuple]:
						commentlist.append("[&colour="+' '.join(map(str,data.comment[commenttype]))+"]")
					elif type(data.comment[commenttype]) in [str, int]:
						commentlist.append("[&colour="+str(data.comment[commenttype])+"]")
		
		return ''.join(commentlist)
	
	def make_info_string(data,terminal=False):
		"""Creates nicely formatted support/branchlengths."""
		if comments:
			commentsline=get_comments_line(data, terminal=terminal)
		else:
			commentsline=''
		# CHECK FORMATTING
		if treeObject.plain: # plain tree only. That's easy.
			return ''
		
		elif treeObject.support_as_branchlengths: # support as branchlengths (eg. PAUP), ignore actual branchlengths
			
			if terminal: # terminal branches have 100% support
				if collapse and treeObject.max_support<cutoff:
					brlength=0.0
				else:
					brlength=treeObject.max_support
				return '%s:%s' % (commentsline, str(brlength))
			else:
				if collapse and data.support<cutoff:
					brlength=0.0
				else:
					brlength=data.support
				return '%s:%s' % (commentsline, str(brlength))
		elif treeObject.branchlengths_only: # write only branchlengths, ignore support
			if collapse and data.branchlength<cutoff:
				brlength=0.0
			else:
				brlength=data.branchlength
			return '%s:%s' % (commentsline, str(brlength))
		else: # write suport and branchlengths (e.g. .con tree of mrbayes)
			
			if terminal:
				if data.branchlength is not None:
					if collapse and data.branchlength<cutoff:
						brlength=0.0
					else:
						brlength=data.branchlength
					return '%s:%s' % (commentsline, str(brlength))
				else:
					return '%s:0' % (commentsline)
			else:
				if data.branchlength is not None and data.support is not None: # we have blen and suppport
					if collapse and data.branchlength<cutoff:
						brlength=0.0
					else:
						brlength=data.branchlength
					return '%s%s:%s' % (commentsline, str(data.support),str(brlength))
				elif data.branchlength is not None: # we have only blen
					if collapse and data.branchlength<cutoff:
						brlength=0.0
					else:
						brlength=data.branchlength
					return '%s0:%s' % (commentsline, str(brlength))
				elif data.support is not None: # we have only support
					#return ':0'
					return '%s%s:0' % (commentsline, str(data.support))
					
				else:
					return '%s:0' % commentsline
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
   
#		if not treeObject.node(node).get_succ():	#if branch is terminal
		if treeObject.is_terminal(node):
			
			return treeObject.node(node).data.taxon+make_info_string(treeObject.node(node).data,terminal=True) 
		else:
			
			succnodes=ladderize_nodes(treeObject.node(node).get_succ(),ladderize=ladderize) 
			
			subtrees=[newickize(sn,ladderize=ladderize) for sn in succnodes]
			
			newsubtrees=[]
			
			for subtree in subtrees:
				if subtree[0]=="(":
					newsubtrees.append(subtree)
				else:
					for bit in subtree.split(","):
						newsubtrees.append(bit)
			newsubtrees.sort()
			
			#print subtrees, newsubtrees
			#print treeObject.node(node).data.branchlength, cutoff, treeObject.node(node).data.branchlength>cutoff
			if not collapse:
				return '(%s)%s' % (','.join(subtrees),make_info_string(treeObject.node(node).data))
			elif treeObject.node(node).data.branchlength is None or treeObject.node(node).data.branchlength>cutoff:
				
				#partstring="".join(["("]*(len(newsubtrees)-1))
				partstring="("
				partstring=partstring+newsubtrees[0]+','
				#print newsubtrees
				for newsubtree in newsubtrees[1:]:
					partstring=partstring+newsubtree+","
				if partstring[-1]==",":
					partstring=partstring[:-1]
				#partstring=":".join(partstring.split(":")[:-1])
				
				partstring=partstring+")"
				#print partstring
#				sys.exit()
				#print '%s%s' % (partstring,make_info_string(treeObject.node(node).data))
				return '%s%s' % (partstring,make_info_string(treeObject.node(node).data))
			else:
				#return '%s,%s' % (','.join(subtrees),make_info_string(treeObject.node(node).data))
				return '%s' % (','.join(newsubtrees))
	treeline=[]
	if treename:		   
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
	
	subtrees=[newickize(sn,ladderize=ladderize) for sn in succnodes] 
	treeline.append('(%s)' % ','.join(subtrees)) 
	if plain_newick: 
		return treeline[-1] 
	else: 
		return ' '.join(treeline)+';'
		





###############################################################################################
# Function to print a tree as a string. Fixes float problem in Bio.Nexus.Trees.Tree.to_string #
###############################################################################################


def tree_to_figtree_string(treeObject, support_as_branchlengths=False,branchlengths_only=False,plain=True,plain_newick=False,ladderize=None):

	"""Return a figtree compatible tree line.
	
	tree_to_string(treeObject,support_as_branchlengths=False,branchlengths_only=False,plain=True)
	"""
	# if there's a conflict in the arguments, we override plain=True
	if support_as_branchlengths or branchlengths_only:
		plain=False
	treeObject.support_as_branchlengths=support_as_branchlengths
	treeObject.branchlengths_only=branchlengths_only
	treeObject.plain=plain

	def make_info_string(node, data,terminal=False):
		"""Creates nicely formatted support/branchlengths."""
		node_data=treeObject.node(node).get_data().comment
		if node_data.has_key("colour"):
			colourbit="[&!color=%s,label=0.0]" % (RGBToHTMLColor(node_data["colour"]))
		else:
			colourbit=""
		
		
		
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

   
	def newickize(node,ladderize=None): 
		"""Convert a node tree to a newick tree recursively.""" 
   
		if not treeObject.node(node).get_succ():	#terminal 
			return treeObject.node(node).data.taxon+make_info_string(node, treeObject.node(node).data,terminal=True) 
		else: 
			succnodes=ladderize_nodes(treeObject.node(node).get_succ(),ladderize=ladderize) 
			subtrees=[newickize(sn,ladderize=ladderize) for sn in succnodes]
			return '(%s)%s' % (','.join(subtrees),make_info_string(node, treeObject.node(node).data))
	
	
	treeline=['#NEXUS\n']
	treeline.append('begin taxa;\n')
	
	treeline.append('\tdimensions ntax='+str(len(treeObject.get_terminals()))+';\n')
	treeline.append('\ttaxlabels\n')
	for taxon in treeObject.get_terminals():
		
		node_data=treeObject.node(taxon).get_data().comment
		if node_data.has_key("colour"):
			treeline.append('\t'+treeObject.node(taxon).get_data().taxon+'[&!color='+RGBToHTMLColor(node_data["colour"])+']\n')
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
	subtrees=[newickize(sn,ladderize=ladderize) for sn in succnodes] 
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
#	max_distance=[0, 0, 0]
#	
#	for x, terminal in enumerate(terminal_nodes):
#		for secondterminal in terminal_nodes[x+1:]:
#			terminaldistance=treeObject.distance(secondterminal, terminal)
#			if float(terminaldistance)>max_distance[0]:
#				max_distance=[terminaldistance, terminal, secondterminal]
#			#distances.append([float(treeObject.distance(secondterminal, terminal)), terminal, secondterminal])
#	
#	print max_distance
#	
	max_distance=[treeObject.distance(terminal_nodes[0], terminal_nodes[1]), terminal_nodes[0], terminal_nodes[1]]
	for x, newterminal in enumerate(terminal_nodes[2:]):
		for terminal in max_distance[1:]:
			terminaldistance=treeObject.distance(newterminal, terminal)
			if float(terminaldistance)>max_distance[0]:
				max_distance=[terminaldistance, terminal, newterminal]
	
#	print max_distance
	
	
	
	
	
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
	
	
	
	#Calculate which branch is the halfway point between the two most distant nodes
	outgroup_brlen=max_distance[0]
	ingroup_brlen=max_distance[0]
	
	for node in nodes_joining_terminals:
		if (max_distance[0]/2)-float(treeObject.distance(node, max_distance[1]))>0:
			if ((max_distance[0]/2)-float(treeObject.distance(node, max_distance[1])))<outgroup_brlen:
				outgroup_node=node
				outgroup_brlen=((max_distance[0]/2)-float(treeObject.distance(node, max_distance[1])))
				
					 
		else:
			if (float(treeObject.distance(node, max_distance[1]))-(max_distance[0]/2))<ingroup_brlen:
				ingroup_node=node
				ingroup_brlen=(float(treeObject.distance(node, max_distance[1]))-(max_distance[0]/2))
				
	
	#identify the terminal-most node on the branch and make it the root node
	
	if outgroup_node in treeObject.node(ingroup_node).succ:
		root_node=outgroup_node
	else:
		root_node=ingroup_node
	
	#root the tree on an outgroup composed of all terminals upstream of the root node
	treeObject.root_with_outgroup(treeObject.get_taxa(root_node))
	
	
	#reset the branchlengths to the two nodes from the new root
	outgroup_data=treeObject.node(outgroup_node).get_data()
	outgroup_data.branchlength=outgroup_brlen
	treeObject.node(outgroup_node).set_data(outgroup_data)
	ingroup_data=treeObject.node(ingroup_node).get_data()
	ingroup_data.branchlength=ingroup_brlen
	treeObject.node(ingroup_node).set_data(ingroup_data)
	
	
	return treeObject

	



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




###################################################################
# Function to create a SNP alignment #
###################################################################

def Create_SNP_alignment(alignment, SNPlocations):
	
	alphabet = Gapped(IUPAC.unambiguous_dna)

	SNPalignment = Alignment(alphabet)

	
	for record in alignment:
		SNPseq=""
		
		for base in SNPlocations:
			SNPseq=SNPseq+record.seq[base].replace("N","?")
		
		SNPalignment.add_sequence(record.id, SNPseq)
		
	
	return SNPalignment






##############################################################
# Function to parsimoneously reconstruct all sites on a tree #
##############################################################
	
	
def parsimonious_sequence_reconstruction(treeObject, alignmentObject, transformation="acctran", sequence_Objecttype="sequence", locations=[], genetic_code_number=1):

	print "Reconstructing SNPs using parsimony"
	print "Using genetic code number", genetic_code_number

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
		#print node, treeObject.node(node).get_data().comment["sequence"][sitenumber]
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
				
				treeObject=sankoff(treeObject, sitenumber, daughter)
			 	
			 	
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
#			if len(min_cost_states)>1 and "N" in min_cost_states:
#				min_cost_states.remove("N")
			
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
#			if len(min_cost_states)>1 and "N" in min_cost_states:
#				min_cost_states.remove("N")
			
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
	
	
	
	
	
	def fix_indels_in_paml_sequences(treeObject, sitenumber, sequence_Objecttype="sequence"):
		
		
		def indel_second_sankoff_traversal(treeObject, sitenumber, node=-1, transformation="acctran", sequence_Objecttype="sequence"):
		
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
				
				if len(min_cost_states)==1 and "-" in min_cost_states:
					
				
					node_data=treeObject.node(node).get_data()
				
					
					node_data.comment["sequence"][sitenumber]="-"
					
					treeObject.node(node).set_data(node_data)
#				if len(min_cost_states)>1 and "N" in min_cost_states:
#					min_cost_states.remove("N")
		
		
		
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
				
				if len(min_cost_states)==1 and "-" in min_cost_states:
					daughter_data.comment["sequence"][sitenumber]="-"
					
					treeObject.node(daughter).set_data(daughter_data)
					
#				if len(min_cost_states)>1 and "N" in min_cost_states:
#					min_cost_states.remove("N")
				
				
				
				treeObject=indel_second_sankoff_traversal(treeObject, sitenumber, daughter, transformation=transformation)
	
			if node==treeObject.node(treeObject.root).succ[0]:
				rootdata=treeObject.node(treeObject.root).get_data()
				rootdata.comment["sequence"][sitenumber]=treeObject.node(node).get_data().comment["sequence"][sitenumber]
				treeObject.node(treeObject.root).set_data(rootdata)
			
		
		
			return treeObject
			
			
			
		treeObject=sankoff(treeObject, sitenumber)
		treeObject=indel_second_sankoff_traversal(treeObject, sitenumber, treeObject.root, transformation=transformation)
		
		return treeObject
	
	
	
	def put_paml_sequences_on_tree(treeObject, paml_tree, paml_alignment, SNPlocations, sequence_Objecttype="sequence"):
		
		
		node=treeObject.root
		
		def replace_sequences_with_paml_reconstruction(node):
			
			daughters=treeObject.node(node).get_succ()
			
			for daughter in daughters:
				
				replace_sequences_with_paml_reconstruction(daughter)
			
			if treeObject.is_internal(node) and node!=treeObject.root:
				
				
				for record in paml_alignment:
					if record.id==str(int(paml_tree.node(node).data.branchlength)):
						
						data=treeObject.node(node).get_data()
						mutable_seq=data.comment[sequence_Objecttype]
								
						for pamllocation, record_base in enumerate(record.seq):
							location=SNPlocations[pamllocation]
							
							mutable_seq[location]=record_base
							
							 	
						
						data.comment[sequence_Objecttype]=mutable_seq
						treeObject.node(node).set_data(data)
			if node==treeObject.node(treeObject.root).succ[0]:
				rootdata=treeObject.node(treeObject.root).get_data()
				rootdata.comment[sequence_Objecttype]=treeObject.node(node).get_data().comment[sequence_Objecttype][:]
				treeObject.node(treeObject.root).set_data(rootdata)
	
		replace_sequences_with_paml_reconstruction(node)
	
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
			
			
			#print node, daughter, node_state, daughter_state
			#print node, node_ambiguity, daughter, daughter_ambiguity, new_daughter_states, daughter_states
			
			#if node_state=="-" or daughter_state=="-":
				#print "here"
			
			if daughter_state!=node_state and daughter_state!="N" and node_state!="N":
#				if daughter==6 or daughter==5:
#					print sitenumber+1, "from node", node, "to node", daughter, node_state, daughter_state
				if node_state=="-" and daughter_state!="-":
					#print "insertion at base", sitenumber+1, "from node", node, "to node", daughter#, node_ambiguity, daughter_ambiguity
					daughter_data=treeObject.node(daughter).get_data()
					if not daughter_data.comment.has_key("insertion_locations"):
						daughter_data.comment["insertion_locations"]=[]
					daughter_data.comment["insertion_locations"].append(sitenumber)
					treeObject.node(daughter).set_data(daughter_data)
				elif node_state!="-" and daughter_state=="-":
					#print "deletion at base", sitenumber+1, "from node", node, "to node", daughter#, node_ambiguity, daughter_ambiguity
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
					SNP.originalSNP=True
					SNP.genetic_code_number=genetic_code_number
					
					
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
		for x, base in enumerate(mutable_seq):
			if base=="N" and treeObject.node(rootnode).get_data().comment[sequence_Objecttype][x]!="N":
				mutable_seq[x]= treeObject.node(rootnode).get_data().comment[sequence_Objecttype][x]
		data.comment[sequence_Objecttype]=mutable_seq
		treeObject.node(node).set_data(data)

	
	#treeObject.display()
	
	if locations==[]:
		locations=xrange(alignmentObject.get_alignment_length())
	
	#print alignmentObject
	if transformation in ["acctran", "deltran"]:
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
			column=alignmentObject[:,columnnumber]
			
			column=column.replace("N","").replace("?","").replace("X","")
			firstbase=column[0]
			for base in column[1:]:
				if base!=firstbase:#can add a speedup here by only reconstructing informative sites?
					
					treeObject=sankoff(treeObject, columnnumber, rootnode)
					treeObject=sankoff_second_traversal(treeObject, columnnumber, rootnode, transformation=transformation)
					treeObject=call_snps(treeObject, rootnode, columnnumber, sequence_Objecttype=sequence_Objecttype, transformation=transformation)
					break
	
	elif transformation=="ML":
		#Make a random name for the PAML run
		chars = string.ascii_letters + string.digits
		tmpname='tmp'+"".join(choice(chars) for x in xrange(randint(8, 10)))
		
		#create temporary alignment and tree for paml
		
		SNPalignment=Create_SNP_alignment(alignmentObject, locations)
		
		#print a phylip file of the SNP alignment
	
		sequencenames={}
		convertnameback={}
		seqnametoindex={}
		
		handle = open(tmpname+".phy", "w")
		print >> handle, len(SNPalignment), SNPalignment.get_alignment_length()
		count=1
		for record in SNPalignment:
	
			name="seq"+str(count)
			
			sequencenames[name]=record.id
			convertnameback[record.id]=name
			seqnametoindex[name]=count-1
			
			print >> handle, name+"  "+record.seq
	# add this back in to split the alignment into blocks of 60 bases
	#		for f in range(0,len(record.seq),60):
	#			if f+60< len(record.seq):
	#				print >> handle, record.seq[f:f+60]
	#			else:
	#				print >> handle, record.seq[f:]
				
			count=count+1
		handle.close()
		
		treestring= tree_to_string(treeObject, False, True, True, True)
		for name in sequencenames.keys():
			treestring=treestring.replace(sequencenames[name]+":", name+":")
		if treestring[-1]!=";":
			treestring=treestring=treestring+";"
		treestring="("+"(".join(treestring.split("(")[1:])
		handle = open(tmpname+".tre", "w")
		print >> handle, treestring
		handle.close()
		
		#create baseml control file for paml
		
		print "Running PAML to reconstruct ancestral states"
		sys.stdout.flush()
	
		create_baseml_control_file(tmpname+".phy", tmpname+".tre", 1)
		
		#run paml
		
		os.system("baseml > "+tmpname+"temp.tmp")
	
		#remove spaces from rst alignment (necessary to allow easier reading of tree and ancestral sequences
		
		os.system("sed 's/node #//g' rst > rstnew")
		os.system('grep -v "^$" rstnew > rstnew2')
		
		#os.system('cat rates')
		
		#extract the tree with all nodes numbered from PAML rst output file
		
		print "Reading PAML tree"
		sys.stdout.flush()
		
		pamltreefile=os.popen('grep -A 1 "tree with node labels for Rod Page\'s TreeView" rstnew | tail -n 1')
		
		pamltreestring=pamltreefile.read().strip()
		
		
		#read paml tree into memory (note that node ids will be stored as branchlengths)
		
		pamltree=Trees.Tree(pamltreestring, rooted=True)
		
		
		for terminal in pamltree.get_terminals():
			terminal_data=pamltree.node(terminal).get_data()
			
			#print sequencenames[terminal_data.taxon.split("_")[1]]
			
			terminal_data.taxon=sequencenames[terminal_data.taxon.split("_")[1]]
			
			pamltree.node(terminal).set_data(terminal_data)
		
		
		#convert branchlengths (node ids) into support attributes
		
		#pamltree.branchlength2support()
		
		#get the root node number from the paml tree (I think this might always be 0)
		
		rootnode=pamltree.root
		
		#extract alignment PAML rst output file
		
		print "Reading ancestral sequences from PAML"
		sys.stdout.flush()
		
		pamlalignstats=os.popen('grep -A 1 "List of extant and reconstructed sequences" rstnew2 | tail -n 1')
		
		try:
			n=int(pamlalignstats.read().split()[0])
		except ValueError:
			DoError("PAML analysis failed")
		
		pamlalignfile=os.popen('grep -A '+str(n+1)+' "List of extant and reconstructed sequences" rstnew2 | tail -n '+str(n+1))
		
		pamlalignment = AlignIO.read(pamlalignfile, "phylip")
		
		#Now we have the ancestral state reconstructions and tree in memory we need to make the input for the recombination detection program
		
		#first create a dictionary of all paml snp sequences
		
		pamlsequences={}
		for record in pamlalignment:
			pamlsequences[record.id]=record.seq
			name="seq"+str(count)
			
			if sequencenames.has_key(record.id):
				record.id=sequencenames[record.id]
				
			
			
		
		
		treeObject=put_paml_sequences_on_tree(treeObject, pamltree, pamlalignment, locations, sequence_Objecttype=sequence_Objecttype)
		
		print "Fixing indels and calling SNPs"
		sys.stdout.flush()

		
		
		#print pamlalignment
		
		os.system("rm baseml.ctl rub rstnew2 rstnew rst1 rst rates mlb lnf 2base.t "+tmpname+"*")
		
		rootnode=treeObject.root
		#print len(locations)
		
		for columnnumber in locations:
			#when we get the new version of biopython, this can be changed to slice notation alignmentObject[:,base]
			column=alignmentObject[:,columnnumber]
#			if columnnumber>500:
#				break
			column=column.replace("N","").replace("?","").replace("X","")
			if "-" in column and (len(column)-1)>len(column.replace("-","")):
				fix_indels_in_paml_sequences(treeObject, columnnumber, sequence_Objecttype="sequence")
			
			firstbase=column[0]
			for base in column[1:]:
				if base!=firstbase:#can add a speedup here by only reconstructing informative sites?
					
					treeObject=call_snps(treeObject, rootnode, columnnumber, sequence_Objecttype=sequence_Objecttype, transformation=transformation)
					break
		
		
		#sys.exit()
	
	
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
				SNP.originalSNP=True
				
				#({"parent":node, "daughter":daughter, "from":ambiguity_to_bases[node_data.comment[sequence_Objecttype][location]][:], "to":ambiguity_to_bases[node_data.comment[sequence_Objecttype][location]][:]})
			
		
		
	return SNP_locations
		
	








	
################################################################
# Function to make branchlength the number of SNPs on a branch #
################################################################
	
	
def branchlengths_to_SNP_count(treeObject, node=-1, lengthtype="SNP_locations", SNP_type=""):
#	print SNP_type
	if node==-1:
		node=treeObject.root
	
	daughters=treeObject.node(node).get_succ()
		
	for daughter in daughters:
		treeObject=branchlengths_to_SNP_count(treeObject, daughter, lengthtype, SNP_type)	
	
	data=treeObject.node(node).get_data()
	if data.comment.has_key(lengthtype):
		SNPcount=0
		for SNP in data.comment[lengthtype]:
		
			if SNP_type=="":
				if not data.comment[lengthtype][SNP].recombination:
					SNPcount+=1
			
			else:
#				print data.comment[lengthtype][SNP].SNP_type, data.comment[lengthtype][SNP].codon_type, SNPcount,
				if not data.comment[lengthtype][SNP].recombination and data.comment[lengthtype][SNP].codon_type==SNP_type:
					SNPcount+=1
#				print SNPcount
				
				
				
				
		data.branchlength=SNPcount
		#data.branchlength=len(data.comment[lengthtype])
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
	tmpname='tmp'+"".join(choice(chars) for x in xrange(randint(8, 10)))

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
		sys.stdout.flush()
		os.system(RAXML+" -f d -s "+tmpname+".phy -m "+model+" -n "+tmpname+" > /dev/null")
		outputname="RAxML_result."+tmpname
	else:
		print "Running RAxML phylogeny with "+model+" model of evolution and "+str(bootstrap)+" bootstrap replicates..."
		os.system(RAXML+" -f a -x "+str(randrange(1,99999))+" -p "+str(randrange(1,99999))+" -# "+str(bootstrap)+" -m "+model+" -s "+tmpname+".phy -n "+tmpname+" > /dev/null")
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
	tmpname='tmp'+"".join(choice(chars) for x in xrange(randint(8, 10)))

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
	
	for site in xrange(alignmentObject.get_alignment_length()):
		if alignmentObject[refnum][site]!="-":
			ref_to_alignment[refpos]=site
			alignment_to_ref[site]=refpos
			refpos+=1
#		else:
#			print site
#			sys.exit()
	
	return ref_to_alignment, alignment_to_ref
	
	
	

	
	
def change_reference_location_to_alignment_location(seqFeature, translation_dict):

	#newFeature=copy.deepcopy(seqFeature)
	
	#print "before",  seqFeature.location
	new_Feature_locations=[]
	for x, part in enumerate(seqFeature.location.parts):
	
		startoffset=translation_dict[int(part.nofuzzy_start)]-int(part.nofuzzy_start)
		endoffset=(translation_dict[int(part.nofuzzy_end)-1]-int(part.nofuzzy_end)+1)
		strand=part.strand
		new_Feature_locations.append(FeatureLocation(start = part.start._shift(startoffset),end = part.end._shift(endoffset), strand=strand))
		
	#print "during", seqFeature.location, seqFeature.location.parts
	new_locs=new_Feature_locations[0]
	for x in new_Feature_locations[1:]:
		new_locs=new_locs+x
	seqFeature.location=new_locs			
	#print "after", seqFeature.location, seqFeature.location.parts
	#what about fuzzy positions???
	
	
	return seqFeature
	




def annotate_SNPs(treeObject, node=-1):

	if node==-1:
		node=treeObject.root
		
	daughters=treeObject.node(node).get_succ()
	
	node_data=treeObject.node(node).get_data()
	
	for daughter in daughters:
		treeObject=annotate_SNPs(treeObject, daughter)
		
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
					
				
				if 'annotation' in node_data.comment and  'annotation' in daughter_data.comment:
					daughter_data.comment["SNP_locations"][SNPlocation].get_annotation_info(node_data.comment['annotation'], daughter_data.comment['annotation'], treeObject.node(node).get_data().comment["sequence"], treeObject.node(daughter).get_data().comment["sequence"], node_node, daughter_node)
					treeObject.node(daughter).set_data(daughter_data)
			
		
		
		
	
	return treeObject








def dNdS_per_branch(treeObject, node=-1, genetic_code_number=1):

	if node==-1:
		node=treeObject.root
	
	
	
	daughters=treeObject.node(node).get_succ()
	node_data=treeObject.node(node).get_data()
	#node_revcompseq=Si_SNPs_temp.revcomp(str(node_data.comment["sequence"][:]))
	for daughter in daughters:
		
		treeObject=dNdS_per_branch(treeObject, node=daughter, genetic_code_number=genetic_code_number)
		
		daughter_data=treeObject.node(daughter).get_data()

		#daughter_revcompseq=Si_SNPs_temp.revcomp(str(daughter_data.comment["sequence"][:]))
		print node, "->", daughter
		sys.stdout.flush()
		for x, feature in enumerate(node_data.comment["annotation"]):
			
			daughter_feature=daughter_data.comment["annotation"][x]
			
			#if feature["frameshift"] or daughter_feature["frameshift"]:
				#print "gene is not original length in node or daughter... skipping"
				#print feature
				#continue
			
			if feature["strand"]==1:
				
				node_seq=str(node_data.comment["sequence"][feature["location"][0]:feature["location"][1]])
				daughter_seq=str(daughter_data.comment["sequence"][daughter_feature["location"][0]:daughter_feature["location"][1]])
				
				#if "-" in node_seq or "-" in daughter_seq:
					#continue
				
				try:
					dndsout=Si_SNPs_temp.dnbyds(node_seq,daughter_seq, genetic_code_number=genetic_code_number)
				except ValueError:
					daughter_data.comment["annotation"][x]["frameshift"]=True
					#print daughter_data.comment["annotation"][x]
				except ZeroDivisionError:
					node_data.comment["annotation"][x]["frameshift"]=True
					#print node_data.comment["annotation"][x]
				except FloatingPointError:
					continue
				else:
					daughter_data.comment["annotation"][x]["dNdS"]=dndsout
				
				#daughter_data.comment["annotation"][x]["dNdS"]=Si_SNPs_temp.dnbyds(node_seq,daughter_seq)
				
		
			else:
				
				node_seqlen=len(str(node_data.comment["sequence"]))
				daughter_seqlen=len(str(daughter_data.comment["sequence"]))
			
				#node_start, node_end=Si_SNPs_temp.find_gene_limits(Si_SNPs_temp.revcomp(str(treeObject.node(node).get_data().comment["sequence"])),  seqlen-feature.location.nofuzzy_end,seqlen-feature.location.nofuzzy_start, 1)#need to work out start in revcomp
				
				#node_seq=node_revcompseq[node_seqlen-feature["location"][1]:node_seqlen-feature["location"][0]]
				node_seq=Si_SNPs_temp.revcomp(str(node_data.comment["sequence"][feature["location"][0]:feature["location"][1]]))
				#daughter_seq=daughter_revcompseq[daughter_seqlen-daughter_feature["location"][1]:daughter_seqlen-daughter_feature["location"][0]]
				daughter_seq=Si_SNPs_temp.revcomp(str(daughter_data.comment["sequence"][daughter_feature["location"][0]:daughter_feature["location"][1]]))
	#			if "-" in node_seq or "-" in daughter_seq:
		#			continue
				

				try:
					dndsout=Si_SNPs_temp.dnbyds(node_seq,daughter_seq, genetic_code_number=genetic_code_number)
				except ValueError:
					daughter_data.comment["annotation"][x]["frameshift"]=True
					#print daughter_data.comment["annotation"][x]
				except ZeroDivisionError:
					node_data.comment["annotation"][x]["frameshift"]=True
					#print node_data.comment["annotation"][x]
				except FloatingPointError:
					continue
				else:
					daughter_data.comment["annotation"][x]["dNdS"]=dndsout

#				 Si_SNPs_temp.dnbyds(Si_SNPs_temp.revcomp(str(treeObject.node(node).get_data().comment["sequence"][feature.location.nofuzzy_start:feature.location.nofuzzy_end])),Si_SNPs_temp.revcomp(str(treeObject.node(daughter).get_data().comment["sequence"][feature.location.nofuzzy_start:feature.location.nofuzzy_end])))
		treeObject.node(daughter).set_data(daughter_data)	
		treeObject.node(node).set_data(node_data)
	
	return treeObject
	
	

def branch_dnds(treeObject, node=-1, genetic_code_number=1):

	if node==-1:
		node=treeObject.root
	
	
	
	daughters=treeObject.node(node).get_succ()
	node_data=treeObject.node(node).get_data()
	node_revcompseq=Si_SNPs_temp.revcomp(str(node_data.comment["sequence"][:]))
	for daughter in daughters:
		
		treeObject=branch_dnds(treeObject, node=daughter, genetic_code_number=genetic_code_number)
		
		daughter_data=treeObject.node(daughter).get_data()

		daughter_revcompseq=Si_SNPs_temp.revcomp(str(daughter_data.comment["sequence"][:]))
		print node, "->", daughter,
		
		total_node_seq=[]
		total_daughter_seq=[]
		for x, feature in enumerate(node_data.comment["annotation"]):
			
			daughter_feature=daughter_data.comment["annotation"][x]
			
			if feature["frameshift"] or daughter_feature["frameshift"]:
				#print "gene is not original length in node or daughter... skipping"
				continue
			
			if feature["strand"]==1:
				
				node_seq=str(node_data.comment["sequence"][feature["location"][0]:feature["location"][1]])
				daughter_seq=str(daughter_data.comment["sequence"][daughter_feature["location"][0]:daughter_feature["location"][1]])
				
				if "-" in node_seq or "-" in daughter_seq:
					continue
				
				total_node_seq.append(node_seq)
				total_daughter_seq.append(daughter_seq)
				
				#daughter_data.comment["annotation"][x]["dNdS"]=Si_SNPs_temp.dnbyds(node_seq,daughter_seq)
				
		
			else:
				
				node_seqlen=len(str(node_data.comment["sequence"]))
				daughter_seqlen=len(str(daughter_data.comment["sequence"]))
			
				#node_start, node_end=Si_SNPs_temp.find_gene_limits(Si_SNPs_temp.revcomp(str(treeObject.node(node).get_data().comment["sequence"])),  seqlen-feature.location.nofuzzy_end,seqlen-feature.location.nofuzzy_start, 1)#need to work out start in revcomp
				
				node_seq=node_revcompseq[node_seqlen-feature["location"][1]:node_seqlen-feature["location"][0]]
				daughter_seq=daughter_revcompseq[daughter_seqlen-daughter_feature["location"][1]:daughter_seqlen-daughter_feature["location"][0]]
				
				if "-" in node_seq or "-" in daughter_seq:
					continue
				total_node_seq.append(node_seq)
				total_daughter_seq.append(daughter_seq)


#				 Si_SNPs_temp.dnbyds(Si_SNPs_temp.revcomp(str(treeObject.node(node).get_data().comment["sequence"][feature.location.nofuzzy_start:feature.location.nofuzzy_end])),Si_SNPs_temp.revcomp(str(treeObject.node(daughter).get_data().comment["sequence"][feature.location.nofuzzy_start:feature.location.nofuzzy_end])))
		try:
			dndsout=Si_SNPs_temp.dnbyds(''.join(total_node_seq),''.join(total_daughter_seq), genetic_code_number=genetic_code_number)
			
			
			N=dndsout["N"]
			S=dndsout["S"]
			pN=dndsout["pN"]
			pS=dndsout["pS"]
			Nd=dndsout["Nd"]
			Sd=dndsout["Sd"]
			CDSlen=dndsout["CDSlen"]
			
			if N!=0:
				pN=Nd/N
			if S!=0:
				pS=Sd/S
			
			
			if pS==0:
				#print "No sites are synonymous."
				
				dN=0.0
				dS=0.0
				dNdS="-"
				z="-"
				
			elif pS<0.75 and pN<0.75:
				dS=(-3*(math.log(1 - ((pS*4)/3))))/4
				dN=(-3*(math.log(1 - ((pN*4)/3))))/4
				if dN==-0.0:
					dN=0.0
				varianceS=(9 * pS * (1 - pS))/(((3 - (4 *pS)) **2) * CDSlen);
				varianceN=(9 * pN * (1 - pN))/(((3 - (4 *pN)) **2) * CDSlen);
				z=(dN - dS) / math.sqrt(varianceS + varianceN)
				dNdS=dN/dS
			else:
				dN=pN
				dS=pS
				dNdS=pN/pS
				z="-"
			print dN, dS, dNdS, z
		except ValueError:
			print "Error calculating dN/dS"
			sys.stdout.flush()
		treeObject.node(daughter).set_data(daughter_data)	
		treeObject.node(node).set_data(node_data)
	
	return treeObject
	

	
def apply_annotation_to_root(treeObject, annotationObject, genetic_code_number=1):

	node=treeObject.root
		
	
	node_data=treeObject.node(node).get_data()

		
	seqlen=len(str(node_data.comment["sequence"]))
	revcompseq=Si_SNPs_temp.revcomp(str(node_data.comment["sequence"][:]))
		
	if not node_data.comment.has_key("annotation"):
		node_data.comment["annotation"]=[]

					
	for feature in annotationObject.features:
		
		if feature.location.strand==1:
			
			node_start, node_end=Si_SNPs_temp.find_gene_limits(str(node_data.comment["sequence"][feature.location.nofuzzy_start:feature.location.nofuzzy_end+9999]),0,feature.location.nofuzzy_end-feature.location.nofuzzy_start, genetic_code_number=genetic_code_number)
			
			#print "f", feature.location.nofuzzy_start, feature.location.nofuzzy_end, feature.location.nofuzzy_end-feature.location.nofuzzy_start
			
			node_start=node_start+feature.location.nofuzzy_start
			node_end=node_end+feature.location.nofuzzy_start
			
#			if feature.location.nofuzzy_end!=node_end:
#				print "f", feature.location.nofuzzy_start, feature.location.nofuzzy_end, feature.location.nofuzzy_end-feature.location.nofuzzy_start, node_start, node_end

		else:
			
			node_start, node_end=Si_SNPs_temp.find_gene_limits(revcompseq[seqlen-feature.location.nofuzzy_end:(seqlen-feature.location.nofuzzy_start)+9999],  0,feature.location.nofuzzy_end-feature.location.nofuzzy_start, genetic_code_number=genetic_code_number)
			
			#print "r", feature.location.nofuzzy_start, feature.location.nofuzzy_end, feature.location.nofuzzy_end-feature.location.nofuzzy_start
			
			node_start=feature.location.nofuzzy_end-node_start
			node_end=feature.location.nofuzzy_end-node_end
			
			node_start_new=node_end
			node_end=node_start
			node_start=node_start_new
			
#			if feature.location.nofuzzy_end!=node_end:
#				print "r", feature.location.nofuzzy_start, feature.location.nofuzzy_end, feature.location.nofuzzy_end-feature.location.nofuzzy_start, node_start, node_end

					
		new_feature={}
		
		new_feature["location"]=[node_start, node_end]
		

		
		for key in ["primary_name", "gene", "systematic_id", "locus_tag"]:
			if key in feature.qualifiers and feature.qualifiers[key][0]!="":
				new_feature["name"]=feature.qualifiers[key][0]
				break
		if not new_feature.has_key("name"):
			new_feature["name"]="no name"
		
		new_feature["strand"]=feature.location.strand
		
		
		
	
		if node_start==feature.location.nofuzzy_start and node_end==feature.location.nofuzzy_end:
			new_feature["frameshift"]=False
		else:
			new_feature["frameshift"]=True
		
		node_data.comment["annotation"].append(new_feature)
		
	
	treeObject.node(node).set_data(node_data)
	
			
	#sys.exit()
	return treeObject






def apply_annotation_to_branches(treeObject, annotationObject, node=-1, genetic_code_number=1):

	if node==-1:
		node=treeObject.root
		apply_annotation_to_root(treeObject, annotationObject, genetic_code_number=genetic_code_number)
	
	def recurse_root_annotation_across_tree(treeObject, node=-1):
		
		daughters=treeObject.node(node).get_succ()
		
		node_data=treeObject.node(node).get_data()
		
			
		
		
		for daughter in daughters:
			
			
			daughter_data=treeObject.node(daughter).get_data()
			seqlen=len(str(daughter_data.comment["sequence"]))
			revcompseq=Si_SNPs_temp.revcomp(str(daughter_data.comment["sequence"][:]))
			
			if not daughter_data.comment.has_key("annotation"):
				daughter_data.comment["annotation"]=[]
			
			
					
			for feature in node_data.comment["annotation"]:

				if feature["strand"]==1:
					
					daughter_start, daughter_end=Si_SNPs_temp.find_gene_limits(str(daughter_data.comment["sequence"][feature["location"][0]:feature["location"][1]+9999]),0,feature["location"][1]-feature["location"][0], genetic_code_number=genetic_code_number)
					daughter_start=daughter_start+feature["location"][0]
					daughter_end=daughter_end+feature["location"][0]
#					if feature["location"][1]!=daughter_end:
#						print node, daughter, "f", feature["location"][0], feature["location"][1], feature["location"][1]-feature["location"][0], daughter_start, daughter_end, daughter_end-daughter_start
				else:
					
					daughter_start, daughter_end=Si_SNPs_temp.find_gene_limits(revcompseq[seqlen-feature["location"][1]:(seqlen-feature["location"][0])+9999],  0,feature["location"][1]-feature["location"][0], genetic_code_number=genetic_code_number)#need to work out start in revcomp?
					
					#daughter_start=daughter_start+feature["location"][0]
					#daughter_end=daughter_end+feature["location"][0]
					
					daughter_start=feature["location"][1]-daughter_start
					daughter_end=feature["location"][1]-daughter_end
					
					
					daughter_start_new=daughter_end
					daughter_end=daughter_start
					daughter_start=daughter_start_new
					
#					if feature["location"][1]!=daughter_end:
#						print node, daughter, "r", feature["location"][0], feature["location"][1], feature["location"][1]-feature["location"][0], daughter_start, daughter_end, daughter_end-daughter_start

					#print daughter_start, daughter_end, feature["location"][0], feature["location"][1]
					
		
			
				
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
			
			recurse_root_annotation_across_tree(treeObject, node=daughter)
		
		return treeObject
	

	
	recurse_root_annotation_across_tree(treeObject, node=node)
	
	return treeObject







def apply_annotation_to_node(treeObject, annotationObject, node, genetic_code_number=1):
		
	
	node_data=treeObject.node(node).get_data()

		
	seqlen=len(str(node_data.comment["sequence"]))
	revcompseq=Si_SNPs_temp.revcomp(str(node_data.comment["sequence"][:]))
		
	if not node_data.comment.has_key("annotation"):
		node_data.comment["annotation"]=[]

					
	for feature in annotationObject.features:
		
		if feature.location.strand==1:
			
			node_start, possible_node_starts, node_end, possible_node_ends=Si_SNPs_temp.find_gene_limits_new(str(node_data.comment["sequence"][feature.location.nofuzzy_start:feature.location.nofuzzy_end+9999]),0,feature.location.nofuzzy_end-feature.location.nofuzzy_start, genetic_code_number=genetic_code_number)
			
			#print "f", feature.location.nofuzzy_start, feature.location.nofuzzy_end, feature.location.nofuzzy_end-feature.location.nofuzzy_start
			
			node_start=node_start+feature.location.nofuzzy_start
			node_end=node_end+feature.location.nofuzzy_start
			
#			if feature.location.nofuzzy_end!=node_end:
#				print "f", feature.location.nofuzzy_start, feature.location.nofuzzy_end, feature.location.nofuzzy_end-feature.location.nofuzzy_start, node_start, node_end, node_end-node_start

		else:
			
			node_start, possible_node_starts, node_end, possible_node_ends=Si_SNPs_temp.find_gene_limits_new(revcompseq[seqlen-feature.location.nofuzzy_end:(seqlen-feature.location.nofuzzy_start)+9999],  0,feature.location.nofuzzy_end-feature.location.nofuzzy_start, genetic_code_number=genetic_code_number)
			
			#print "r", feature.location.nofuzzy_start, feature.location.nofuzzy_end, feature.location.nofuzzy_end-feature.location.nofuzzy_start
			
			
			node_start=feature.location.nofuzzy_end-node_start
			node_end=feature.location.nofuzzy_end-node_end
			
			node_start_new=node_end
			node_end=node_start
			node_start=node_start_new
			
#			if feature.location.nofuzzy_end!=node_end:
#				print "r", feature.location.nofuzzy_start, feature.location.nofuzzy_end, feature.location.nofuzzy_end-feature.location.nofuzzy_start, node_start, node_end, node_end-node_start

					
		new_feature={}
		
		new_feature["location"]=[node_start, node_end]
		
#		if len(possible_node_ends)>0:
#			print node, feature, feature.location.nofuzzy_start, feature.location.nofuzzy_end
		
		for key in ["primary_name", "gene", "systematic_id", "locus_tag"]:
			if feature.qualifiers.has_key(key):
				new_feature["name"]=feature.qualifiers[key][0]
				break
		if not new_feature.has_key("name"):
			new_feature["name"]="no name"
		
		new_feature["strand"]=feature.location.strand
		
		
		
	
		if node_start==feature.location.nofuzzy_start and node_end==feature.location.nofuzzy_end:
			new_feature["frameshift"]=False
		else:
			new_feature["frameshift"]=True
		
		node_data.comment["annotation"].append(new_feature)
		
	
	treeObject.node(node).set_data(node_data)
	
			
	#sys.exit()
	return treeObject



def apply_annotation_to_branches_new(treeObject, annotationObject, node=-1):
	
	def recurse_annotation_across_tree(treeObject, node=-1):
		
		daughters=treeObject.node(node).get_succ()

		
		for daughter in daughters:
			
			recurse_annotation_across_tree(treeObject, node=daughter)
		apply_annotation_to_node(treeObject, annotationObject, node)
		
#			daughter_data=treeObject.node(daughter).get_data()
#			seqlen=len(str(daughter_data.comment["sequence"]))
#			revcompseq=Si_SNPs_temp.revcomp(str(daughter_data.comment["sequence"][:]))
#			
#			if not daughter_data.comment.has_key("annotation"):
#				daughter_data.comment["annotation"]=[]
#			
#			
#					
#			for feature in node_data.comment["annotation"]:
#
#				if feature["strand"]==1:
#					
#					daughter_start, daughter_end=Si_SNPs_temp.find_gene_limits(str(daughter_data.comment["sequence"][feature["location"][0]:feature["location"][1]+9999]),0,feature["location"][1]-feature["location"][0])
#					daughter_start=daughter_start+feature["location"][0]
#					daughter_end=daughter_end+feature["location"][0]
#					if feature["location"][1]!=daughter_end:
#						print node, daughter, "f", feature["location"][0], feature["location"][1], feature["location"][1]-feature["location"][0], daughter_start, daughter_end, daughter_end-daughter_start
#				else:
#					
#					daughter_start, daughter_end=Si_SNPs_temp.find_gene_limits(revcompseq[seqlen-feature["location"][1]:(seqlen-feature["location"][0])+9999],  0,feature["location"][1]-feature["location"][0])#need to work out start in revcomp?
#					
#					#daughter_start=daughter_start+feature["location"][0]
#					#daughter_end=daughter_end+feature["location"][0]
#					
#					daughter_start=feature["location"][1]-daughter_start
#					daughter_end=feature["location"][1]-daughter_end
#					
#					
#					daughter_start_new=daughter_end
#					daughter_end=daughter_start
#					daughter_start=daughter_start_new
#					
#					if feature["location"][1]!=daughter_end:
#						print node, daughter, "r", feature["location"][0], feature["location"][1], feature["location"][1]-feature["location"][0], daughter_start, daughter_end, daughter_end-daughter_start
#
#					#print daughter_start, daughter_end, feature["location"][0], feature["location"][1]
#					
#		
#			
#				
#				new_feature={}
#				
#				new_feature["location"]=[daughter_start, daughter_end]
#				
#				new_feature["name"]=feature["name"]
#				
#				new_feature["strand"]=feature["strand"]
#				
#				
#				
#				#print node_start, feature.location.nofuzzy_start, node_end, feature.location.nofuzzy_end
#				
#				if daughter_start==feature["location"][0] and daughter_end==feature["location"][1]:
#					new_feature["frameshift"]=False
#				else:
#					new_feature["frameshift"]=True
#				
#				daughter_data.comment["annotation"].append(new_feature)
#				
#			
#			treeObject.node(daughter).set_data(daughter_data)
			
			
		
		return treeObject
	

	
	recurse_annotation_across_tree(treeObject, node=treeObject.root)
	
	return treeObject



















def print_changed_genes(treeObject, node=-1, handle=False, summary_handle=False):

	if node==-1:
		node=treeObject.root
		
	
	def print_node_changed_genes(node=-1):

		
		daughters=treeObject.node(node).get_succ()
		
		node_data=treeObject.node(node).get_data()
		
			
		if node==treeObject.root:
			parent_name="root"
		else:
			parent_name=node
		
		for daughter in daughters:
			
			
			daughter_data=treeObject.node(daughter).get_data()
			
			if "annotation" in daughter_data.comment:
			
				if treeObject.is_internal(daughter):
					daughter_name=daughter
				else:
					daughter_name=daughter_data.taxon
					
				for featurenum, feature in enumerate(node_data.comment["annotation"]):

					
					daughter_start=daughter_data.comment["annotation"][featurenum]["location"][0]
					daughter_end=daughter_data.comment["annotation"][featurenum]["location"][1]
					

				
					if handle and ((feature["strand"]==1 and daughter_end!=feature["location"][1]) or (feature["strand"]==-1 and daughter_start!=feature["location"][0])):
			
					
						if feature["strand"]==1:
							print >> handle, "FT   CDS             "+str(daughter_start+1)+".."+str(daughter_end)
						else:
							print >> handle, "FT   CDS             complement("+str(daughter_start+1)+".."+str(daughter_end)+")"
						
						print >>handle, 'FT                   /primary_name="'+feature['name']+'"'
						print >>handle, 'FT                   /node="'+str(parent_name)+'->'+str(daughter_name)+'"'
						if (feature["strand"]==1 and daughter_end>feature["location"][1]) or (feature["strand"]==-1 and daughter_start<feature["location"][0]):
							print >> handle, 'FT                   /colour=4'
						else:
							print >> handle, 'FT                   /colour=11'
					
						if summary_handle:
					
					
							tostop=[]
							fromstop=[]
							insertions=[]
							deletions=[]
						
							if "SNP_locations" in daughter_data.comment:
								for y in daughter_data.comment['SNP_locations']:
									SNP=daughter_data.comment['SNP_locations'][y]
									if y<=daughter_end and y>daughter_start:
										if SNP.codon_type=="M":
											fromstop.append(str(y))
										elif SNP.codon_type=="O":
											tostop.append(str(y))
										
						
							if "insertion_locations" in daughter_data.comment:
								for insertion in daughter_data.comment['insertion_locations']:
									if insertion[1]<=daughter_end and insertion[0]>daughter_start:
										insertions.append(str(insertion[0])+".."+str(insertion[1]))
						
							if "deletion_locations" in daughter_data.comment:
								for deletion in daughter_data.comment['deletion_locations']:
									if deletion[1]<=daughter_end and deletion[0]>daughter_start:
										deletions.append(str(deletion[0])+".."+str(deletion[1]))
						
					
							if feature["strand"]==1:
						
								if daughter_end>feature["location"][1]:
									changetype="increase"
									change=daughter_end-feature["location"][1]
								else:
									changetype="decrease"
									change=feature["location"][1]-daughter_end
								
							elif feature["strand"]==-1:
								if daughter_start<feature["location"][0]:
									changetype="increase"
									change=feature["location"][0]-daughter_start
								else:
									changetype="decrease"
									change=daughter_start-feature["location"][0]
						
							try: percentchange=(float(change)/(feature["location"][1]-feature["location"][0]))*100
							except ZeroDivisionError:
								percentchange="-"
						
							print >> summary_handle, '\t'.join(map(str,[ feature["name"], feature["strand"], parent_name, daughter_name, feature["location"][1]-feature["location"][0] , daughter_end-daughter_start, changetype, change, percentchange, ', '.join(fromstop), ', '.join(tostop), ', '.join(insertions), ', '.join(deletions)]))
				
				
			
			print_node_changed_genes(node=daughter)
		
		
	

	print >> summary_handle, '\t'.join([ "CDS", "strand", "parent_node", "daughter_node", "parent CDS length", "daughter CDS length", "Increase/Decrease", "length change", "percent change", "STOP->non-STOP", "Non-STOP->STOP", "Insertions", "Deletions"])
	print_node_changed_genes(node=node)
	
	return











	
	
def calculate_branch_dNdS(treeObject, handle, node=-1, genetic_code_number=1):

	if node==-1:
		node=treeObject.root
		print >> handle, '\t'.join(["Branch", 'N', 'S', 'Nd', 'Sd', 'pN', 'pS', 'dN', 'dS', 'varianceN', 'varianceS', 'z', 'dN/dS', '2-tailed Fisher\'s exact test p-value'])
	
	daughters=treeObject.node(node).get_succ()
	
	
	
	for daughter in daughters:
		
		treeObject=calculate_branch_dNdS(treeObject, handle, node=daughter, genetic_code_number=genetic_code_number)
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
		Fisher=0
		
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
		
		Fisher=fisher.pvalue(N,S,Nd,Sd).two_tail
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
			
			node_data.comment["dNdS"]={'N':N, 'S':S, 'dN':dN, 'dS':dS, 'pN':pN, 'pS':pS, 'varS':varianceS, 'varN':varianceN, 'z':z, 'CDSlen':CDSlen, 'Fisher':Fisher}
			
		else:
			#print "Too divergent for JC! Using pN/pS instead."
			dS=pS
			dN=pN
		
		
		
		
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

		node_data.comment["dNdS"]={'N':N, 'S':S, 'dN':dN, 'dS':dS, 'pN':pN, 'pS':pS, 'varS':varianceS, 'varN':varianceN, 'z':z, 'CDSlen':CDSlen, 'dNdS':dnds, 'Fisher':Fisher}
		
		print >> handle, '\t'.join(map(str, [node_node+"->"+daughter_node, N, S, Nd, Sd, pN, pS, dN, dS, varianceN, varianceS, z, dnds, Fisher]))
		print node_node+"->"+daughter_node, dnds, Fisher
		treeObject.node(daughter).set_data(node_data)
		
	return treeObject
			




def print_summary(treeObject, handle, node=-1):

	if node==-1:
		node=treeObject.root
		print >> handle, '\t'.join(["Position","Variant type","Length", "Parent node", "Daughter node", "Parent base", "Daughter base", "Parent amino acid", "Daughter amino acid", "codon type", "CDS name", "Strand", "Position_in_CDS", "Homoplasy", "Removed Homoplasy", "Removed Recombination donor"])

	daughters=treeObject.node(node).get_succ()
	node_data=treeObject.node(node).get_data()
	if node==treeObject.root:
		node_node="root"
	else:
		node_node=str(node)

	for daughter in daughters:

		print_summary(treeObject, handle, daughter)

		daughter_data=treeObject.node(daughter).get_data()
		
		if not treeObject.is_internal(daughter):
			daughter_node=daughter_data.taxon
		else:
			daughter_node=str(daughter)
		
		
		locandtype=[]
		if "SNP_locations" in daughter_data.comment:
			for SNP_location in daughter_data.comment["SNP_locations"]:
				if daughter_data.comment["SNP_locations"][SNP_location].originalSNP:
					locandtype.append([SNP_location, "SNP_locations"])
		if "old_SNP_locations" in daughter_data.comment:
			for SNP_location in daughter_data.comment["old_SNP_locations"]:
				if daughter_data.comment["old_SNP_locations"][SNP_location].originalSNP:
					locandtype.append([SNP_location, "old_SNP_locations"])
			
		locandtype.sort()	
			
		for x in xrange(0,len(locandtype)):
				
				SNP_location=locandtype[x][0]
				SNP_type=locandtype[x][1]
				
				SNP=daughter_data.comment[SNP_type][SNP_location]

				homoplasyline=[]
				oldhomoplasyline=[]

				
				
				if SNP.homoplasy:
					for htype in SNP.homoplasies:
						if htype[0]=="c":
							homoplasyline.append("c with "+str(htype[2]))
						elif htype[0]=="r":
							homoplasyline.append("r from "+str(htype[2]))
						elif htype[0]=="d":
							homoplasyline.append("r in "+str(htype[2]))
				if SNP.oldhomoplasy:
					for htype in SNP.oldhomoplasies:
						if not htype in SNP.homoplasies:
							if htype[0]=="c":
								oldhomoplasyline.append("c with "+str(htype[2])+"")
							elif htype[0]=="r":
								oldhomoplasyline.append("r from "+str(htype[2])+"")
							elif htype[0]=="d":
								oldhomoplasyline.append("r in "+str(htype[2])+"")

				recs=[]
				
				if SNP_type=="old_SNP_locations" and SNP.recombination and len(SNP.recombination):

					if not treeObject.is_internal(SNP.recombination["donor"]):
						donor_node_name=treeObject.node(SNP.recombination["donor"]).data.taxon
					else:
						donor_node_name=str(SNP.recombination["donor"])
					
					recs.append(donor_node_name)
#				elif SNP_type=="old_SNP_locations" and not SNP.recombination:
#					print "Shouldn't get here!"
				
				print >> handle, '\t'.join(map(str,[SNP.position+1,"SNP","1", node_node, daughter_node,SNP.parent_base[0], SNP.daughter_base[0], SNP.parent_aminoacid, SNP.daughter_aminoacid, SNP.codon_type, SNP.CDSname, SNP.strand, SNP.position_in_CDS+1, ', '.join(homoplasyline), ', '.join(oldhomoplasyline), ', '.join(recs)]))


			
		if "insertion_locations" in daughter_data.comment:
			for insertion in daughter_data.comment['insertion_locations']:
				print >> handle, '\t'.join(map(str,[insertion[0]+1,"Insertion",(insertion[1]+1)-insertion[0],node_node, daughter_node,"", "", "", "", "", "", "", "", "", "", ""]))
						
				#insertions.append(str(insertion[0])+".."+str(insertion[1]))
						
		if "deletion_locations" in daughter_data.comment:
			for deletion in daughter_data.comment['deletion_locations']:
				#print deletion
				print >> handle, '\t'.join(map(str,[deletion[0]+1,"Deletion",(deletion[1]+1)-deletion[0],node_node,daughter_node,"", "", "", "", "", "", "", "", "", "", ""]))
				#deletions.append(str(deletion[0])+".."+str(deletion[1]))








def write_tab_output(treeObject, handle, node=-1, colour_snps_by="synonymous"):

	if node==-1:
		node=treeObject.root
	   
	root_data=treeObject.node(treeObject.root).get_data()
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
				if colour_snps_by in ["base", "homoplasy", "homoplasy_bases"]:
					  
					rootbase=root_data.comment["sequence"][daughter_data.comment["SNP_locations"][SNP_location].position].upper()
					node_base=daughter_data.comment["sequence"][daughter_data.comment["SNP_locations"][SNP_location].position].upper()

					downstream_nodes=get_downstream_nodes(treeObject,daughter)
					taxon_list=[]
					for downstream_node in downstream_nodes:
						if treeObject.is_terminal(downstream_node):
							taxon_data=treeObject.node(downstream_node).get_data()
							taxon_base=taxon_data.comment["sequence"][daughter_data.comment["SNP_locations"][SNP_location].position].upper()
							if taxon_base==node_base:
								taxon_list.append(taxon_data.taxon)
					daughter_data.comment["SNP_locations"][SNP_location].write_tab_format(handle, taxon_list ,colourby=colour_snps_by, nodebase=node_base, rootbase=rootbase)
				else:
					daughter_data.comment["SNP_locations"][SNP_location].write_tab_format(handle, treeObject.get_taxa(node_id=daughter) ,colourby=colour_snps_by)
	
		
		if colour_snps_by=="synonymous" and daughter_data.comment.has_key("insertion_locations"):
			for insertion_location in daughter_data.comment["insertion_locations"]:
				
				print >> handle, 'FT   insertion       '+str(insertion_location[0]+1)+".."+str(insertion_location[1]+1)
				print >> handle, 'FT                   /node="'+str(node_node)+'->'+str(daughter_node)+'"'
				print >> handle, 'FT                   /colour=7'
				strain_list=treeObject.get_taxa(node_id=daughter)
				if len(strain_list)>0:
					print >> handle, 'FT                   /taxa="'+', '.join(strain_list)+'"'
					
					
		if colour_snps_by=="synonymous" and daughter_data.comment.has_key("deletion_locations"):
			for deletion_location in daughter_data.comment["deletion_locations"]:
				print >> handle, 'FT   deletion        '+str(deletion_location[0]+1)+".."+str(deletion_location[1]+1)
				print >> handle, 'FT                   /node="'+str(node_node)+'->'+str(daughter_node)+'"'
				print >> handle, 'FT                   /colour=10'
				strain_list=treeObject.get_taxa(node_id=daughter)
				if len(strain_list)>0:
					print >> handle, 'FT                   /taxa="'+', '.join(strain_list)+'"'


		
		if colour_snps_by=="synonymous" and daughter_data.comment.has_key("annotation") and node_data.comment.has_key("annotation"):
		
			for x, feature in enumerate(daughter_data.comment["annotation"]):
				if (feature["strand"]==1 and feature["location"][1]!=node_data.comment["annotation"][x]["location"][1]) or (feature["strand"]==-1 and feature["location"][0]!=node_data.comment["annotation"][x]["location"][0]):

#				if feature["location"][1]!=node_data.comment["annotation"][x]["location"][1]:
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




def reset_homoplasies(treeObject, node=-1):
	
	if node==-1:
		node=treeObject.root
	
	nodedata=treeObject.node(node).get_data()
	
	if "SNP_locations" in nodedata.comment:
		for loc in nodedata.comment["SNP_locations"]:
			if nodedata.comment["SNP_locations"][loc].homoplasy:
				nodedata.comment["SNP_locations"][loc].homoplasy=False
				nodedata.comment["SNP_locations"][loc].homoplasies=[]
		treeObject.node(node).set_data(nodedata)
	
	daughters=treeObject.node(node).get_succ()
	for daughter in daughters:
		treeObject=reset_homoplasies(treeObject, daughter)
	
	return treeObject




def identify_homoplasies(treeObject, locations=[], original=True):
	
	
	def add_SNP_locations_to_list(treeObject, node, SNPlist, locations):
		
		daughters=treeObject.node(node).get_succ()
		
		for daughter in daughters:
		
			data=treeObject.node(daughter).get_data().comment
			
			if data.has_key("SNP_locations"):
			
				if locations==[]:
					SNPlocations=data["SNP_locations"]
				else:
					SNPlocations=[]
					for location in locations:	
						if location in data["SNP_locations"]:
							SNPlocations.append(location)
			
				for location in SNPlocations:
					if not SNPlist.has_key(location):
						SNPlist[location]=[]
					SNPlist[location].append([node, daughter, data["SNP_locations"][location].parent_base, data["SNP_locations"][location].daughter_base])
			
			if treeObject.is_internal(daughter):
				SNPlist=add_SNP_locations_to_list(treeObject, daughter, SNPlist, locations)
	
		
	
		return SNPlist
	
	
	SNPlist={}
	SNPlist=add_SNP_locations_to_list(treeObject, treeObject.root, SNPlist, locations)
	#homoplasies={}

	for key in SNPlist.keys():
		if len(SNPlist[key])>1:
			for x, SNP in enumerate(SNPlist[key]):
				for SNPb in SNPlist[key][x+1:]:
				
					if bases_to_ambiguity[''.join(SNP[3])]==bases_to_ambiguity[''.join(SNPb[3])]:
						#print "Base", key, "convergence", str(SNP[0])+"->"+str(SNP[1]), "and", str(SNPb[0])+"->"+str(SNPb[1]), "to", bases_to_ambiguity[''.join(SNP[3])]
						
						node_data=treeObject.node(SNP[1]).get_data()
						if treeObject.is_terminal(SNPb[1]):
							nodename=treeObject.node(SNPb[1]).data.taxon
						else:
							nodename=str(SNPb[1])
						if original:
							node_data.comment["SNP_locations"][key].oldhomoplasy=True
							node_data.comment["SNP_locations"][key].oldhomoplasies.append(["c", SNPb[1], nodename])
						node_data.comment["SNP_locations"][key].homoplasy=True
						node_data.comment["SNP_locations"][key].homoplasies.append(["c", SNPb[1], nodename])
						treeObject.node(SNP[1]).set_data(node_data)
						node_data=treeObject.node(SNPb[1]).get_data()
						if treeObject.is_terminal(SNP[1]):
							nodename=treeObject.node(SNP[1]).data.taxon
						else:
							nodename=str(SNP[1])
						if original:
							node_data.comment["SNP_locations"][key].oldhomoplasy=True
							node_data.comment["SNP_locations"][key].oldhomoplasies.append(["c", SNP[1], nodename])
						node_data.comment["SNP_locations"][key].homoplasy=True
						node_data.comment["SNP_locations"][key].homoplasies.append(["c", SNP[1], nodename])
						treeObject.node(SNPb[1]).set_data(node_data)
						
						#homoplasies[key]=[[SNP[1]]+get_downstream_nodes(treeObject, SNP[1]),[SNPb[1]]+get_downstream_nodes(treeObject, SNPb[1])]
						
						
						
					elif bases_to_ambiguity[''.join(SNP[3])]==bases_to_ambiguity[''.join(SNPb[2])] and  treeObject.common_ancestor(SNP[1], SNPb[1])==SNPb[1]:
#						if str(SNP[1])=="12" or str(SNPb[1])=="12":
#							print "Base", key, "reversal", str(SNP[0])+"->"+str(SNP[1]), "to", bases_to_ambiguity[''.join(SNP[3])], "and", str(SNPb[0])+"->"+str(SNPb[1]), "to", bases_to_ambiguity[''.join(SNPb[3])]


						
						node_data=treeObject.node(SNP[1]).get_data()
						if treeObject.is_terminal(SNPb[1]):
							nodename=treeObject.node(SNPb[1]).data.taxon
						else:
							nodename=str(SNPb[1])
						if original:
							node_data.comment["SNP_locations"][key].oldhomoplasy=True
							node_data.comment["SNP_locations"][key].oldhomoplasies.append(["r", SNPb[1], nodename])
						
						node_data.comment["SNP_locations"][key].homoplasy=True
						node_data.comment["SNP_locations"][key].homoplasies.append(["r", SNPb[1], nodename])
						treeObject.node(SNP[1]).set_data(node_data)
						node_data=treeObject.node(SNPb[1]).get_data()
						if treeObject.is_terminal(SNP[1]):
							nodename=treeObject.node(SNP[1]).data.taxon
						else:
							nodename=str(SNb[1])
						if original:
							node_data.comment["SNP_locations"][key].oldhomoplasy=True
							node_data.comment["SNP_locations"][key].oldhomoplasies.append(["d", SNP[1], nodename])
						node_data.comment["SNP_locations"][key].homoplasy=True
						node_data.comment["SNP_locations"][key].homoplasies.append(["d", SNP[1], nodename])
						treeObject.node(SNPb[1]).set_data(node_data)
						
					elif bases_to_ambiguity[''.join(SNP[2])]==bases_to_ambiguity[''.join(SNPb[3])] and  treeObject.common_ancestor(SNP[1], SNPb[1])==SNP[1]:
#						if str(SNP[1])=="12" or str(SNPb[1])=="12":
#							print "Base", key, "reversal", str(SNP[0])+"->"+str(SNP[1]), "to", bases_to_ambiguity[''.join(SNP[3])], "and", str(SNPb[0])+"->"+str(SNPb[1]), "to", bases_to_ambiguity[''.join(SNPb[3])]

						node_data=treeObject.node(SNP[1]).get_data()
						if treeObject.is_terminal(SNPb[1]):
							nodename=treeObject.node(SNPb[1]).data.taxon
						else:
							nodename=str(SNPb[1])
						if original:
							node_data.comment["SNP_locations"][key].oldhomoplasy=True
							node_data.comment["SNP_locations"][key].oldhomoplasies.append(["d", SNPb[1], nodename])
						node_data.comment["SNP_locations"][key].homoplasy=True
						node_data.comment["SNP_locations"][key].homoplasies.append(["d", SNPb[1], nodename])
						treeObject.node(SNP[1]).set_data(node_data)
						node_data=treeObject.node(SNPb[1]).get_data()
						if treeObject.is_terminal(SNP[1]):
							nodename=treeObject.node(SNP[1]).data.taxon
						else:
							nodename=str(SNP[1])
						if original:
							node_data.comment["SNP_locations"][key].oldhomoplasy=True
							node_data.comment["SNP_locations"][key].oldhomoplasies.append(["r", SNP[1], nodename])
						node_data.comment["SNP_locations"][key].homoplasy=True
						node_data.comment["SNP_locations"][key].homoplasies.append(["r", SNP[1], nodename])
						treeObject.node(SNPb[1]).set_data(node_data)
						
			#print homoplasies			
				
	return treeObject
			







def get_homoplasic_sites(treeObject):
	
	rootnode=treeObject.root
	homoplasic_sites=set()
	
	def get_node_homoplasies(treeObject, node, homoplasic_sites):
		
		daughters=treeObject.node(node).get_succ()
		
		for daughter in daughters:
			homoplasic_sites=get_node_homoplasies(treeObject, daughter, homoplasic_sites)
		
		node_data=treeObject.node(node).get_data()
		
		if node_data.comment.has_key("SNP_locations"):
			locations=node_data.comment["SNP_locations"]
		else:
			locations=[]
		
		for location in locations:
			if node_data.comment["SNP_locations"][location].homoplasy and location not in homoplasic_sites:
				#print location, node_data.comment["SNP_locations"][location].homoplasy
				homoplasic_sites.add(location)
	
		return homoplasic_sites
	
	homoplasic_sites=get_node_homoplasies(treeObject, rootnode, homoplasic_sites)

	return list(homoplasic_sites)








def genome_diagram_for_tree(treeObject,filename,ladderize=None, printtype="SNPs", colourby="SNPtype", referenceObject=None, fragments=1, locations=[]):
   
   
	def SNPs_to_features(node, features=[], depth=0, colourby="SNPtype", taxonsets=[], locations=[]):
		
		node_data=treeObject.node(node).get_data()
		if colourby=="base":
			root_data=treeObject.node(treeObject.root).get_data()
			base_colour={"A": "2", "C": "3", "G": "10", "T": "4", "-":"1", "N":"13"}
		elif colourby=="taxa":
			colour_list=['1','2','3','4','5','6','8','9','10','11','12','13','14','15']
		
	
		if locations==[] or colourby in ["SNPtype", "homoplasy"]:
			if node_data.comment.has_key("SNP_locations"):
				locations=node_data.comment["SNP_locations"]
		
		for SNP in locations:
			
			if colourby=="SNPtype":
				feature=SeqFeature(FeatureLocation(SNP,SNP), strand=1, type="SNP", qualifiers={"colour":node_data.comment["SNP_locations"][SNP].colour})
				features.append([feature,depth])
			
			elif colourby=="homoplasy":
				if node_data.comment["SNP_locations"][SNP].homoplasy:
					feature=SeqFeature(FeatureLocation(SNP,SNP), strand=1, type="SNP", qualifiers={"colour":"2"})
					features.append([feature,depth])
#					else:
#						feature=SeqFeature(FeatureLocation(node_data.comment["SNP_locations"][SNP].position,node_data.comment["SNP_locations"][SNP].position), strand=1, type="SNP", qualifiers={"colour":"1"})


			elif colourby=="base":
				#if node_data.comment["SNP_locations"][SNP].homoplasy:
				rootbase=root_data.comment["sequence"][SNP].upper()
				node_base=node_data.comment["sequence"][SNP].upper()
				#print node_base, rootbase
				
				if rootbase!=node_base and rootbase in ["A", "C", "G", "T", "-", "N"] and node_base in ["A", "C", "G", "T", "-", "N"]:#node_base in ["A", "C", "G", "T"]:#
					feature=SeqFeature(FeatureLocation(SNP,SNP), strand=0, type="SNP", qualifiers={"colour":base_colour[node_base]})
					features.append([feature,depth])
			elif colourby=="taxa":
				#if node_data.comment["SNP_locations"][SNP].homoplasy:
				if treeObject.node(node).get_data().comment["sequence"][SNP]==treeObject.node(treeObject.root).get_data().comment["sequence"][SNP]:
					continue
				taxaset=set()
				for taxon in treeObject.get_terminals():
					if treeObject.node(taxon).get_data().comment["sequence"][SNP]==treeObject.node(node).get_data().comment["sequence"][SNP]:
						taxaset.add(taxon)
						
				snpcolour=-1
				x=0
				for x, tset in enumerate(taxonsets):
					
					if len(tset.symmetric_difference(taxaset))==0:
						snpcolour=x
						break
			
				if snpcolour==-1:
					snpcolour=x+1
					taxonsets.append(taxaset)

				while snpcolour>=len(colour_list):
					snpcolour=snpcolour-len(colour_list)

			#else:
				#snpcolour=0
			
			
				feature=SeqFeature(FeatureLocation(SNP,SNP), strand=0, type="SNP", qualifiers={"colour":colour_list[snpcolour]})
				features.append([feature,depth])
			
				
		
		return features, taxonsets
		
		
   
	def add_SNP_tracks(node,track_number=1,ladderize=None, previous_features=[], depth=0, colourby="SNPtype", taxonsets=[], locations=[]): 
		"""Convert a node tree to a newick tree recursively."""
   			
   		
		if treeObject.node(node).get_succ():
   			succnodes=ladderize_nodes(treeObject.node(node).get_succ(),ladderize=ladderize)
   			succnodes.reverse()
   			if locations==[] or colourby in ["SNPtype", "homoplasy"]:
				features, taxonsets=SNPs_to_features(node, previous_features[:], depth, colourby=colourby, taxonsets=taxonsets, locations=locations)
			else:
				features, taxonsets=SNPs_to_features(node, [], depth, colourby=colourby, taxonsets=taxonsets, locations=locations)
	   		for succnode in succnodes:
	   			track_number, taxonsets=add_SNP_tracks(succnode,track_number,ladderize=ladderize,previous_features=features, depth=depth+1, colourby=colourby, taxonsets=taxonsets, locations=locations)
   		else:
			gd_track_for_features = gd_diagram.new_track(track_number, scale=0, height=1)
			gd_feature_set = gd_track_for_features.new_set()
	
	
			if locations==[] or colourby in ["SNPtype", "homoplasy"]:
				features, taxonsets=SNPs_to_features(node, previous_features[:], depth, colourby=colourby, taxonsets=taxonsets, locations=locations)
			else:
				features, taxonsets=SNPs_to_features(node, [], depth, colourby=colourby, taxonsets=taxonsets, locations=locations)
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
#				if colourby=="base":
#					print feature, color
				gd_feature_set.add_feature(feature[0], color=color, label=0, arrowhead_length=0.25)

			track_number+=1
   
   		return track_number, taxonsets
   		
   	
   	
   	
   	
   	def changed_genes_to_features(node, features=[], depth=0):
		
		node_data=treeObject.node(node).get_data()
		
		
		
		if node_data.comment.has_key("annotation") and node!=treeObject.root:
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
   			
   		
		if treeObject.node(node).get_succ():
   			succnodes=ladderize_nodes(treeObject.node(node).get_succ(),ladderize=ladderize)
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
			if feature.type not in  ["CDS", "feature"]:#,"tRNA","repeat_unit"] :
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
	
#	succnodes=ladderize_nodes(treeObject.node(treeObject.root).succ)
#	succnodes.reverse()
#	for node in succnodes:
#		if printtype in ["SNPs", "taxa", "base"]:
#			track_number, taxonsets=add_SNP_tracks(node, track_number=track_number,ladderize=ladderize, colourby=colourby, locations=locations)
#			print taxonsets, len(taxonsets)
#		elif printtype=="changed_genes":
#			track_number=add_changed_gene_tracks(node, track_number=track_number,ladderize=ladderize)
#	for tset in taxonsets:		
#		print tset
	if printtype in ["SNPs", "taxa", "base"]:
		track_number, taxonsets=add_SNP_tracks(treeObject.root, track_number=track_number,ladderize=ladderize, colourby=colourby, locations=locations)
	elif printtype=="changed_genes":
		track_number=add_changed_gene_tracks(treeObject.root, track_number=track_number,ladderize=ladderize)
	
	gd_diagram.draw(format="linear", orientation="landscape", pagesize='A4', fragments=fragments)#, start=0, end=len(record))
	
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
		
			daughters=tree.node(node).succ
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
		
			daughters=tree.node(node).succ
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
			daughters=tree.node(node).succ
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
		
		
		for col in xrange(startcol, thiscol): 
			char_matrix[thisrow][col] = '-' 
			
		if tree.is_internal(node):
			daughters=tree.node(node).succ
			# Draw a vertical line 
			toprow = row_positions[daughters[0]] 
			botrow = row_positions[daughters[-1]] 
			for row in xrange(toprow+1, botrow): 
				char_matrix[row][thiscol] = '|'
				if row==thisrow and show_nodes:
					for x in xrange(len(str(node))):
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
			for x in xrange(len(taxon)):
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
		
		if len(bases1)>0 and len(bases2)>0 and len(bases1.intersection(bases2))==0:
			numdiffs+=1
	
	#print numdiffs
	return numdiffs


	
	



#
#def get_max_depth(treeObject, ladderize=None):
#	node=treeObject.root
#
#	def recurse(node, depth, right, left, max_depth, max_right, max_left, ladderize=None):
#	
#		if depth>max_depth:
#			max_depth=depth
#		if left>max_left:
#			max_left=left
#		if right>max_right:
#			max_right=right
#		
#		if treeObject.node(node).get_succ():
#   			daughtersnew=ladderize_nodes(treeObject.node(node).get_succ(),ladderize=ladderize)
#   			daughtersnew.reverse()
#   			for daughter in daughtersnew:
#   				
#	   			if daughter==daughtersnew[0]:
#	   				newdepth=depth+1
#	   				newleft=left+1
#					max_depth, max_right, max_left=recurse(daughter, newdepth, right, newleft, max_depth, max_right, max_left, ladderize=ladderize)
#				else:
#					newdepth=depth+1
#					newright=right+1
#					max_depth, max_right, max_left=recurse(daughter, newdepth, newright, left, max_depth, max_right, max_left, ladderize=ladderize)
#				
#		
#		return max_depth, max_right, max_left
#	
#	
#	max_depth, max_right, max_left=recurse(node, 0, 0, 0, 0, 0, 0, ladderize=ladderize)
#	
#	return max_depth, max_right, max_left





#def get_max_distances(treeObject, ladderize=None):
#	
#	def find_max_distance_to_terminal(node1):
#		
#		max_distance=0
#		
#		for node2 in treeObject.get_terminals():
#			distance=treeObject.distance(node1, node2)
#			if distance>max_distance:
#				max_distance=distance
#		
#		return max_distance
#		
#	
#	def find_top(ladderize=None):
#		
#		
#		top=treeObject.root
#		
#		while treeObject.node(top).succ:
#			daughtersnew=ladderize_nodes(treeObject.node(top).succ,ladderize=ladderize)
#   			daughtersnew.reverse()
#   			top=daughtersnew[0]
#   	
#   		return top
#   	
#   	def find_bottom(ladderize=None):
#		
#		
#		bottom=treeObject.root
#		
#		while treeObject.node(bottom).succ:
#			daughtersnew=ladderize_nodes(treeObject.node(bottom).succ,ladderize=ladderize)
#   			daughtersnew.reverse()
#   			bottom=daughtersnew[1]
#   	
#   		return bottom
#   	
#   	top=find_top(ladderize=ladderize)
#	max_top_distance=find_max_distance_to_terminal(top)
#	bottom=find_bottom(ladderize=ladderize)
#	max_bottom_distance=find_max_distance_to_terminal(bottom)
#	max_root_distance=find_max_distance_to_terminal(treeObject.root)
#	
#	
#	return top, bottom, max_top_distance, max_bottom_distance, max_root_distance


#def get_max_distances(treeObject, ladderize=None):
#	
#	def find_max_distance_to_terminal(node1):
#		
#		max_distance=0
#		
#		for node2 in treeObject.get_terminals():
#			distance=treeObject.distance(node1, node2)
#			if distance>max_distance:
#				max_distance=distance
#		
#		return max_distance		
#		
#	
#	def find_top(ladderize=None):
#		
#		
#		top=treeObject.root
#		
#		while treeObject.node(top).succ:
#			daughtersnew=ladderize_nodes(treeObject.node(top).succ,ladderize=ladderize)
#   			daughtersnew.reverse()
#   			top=daughtersnew[0]
#   	
#   		return top
#   	
#   	def find_bottom(ladderize=None):
#		
#		
#		bottom=treeObject.root
#		
#		while treeObject.node(bottom).succ:
#			daughtersnew=ladderize_nodes(treeObject.node(bottom).succ,ladderize=ladderize)
#   			daughtersnew.reverse()
#   			bottom=daughtersnew[1]
#   	
#   		return bottom
#   	
#   	top=find_top(ladderize=ladderize)
#	max_top_distance=find_max_distance_to_terminal(top)
#	bottom=find_bottom(ladderize=ladderize)
#	max_bottom_distance=find_max_distance_to_terminal(bottom)
#	max_root_distance=find_max_distance_to_terminal(treeObject.root)
#	
#	
#	return top, bottom, max_top_distance, max_bottom_distance, max_root_distance



def get_max_distances(treeObject, ladderize=None):
	
	
	max_distance=0
	max1=0
	max2=0
	for node1 in treeObject.get_terminals():
		for node2 in treeObject.get_terminals():
			distance=treeObject.distance(node1, node2)
			if distance>max_distance:
				max_distance=distance
				max1=node1
				max2=node2
	
	max_distance=0
	for node3 in treeObject.get_terminals():
		distance1=treeObject.distance(max1, node3)
		distance2=treeObject.distance(max2, node3)
		distance=distance1+distance2
		if distance>max_distance:
			max_distance=distance
			max3=node3
	
#	CA=treeObject.common_ancestor(max1,max2)
#	if CA==0:
#		CA=treeObject.common_ancestor(max1,max3)
#	if CA==0:
#		CA=treeObject.common_ancestor(max2,max3)
#	
#   	
#   	max_distance1=treeObject.distance(CA, max1)
#   	max_distance2=treeObject.distance(CA, max2)
#   	max_distance3=treeObject.distance(CA, max3)
   	distance12=treeObject.distance(max1, max2)
   	distance13=treeObject.distance(max1, max3)
   	distance23=treeObject.distance(max2, max3)
	
	#print max1, max2, max3, distance12, distance13, distance23
	
	return max1, max2, max3, distance12, distance13, distance23
			
			




def find_clusters_on_all_branches(treeObject, prefix="MW", homoplasy=False, locations=[]):
	

	handle=open(prefix+"_recombination_key.tre","w")
	
	
	print >> handle, tree_to_figtree_string(treeObject, False, False, False, False)
	
	handle.close()
	
	missing_set=set(["N", "?"])
	missing_and_gaps_set=set(["N", "?", "-"])
	
	def count_sequence_diffs(node1, node2, start, end):

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
			
			if len(bases1)>0 and len(bases2)>0 and len(bases1.intersection(bases2))==0:
				numdiffs+=1
		#print numdiffs
		return numdiffs
	def count_diffs_for_seqsb(seq1, seq2):
		numdiffs=0
		mutations=set()
		for x, ambiguity1 in enumerate(seq1):
	
			bases1=set(ambiguity_to_bases[ambiguity1])
			bases2=set(ambiguity_to_bases[seq2[x]])
	
			if "-" in bases1:
				bases1.remove("-")
			elif "-" in bases2:
				bases2.remove("-")
			
			#print bases1, bases2, bases1.union(bases2)
			
			if len(bases1)>0 and len(bases2)>0 and len(bases1.intersection(bases2))==0:
				numdiffs+=1
				mutations.add(x)
		#print numdiffs
		return mutations
		return numdiffs
	
	def get_consensus_for_seqs(seq1, seq2):

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
	
	
	def get_consensus(node1, node2, start, end):

		seq1=treeObject.node(node1).get_data().comment["sequence"][start:end]
	
		seq2=treeObject.node(node2).get_data().comment["sequence"][start:end]
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
	
	
	
	def node_donor_homoplasies(node, locations=[], start=0, end=0, new=True):
	
		
		daughters=treeObject.node(node).get_succ()
		for daughter in daughters:
			node_donor_homoplasies(daughter, locations)
		
		potential_donors=get_nodes_not_downstream(treeObject, node)

		
		node_seq=treeObject.node(node).get_data().comment["sequence"]
		node_seq_len=len(node_seq)
		if end==0:
			end=node_seq_len
		
		if locations==[]:
				locations=xrange(start,end)
		
		if treeObject.node(node).get_data().comment.has_key("SNP_locations"):
			SNPlocations=treeObject.node(node).get_data().comment["SNP_locations"]
		else:	
			return
		
		#if len(SNPlocations)<minsnps:
			#return
		
		if new or not node in homoplasyposns:
			homoplasyposns[node]={}
			seqposntonogapposn[node]={}
			nogapposntoseqposn[node]={}
			alignmentlengths[node]={}
		else:
			for x in xrange(start,end):
				for donor in potential_donors:
					if donor in homoplasyposns[node] and x in homoplasyposns[node][donor]:
						del homoplasyposns[node][donor][x]
						del seqposntonogapposn[node][donor][x]
						del nogapposntoseqposn[node][donor][x]
						del alignmentlengths[node][donor][x]
		
		
		for donor in potential_donors:
			if new or not donor in homoplasyposns[node]:
				homoplasyposns[node][donor]=[]
				seqposntonogapposn[node][donor]={}
				nogapposntoseqposn[node][donor]={}
			donor_seq=treeObject.node(donor).get_data().comment["sequence"]
			if donor==treeObject.root:
				donor_parent_seq=donor_seq
			else:
				donor_parent_seq=treeObject.node(treeObject.node(donor).get_prev()).get_data().comment["sequence"]

			gapcount=0

				

			
			#speedup appending by avioding dots
			append=homoplasyposns[node][donor].append
			
			for x in locations:
				
				if x<start or x>end:
					continue
				
#				if node==53 and node_seq[x]!="-" and treeObject.node(treeObject.node(node).get_prev()).get_data().comment["sequence"][x]!="-" and x>256 and x<1032 and ((node_seq[x]!=treeObject.node(treeObject.node(node).get_prev()).get_data().comment["sequence"][x] and not x in SNPlocations) or (node_seq[x]==treeObject.node(treeObject.node(node).get_prev()).get_data().comment["sequence"][x] and x in SNPlocations)):
##				if node==53 and x==520:
#					print x, donor, node, donor_seq[x], donor_parent_seq[x], node_seq[x], treeObject.node(treeObject.node(node).get_prev()).get_data().comment["sequence"][x], x in SNPlocations
				
				if node_seq[x] in missing_and_gaps_set or (donor_seq[x] in missing_set and donor_parent_seq[x] in missing_set):
					gapcount+=1
					
					
				#elif x in SNPlocations and SNPlocations[x].homoplasy and (node_seq[x]==donor_seq[x] or node_seq[x]==donor_parent_seq[x]):
				elif x in SNPlocations and node_seq[x]!=treeObject.node(treeObject.node(node).get_prev()).get_data().comment["sequence"][x] and (node_seq[x]==donor_seq[x] or node_seq[x]==donor_parent_seq[x]):
				
				
					append(x-gapcount)
					seqposntonogapposn[node][donor][x]=x-gapcount
					nogapposntoseqposn[node][donor][x-gapcount]=x
					
			alignmentlengths[node][donor]=node_seq_len-gapcount

#			if len(homoplasyposns[node][donor])<minsnps:
#				del homoplasyposns[node][donor]
#				del seqposntonogapposn[node][donor]
#				del nogapposntoseqposn[node][donor]
#				del alignmentlengths[node][donor]
#		
#		if len(homoplasyposns[node])==0:
#			del homoplasyposns[node]
#			del seqposntonogapposn[node]
#			del nogapposntoseqposn[node]
#			del alignmentlengths[node]
		
		
		
		
	def print_blocks_to_tab(treeObject, blocks, handle):
		
		for block in blocks:
			node=block["recipient"]
			
			colour=str(int(treeObject.node(block["common_ancestor"]).data.comment["colour"][0]))+" "+str(int(treeObject.node(block["common_ancestor"]).data.comment["colour"][1]))+" "+str(int(treeObject.node(block["common_ancestor"]).data.comment["colour"][2]))
			
			donorlist=[]
			for donor in block["donors"]:
				if treeObject.is_terminal(donor):
					donorname=treeObject.node(donor).data.taxon
				else:
					donorname=donor
				donorlist.append(str(donorname))
			
			downstreamnamelist=treeObject.get_taxa(node)
			
			if treeObject.is_terminal(node):
				nodename=treeObject.node(node).data.taxon
			else:
				nodename=node
			
			if treeObject.is_terminal(block["common_ancestor"]):
				ancestor=treeObject.node(block["common_ancestor"]).data.taxon
			else:
				ancestor=block["common_ancestor"]
			
			print >> handle, "FT   misc_feature    "+str(block["start"]+1)+".."+str(block["end"]+1)
			print >> handle, "FT                   /node="+str(nodename)
			print >> handle, 'FT                   /possible_donors="'+', '.join(donorlist)+'"'
			print >> handle, "FT                   /donor_common_ancestor="+str(ancestor)
			print >> handle, 'FT                   /rank='+str(block["removal_round"])
			print >> handle, 'FT                   /neg_log_likelihood='+str(block["ll"])
			print >> handle, 'FT                   /recombined_snps='+str(block["snpcount"])	
			print >> handle, 'FT                   /pvalue='+str(block["p"])
			print >> handle, 'FT                   /taxa="'+', '.join(downstreamnamelist)+'"'
			print >> handle, 'FT                   /colour='+colour
	
	
	
	def get_block_likelihood(n, c, N, C):
	
		n=float(n)
		c=float(c)
		N=float(N)
		C=float(C)
	#	for x in snpposns:
	#		if x>=start and x<=end:
	#			c+=1
	
		#print c, n, C, N
		part1=math.log((c/n),10)*c
		
		try:
			part2=math.log((((n-c)/n)),10)*(n-c)
		except (ValueError, ZeroDivisionError):
			part2=0
		try:
			part3=math.log((((C-c)/(N-n))),10)*(C-c)
		except (ValueError, ZeroDivisionError):
			part3=0
		try:
			part4=math.log(((((N-n)-(C-c))/(N-n))),10)*((N-n)-(C-c))
		except (ValueError, ZeroDivisionError):
			part4=0
		
		likelihood=(part1+part2+part3+part4)*-1
		
		#print start, end, c, n, C, N, likelihood
		#print likelihood
		
		return likelihood
	
	
	
	
	
	def reduce_factorial(l,i):
		
		
		if l==i and l<1000:
			return math.log(math.factorial(l),10)
		
		
		factorial=math.log(1.0,10)
		
		
		for x in xrange(int(l)-(int(i)-1),int(l)+1):
			#print x, factorial
			factorial=factorial+math.log(x,10)
		
		
		return factorial
		
		
	
#	def refine_blocks(recipient, donor, snpposns, tempblocks, cutoff, lennogaps, minsnps, N, C):
#	
#			
#		newblocks=[]
#		
#		for block in tempblocks:
#			
#			
#			#trim to first and last SNPs in block
#			start=int(block[1])
#			end=int(block[2])
#			#print "a", start, end
#			x=0
#			
#			while x<(len(snpposns)-1) and snpposns[x]<start:
#				x+=1
#				
#			start=snpposns[x]
#			snpposnstart=x			
#			
#			
#			x=len(snpposns)-1
#			
#			while x>0 and snpposns[x]>end:
#				x-=1
#				
#			end=snpposns[x]
#			snpposnend=x
#			
#
##			for x in xrange(0,len(snpposns)):
##				if snpposns[x]==start:
##					snpposnstart=x
##				if snpposns[x]==end:
##					snpposnend=x
##					break
#			
#	
#			if start>end:
#				continue
#			
#			old_like=get_block_likelihood((end+1)-start, (snpposnend+1)-snpposnstart, N, C)
#			#print "b", start, end, snpposnstart, snpposnend, old_like
#			#print "initial likelihood=", old_like, "start=", start, "end=", end, snpposnstart, snpposnend
#			
#			newstart=start
#			oldstart=start
#			newsnpposnstart=snpposnstart
#			
#			if snpposnstart+1<(snpposnend+1)-minsnps:
#				for y, snpposnloc in enumerate(snpposns[snpposnstart+1: (snpposnend+1)-minsnps]):
#					x=y+1
#					new_like=get_block_likelihood((end+1)-snpposnloc, (snpposnend+1)-(snpposnstart+x), N, C)
#					#print snpposnstart, newsnpposnstart, x, new_like	
#					if new_like>old_like:
#						break
#					else:
#						newstart=snpposnloc
#						newsnpposnstart=snpposnstart+x
#						old_like=new_like
#			
#				start=newstart
#				snpposnstart=newsnpposnstart
#				#print "c",x, start, end, snpposnstart, snpposnend, old_like
#			
#			if start==oldstart and snpposnstart>0:
#				newstart=start
#				
#				for y, snpposnloc  in enumerate( snpposns[snpposnstart-1:0:-1]):
#					x=y+1
#					new_like=get_block_likelihood((end+1)-snpposnloc, (snpposnend+1)-(snpposnstart-x), N, C)
#						
#					if new_like>old_like:
#						break
#					else:
#						newstart=snpposnloc
#						newsnpposnstart=snpposnstart-x
#						old_like=new_like
#			
#				start=newstart
#				snpposnstart=newsnpposnstart
#				#print "d", start, end, snpposnstart, snpposnend, old_like
#			
#			newsnpposnend=snpposnend
#			newend=end
#			oldend=end
#			if snpposnend>=minsnps and (snpposnend-1)>(snpposnstart+minsnps):
#				for y, snpposnloc in enumerate( snpposns[ snpposnend-1: snpposnstart+minsnps: -1]):
#					x=y+1
#					new_like=get_block_likelihood((snpposnloc+1)-start, (snpposnend+1-x)-snpposnstart, N, C)
#					
#					if new_like>old_like:
#						break
#					else:
#						newend=snpposnloc
#						newsnpposnend=snpposnend-x
#						old_like=new_like
#			
#				end=newend
#				snpposnend=newsnpposnend
#				#print "e", start, end, snpposnstart, snpposnend, old_like
#			
#			if end==oldend and snpposnend+1<len(snpposns):
#				newend=end
#				for y, snpposnloc in enumerate(snpposns[ snpposnend+1:]):
#					x=y+1
#					new_like=get_block_likelihood((snpposnloc+1)-start, (snpposnend+1+x)-snpposnstart, N, C)
#					
#					if new_like>old_like:
#						break
#					else:
#						newend=snpposnloc
#						newsnpposnend=snpposnend+x
#						old_like=new_like
#			
#				end=newend
#				snpposnend=newsnpposnend
#				#print "f", start, end, snpposnstart, snpposnend, old_like
#				
#			likelihood=old_like
#			#print "g", start, end, snpposnstart, snpposnend, old_like
#			#print "new likelihood=", old_like, "start=", start, "end=", end
#			
#			snpcount=(snpposnend+1)-snpposnstart
#			if snpcount>=minsnps:
#				newblocks.append([likelihood, start, end, snpposnstart, snpposnend])
#	
#			
#		return newblocks
	
	
		
	
	
	
	
	
	def get_mutations_in_region(node, start, end):
		node_data=treeObject.node(node).get_data().comment
		count=0
		
		SNPlocations=set()
		append=SNPlocations.add
		if "SNP_locations" in node_data:
			for SNP in node_data["SNP_locations"]:
				if SNP>=start and SNP<=end:
					append(SNP)
		return SNPlocations
	
	def get_mutations_between_nodes_oldb(node1, node2, start, end):
	
		intersection=treeObject.common_ancestor(node1, node2)
		
		#print intersection
		
		locations=set()
		
		node=node1
		
		while node!=intersection:
			locations.update(get_mutations_in_region(node, start, end))
			node=treeObject.node(node).get_prev()
			#print locations
		
		
		node=node2
		
		while node!=intersection:
			locations.update(get_mutations_in_region(node, start, end))
			node=treeObject.node(node).get_prev()
			#print locations
	
		node1_seq=treeObject.node(node1).get_data().comment["sequence"]
		if node1!=0:
			node3_seq=treeObject.node(treeObject.node(node1).get_prev()).get_data().comment["sequence"]
		node2_seq=treeObject.node(node2).get_data().comment["sequence"]
		if node2!=0:
			node4_seq=treeObject.node(treeObject.node(node2).get_prev()).get_data().comment["sequence"]
		
		
		numdiffs=0
		
		#print locations
		
		#mutations=set()
		
		for location in locations:
			
			
			if node1!=0:
				bases1=set(ambiguity_to_bases[node1_seq[location]])
				bases3=set(ambiguity_to_bases[node3_seq[location]])
				union1=bases1.union(bases3)
#				if node1 in [2, 83] and node2 in [2, 83]:
#					print bases1, bases3, union1
			else:
				union1=set(ambiguity_to_bases[node1_seq[location]])
			
			if node2!=0:
				bases2=set(ambiguity_to_bases[node2_seq[location]])
				bases4=set(ambiguity_to_bases[node4_seq[location]])
				union2=bases2.union(bases4)
#				if node1 in [2, 83] and node2 in [2, 83]:
#					print bases2, bases4, union2
			else:
				union2=set(ambiguity_to_bases[node2_seq[location]])
			
#			union_list=list(union)
#			union_list.sort()
#			
#			consensus_base=set(bases_to_ambiguity[''.join(union_list)]
	
			if "-" in union1:
				union1.remove("-")
			if "-" in union2:
				union2.remove("-")
			
			#print bases1, bases2, bases1.union(bases2)
			
			if len(union1)>0 and len(union2)>0 and len(union1.intersection(union2))==0:
				numdiffs+=1
				#mutations.add(location)
		#return mutations
		return numdiffs
	
	
	
	
	
	def get_mutations_between_nodes_old(node1, node2, start, end, final=False):
	
		numdiffs=0
		if node1==treeObject.root:
			recipient_parent_seq=treeObject.node(node1).get_data().comment["sequence"]
		else:
			recipient_parent_seq=treeObject.node(treeObject.node(node1).get_prev()).get_data().comment["sequence"]
		recipient_seq=treeObject.node(node1).get_data().comment["sequence"]
		if node2==treeObject.root:
			donor_parent_seq=treeObject.node(node2).get_data().comment["sequence"]
		else:
			donor_parent_seq=treeObject.node(treeObject.node(node2).get_prev()).get_data().comment["sequence"]
		donor_daughter_seq=treeObject.node(node2).get_data().comment["sequence"]
		
		for x in xrange(start,end+1):
			if not recipient_seq[x]==donor_parent_seq[x] and not recipient_seq[x]==donor_daughter_seq[x]:
				numdiffs+=1
				if final:
					print recipient_parent_seq[x], recipient_seq[x], donor_parent_seq[x], donor_daughter_seq[x] 
		

		return numdiffs
	
	
	
	def get_mutations_between_nodes(node1, node2, start, end, snpposns, getposns=False, getimprovement=False):
		
		if getposns:
			mutposns=[]
		
		numdiffs=0
		improvement=0
		if node1==treeObject.root:
			recipient_parent_seq=treeObject.node(node1).get_data().comment["sequence"]
		else:
			recipient_parent_seq=treeObject.node(treeObject.node(node1).get_prev()).get_data().comment["sequence"]
		recipient_seq=treeObject.node(node1).get_data().comment["sequence"]
		if node2==treeObject.root:
			donor_parent_seq=treeObject.node(node2).get_data().comment["sequence"]
		else:
			donor_parent_seq=treeObject.node(treeObject.node(node2).get_prev()).get_data().comment["sequence"]
		donor_daughter_seq=treeObject.node(node2).get_data().comment["sequence"]
		
		for x in snpposns:
			if x<=end and x>=start:
				if not recipient_seq[x]==donor_parent_seq[x] and not recipient_seq[x]==donor_daughter_seq[x] and recipient_seq[x] not in missing_and_gaps_set and (donor_parent_seq[x] not in missing_and_gaps_set or donor_daughter_seq[x] not in missing_and_gaps_set):
					numdiffs+=1
					improvement+=1
					if getposns:
						mutposns.append(x)
					if getimprovement and "SNP_locations" in treeObject.node(node1).get_data().comment and x in treeObject.node(node1).get_data().comment["SNP_locations"]:
						improvement-=1
		
		if getposns:
			return mutposns
		elif getimprovement:
			return improvement
		else:
			return numdiffs
	
	
	
	
	
	def refine_blocks(recipient, donor, snpposns, tempblocks, cutoff, lennogaps, minsnps, N, C, snplocations):
	
			
		newblocks=[]
		
		for block in tempblocks:
			
			
			#trim to first and last SNPs in block
			start=int(block[1])
			end=int(block[2])

			
			
			x=0
			
			while x<(len(snpposns)-1) and snpposns[x]<start:
				x+=1
				
			start=snpposns[x]
			snpposnstart=x			
			
			
			x=len(snpposns)-1
			
			while x>0 and snpposns[x]>end:
				x-=1
				
			end=snpposns[x]
			snpposnend=x

	
			if start>end:
				continue

			
			homoplasies=(snpposnend+1)-snpposnstart
			mutations=get_mutations_between_nodes(recipient, donor, nogapposntoseqposn[recipient][donor][start], nogapposntoseqposn[recipient][donor][end]+1, snplocations, getimprovement=True)
			#print start, end, recipient, donor, mutations, homoplasies
			
			if mutations>=homoplasies:
				continue
			
			old_like=get_block_likelihood((end+1)-start, homoplasies-mutations, N, C)
			toprint={}
			toprint["start"]=[map(str,[old_like, (end+1)-start, (snpposnend+1)-snpposnstart, homoplasies, mutations])]
			printlist=["start", "a","aend","b","bend","c","cend","d","dend","end"]
			newstart=start
			oldstart=start
			newhomoplasies=homoplasies
			newmutations=mutations
			newsnpposnstart=snpposnstart
			
			if snpposnstart+1<(snpposnend+1)-minsnps:
				toprint["a"]=[]
				for y, snpposnloc in enumerate(snpposns[snpposnstart+1: (snpposnend+1)-minsnps]):
					x=y+1

					newhomoplasies-=1
					#newmutations-=get_mutations_between_nodes(recipient, donor, newstart, snpposnloc)
					newmutations=get_mutations_between_nodes(recipient, donor, nogapposntoseqposn[recipient][donor][snpposnloc], nogapposntoseqposn[recipient][donor][end]+1, snplocations, getimprovement=True)
					toprint["a"].append(map(str,[ recipient, donor, snpposnloc, end+1, (end+1)-snpposnloc, newhomoplasies, newmutations]))
					if newmutations>=newhomoplasies:
						break
					
					new_like=get_block_likelihood((end+1)-snpposnloc, newhomoplasies-newmutations, N, C)

					
					if new_like>old_like:
						break
					else:
						newstart=snpposnloc
						newsnpposnstart=snpposnstart+x
						old_like=new_like
						homoplasies=newhomoplasies
						mutations=newmutations
			
				start=newstart
				snpposnstart=newsnpposnstart
			toprint["aend"]=[map(str,[x, start, end, snpposnstart, snpposnend, (snpposnend+1)-snpposnstart, old_like, homoplasies, mutations])]
			
			if start==oldstart and snpposnstart>0:
				toprint["b"]=[]
				newstart=start
				newhomoplasies=homoplasies
				newmutations=mutations
				for y, snpposnloc  in enumerate( snpposns[snpposnstart-1:0:-1]):
					x=y+1
					
					
					newhomoplasies+=1
					
					new_like=get_block_likelihood((end+1)-snpposnloc, newhomoplasies-newmutations, N, C)
						
					
					if new_like>old_like:
						break
						
					
					#newmutations+=get_mutations_between_nodes(recipient, donor, snpposnloc, newstart)
					newmutations=get_mutations_between_nodes(recipient, donor, nogapposntoseqposn[recipient][donor][snpposnloc], nogapposntoseqposn[recipient][donor][end]+1, snplocations, getimprovement=True)
					toprint["b"].append(map(str,[ recipient, donor, snpposnloc, end+1, (end+1)-snpposnloc, newhomoplasies, newmutations]))
					if newmutations>=newhomoplasies:
						break
					
					new_like=get_block_likelihood((end+1)-snpposnloc, newhomoplasies-newmutations, N, C)
					
					if new_like>old_like:
						break
					else:
						newstart=snpposnloc
						newsnpposnstart=snpposnstart-x
						old_like=new_like
						homoplasies=newhomoplasies
						mutations=newmutations
			
				start=newstart
				snpposnstart=newsnpposnstart
			toprint["bend"]=[map(str,[start, end, snpposnstart, snpposnend, (snpposnend+1)-snpposnstart, old_like, homoplasies, mutations])]
			
			newsnpposnend=snpposnend
			newend=end
			oldend=end
			newhomoplasies=homoplasies
			newmutations=mutations
			if snpposnend>=minsnps and (snpposnend-1)>(snpposnstart+minsnps):
				toprint["c"]=[]
				for y, snpposnloc in enumerate( snpposns[ snpposnend-1: snpposnstart+minsnps:   -1]):
					x=y+1
					
					
					newhomoplasies-=1
					#newmutations-=get_mutations_between_nodes(recipient, donor, snpposnloc+1, newend+1)
					newmutations=get_mutations_between_nodes(recipient, donor, nogapposntoseqposn[recipient][donor][start], nogapposntoseqposn[recipient][donor][snpposnloc]+1, snplocations, getimprovement=True)
					toprint["c"].append(map(str,[ recipient, donor, start, snpposnloc+1, (snpposnloc+1)-start, newhomoplasies, newmutations]))
					
					if newmutations>=newhomoplasies or newhomoplasies<1:
						break
					
					new_like=get_block_likelihood((snpposnloc+1)-start, newhomoplasies-newmutations, N, C)
					
										
					if new_like>old_like:
						break
					else:
						newend=snpposnloc
						newsnpposnend=snpposnend-x
						old_like=new_like
						homoplasies=newhomoplasies
						mutations=newmutations
			
				end=newend
				snpposnend=newsnpposnend
				
			toprint["cend"]=[map(str,[start, end, snpposnstart, snpposnend, old_like, homoplasies, mutations])]
			
			if end==oldend and snpposnend+1<len(snpposns):
				toprint["d"]=[]
				newend=end
				newhomoplasies=homoplasies
				newmutations=mutations
				for y, snpposnloc in enumerate(snpposns[ snpposnend+1:]):
					x=y+1
					
					newhomoplasies+=1
					
					new_like=get_block_likelihood((snpposnloc+1)-start, newhomoplasies-newmutations, N, C)
					
					if new_like>old_like:
						break
					
					#newmutations+=get_mutations_between_nodes(recipient, donor, newend+1, snpposnloc+1)
					newmutations=get_mutations_between_nodes(recipient, donor, nogapposntoseqposn[recipient][donor][start], nogapposntoseqposn[recipient][donor][snpposnloc]+1, snplocations, getimprovement=True)
					toprint["d"].append(map(str,[ recipient, donor, start, snpposnloc+1, (snpposnloc+1)-start, newhomoplasies, newmutations]))
					
					if newmutations>=newhomoplasies:
						break
					
					new_like=get_block_likelihood((snpposnloc+1)-start, newhomoplasies-newmutations, N, C)
					
					
					if new_like>old_like:
						break
					else:
						newend=snpposnloc
						newsnpposnend=snpposnend+x
						old_like=new_like
						homoplasies=newhomoplasies
						mutations=newmutations
			
				end=newend
				snpposnend=newsnpposnend
			toprint["dend"]=[map(str,[start, end, snpposnstart, snpposnend, (snpposnend+1)-snpposnstart, old_like, homoplasies, mutations])]
				
			likelihood=old_like
			#print "g", start, end, snpposnstart, snpposnend, old_like, homoplasies, mutations
			#print "new likelihood=", old_like, "start=", start, "end=", end
			
			#print old_like, (end+1)-start, (snpposnend+1)-snpposnstart, homoplasies-mutations
			
			toprint["end"]=[map(str,[old_like, (end+1)-start, (snpposnend+1)-snpposnstart, homoplasies-mutations, homoplasies, mutations])]
			
			if (homoplasies-mutations)>(snpposnend+1)-snpposnstart:
				for f in printlist:
					if f in toprint:
						for y in toprint[f]:
							print f, "\t".join(toprint[y])
			
			snpcount=(snpposnend+1)-snpposnstart
			#print start, end, recipient, donor, mutations, homoplasies, likelihood
			if snpcount>=minsnps:
				newblocks.append([likelihood, start, end, snpposnstart, snpposnend, N, C])
			
			#if recipient in [2, 83] and donor in [2, 83]:
			if start>end:
				print donor, "->", recipient
				for x in printlist:
					if x in toprint:
						for y in toprint[x]:
							print x, "\t".join(y)
	
			
		return newblocks

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	
	
	################################################################
	# Function to detect recombination regions using moving window #
	################################################################
	
	
	def detect_blocks_using_moving_windows_new(recipient, donor, snpposns, lennogaps, minsnps, N, C, snplocations):
		
		totalsnps=float(len(snpposns))
		
		if totalsnps<minsnps:
			return []
		
		windowcounts=[0]*(int(lennogaps))
		
		#window=int(float(lennogaps)/(totalsnps/10))
		window=int(float(N)/(float(C)/10))
		if window<100:
			window=100
		elif window>lennogaps:
			window=lennogaps/2
		if window>10000:
			window=10000
		if window<100:
			window=100
		
		#starttime=time.clock()	
		
		
		#N=float(lennogaps)
		

		#print totalsnps, window, N, C

		#add snp posistions and midway points between two snps to measuredlocs, with dict of those that are real snp posns
		snpwindowposns={}
		measuredlocs=[]
		snpposnstomeasuredlocs={}
		intlocs=[]
		for x, position in enumerate(snpposns):
			snpwindowposns[position]=0
			snpposnstomeasuredlocs[x]=len(measuredlocs)
			measuredlocs.append(position)
			if x+1<len(snpposns):
				interposn=position+((snpposns[x+1]-position)/2)
				if not interposn in [position, snpposns[x+1]]:
					snpwindowposns[interposn]=0
					measuredlocs.append(interposn)
					intlocs.append(interposn)
		
		#
		halfwindow=window/2
		
		for x, position in enumerate(snpposns):
			y=snpposnstomeasuredlocs[x]
			
			while y<len(measuredlocs) and measuredlocs[y]<position+halfwindow:
				snpwindowposns[measuredlocs[y]]+=1
				y+=1
			
			if y==len(measuredlocs) and (lennogaps-position)<halfwindow:
				y=0
				newposition=halfwindow-(lennogaps-position)
				while y<len(measuredlocs) and snpposns[y]<newposition:
					snpwindowposns[measuredlocs[y]]+=1
					y+=1
			
			y=snpposnstomeasuredlocs[x]-1
			
			while y>=0 and measuredlocs[y]>position-halfwindow:
				snpwindowposns[measuredlocs[y]]+=1
				y-=1
			
			if y==-1 and (halfwindow-position)<0:
				y=len(measuredlocs)-1
				newposition=lennogaps-(halfwindow-position)
				while y>=0 and measuredlocs[y]>newposition:
					snpwindowposns[measuredlocs[y]]+=1
					y-=1
			
		#Bit to account for mutations in moving window
		#print snpwindowposns, measuredlocs, snpposns
		intloc=-1
		surroundingsnps=-1
		locvalue=0
		for x, loc in enumerate(measuredlocs):
			
			
			if not loc in intlocs:
				if surroundingsnps==-1:
					surroundingsnps=loc
				elif intloc!=-1 and locvalue>0:
					nummuts=get_mutations_between_nodes(recipient, donor, nogapposntoseqposn[recipient][donor][surroundingsnps], nogapposntoseqposn[recipient][donor][loc]+1, snplocations, getposns=False, getimprovement=True)
					snpwindowposns[measuredlocs[intloc]]-=nummuts
					if snpwindowposns[measuredlocs[intloc]]<0:
						snpwindowposns[measuredlocs[intloc]]=0
					intloc=-1
					surroundingsnps=loc
				elif intloc!=-1 or loc==surroundingsnps+1:
					intloc=-1
					surroundingsnps=loc
				else:
					print "Impossible"
					print intloc, nummuts, surroundingsnps
					for x in measuredlocs:
						print x, snpwindowposns[x]
					print intlocs
					sys.exit()
			else:
				intloc=x
				locvalue=snpwindowposns[measuredlocs[x]]
		
		#print snpwindowposns
		#sys.exit()
				
#		print time.clock()-starttime
#		
		#print snpwindowposns
#		
#		sys.exit()
		
#		threshold=1-(0.05/(float(lennogaps)/(window/10)))
#		cutoff=0
#		pvalue=0.0
#		while pvalue<=threshold:
#			part1=reduce_factorial(window,cutoff)-reduce_factorial(cutoff,cutoff)
#			part2=math.log((float(totalsnps)/lennogaps),10)*cutoff
#			part3=math.log((1.0-(float(totalsnps)/lennogaps)),10)*(window-cutoff)
#			
#			logthing=part1 + part2 + part3
#			
#			pvalue+=10**logthing
#			
#			cutoff+=1


#		threshold=1-(0.05/(float(N)/(window/10)))
#		cutoff=0
#		pvalue=0.0
#		while pvalue<=threshold:
#			part1=reduce_factorial(window,cutoff)-reduce_factorial(cutoff,cutoff)
#			part2=math.log((C/N),10)*cutoff
#			part3=math.log((1.0-(C/N)),10)*(window-cutoff)
#			
#			logthing=part1 + part2 + part3
#			
#			pvalue+=10**logthing
#			
#			cutoff+=1
#
#			
#			
#		cutoff-=1
#		#print cutoff
#		if cutoff<minsnps:
#			cutoff=minsnps
		#print cutoff
		cutoff=minsnps
		newblocks=[]
		if (cutoff<=totalsnps and cutoff>=minsnps) or minsnps<=totalsnps:

#			tempblocks=[]
#			blockstart=0
#			inblock=False
#		
#			x=0
#			y=0
#			
#			while x<len(windowcounts):
#			#for x, pvalue in enumerate(pvalues):
#				#if count>0:
#				#	print >> output, nongapposns[x]+1, windowcounts[x], X2[x], pvalues[x]
#				value=windowcounts[x]
#		
#				
#				#print x, pvalue, inblock
#				
#		#		if value > cutoff:
#		#			print x, value
#				
#				if value>cutoff and not inblock:
#					
#					blockstart=x+1
#					inblock=True
#				elif value<=cutoff and inblock:
#						
#					tempblocks.append([0.0,blockstart,x+1])
#					inblock=False
#			
#				x=x+1
#				
#			if inblock:
#				tempblocks.append([0.0,blockstart,len(windowcounts)])
		
			tempblocks=[]
			blockstart=0
			inblock=False
			
			for x, position in enumerate(measuredlocs):
				if not inblock and snpwindowposns[position]>cutoff:
					blockstart=position
					inblock=True
				elif inblock and snpwindowposns[position]<=cutoff:
					tempblocks.append([0.0,blockstart,measuredlocs[x-1]])
					inblock=False
			if inblock:
				tempblocks.append([0.0,blockstart,measuredlocs[-1]])
			
#			if recipient==47:
#				print "a", tempblocks
			#print tempblocks
			newblocks=refine_blocks(recipient, donor, snpposns, tempblocks, cutoff, lennogaps, minsnps, N, C, snplocations)
			#print newblocks
#			if recipient==47:
#				print "b", newblocks
		
		return newblocks

	
	
	
	def update_node_donor_homoplasies(recipient, donor, snplocations, block, locations, remround, snpcount):
		
		
		def remove_SNP(location, node):
			snpnodedata=treeObject.node(node).get_data()
			if "SNP_locations" in snpnodedata.comment and location in snpnodedata.comment["SNP_locations"]:
			
#				if location==908:
#					print "HERE!!!", location, node, parent_base, daughter_base, donor_parent_base, donor_base
			

				if not "old_SNP_locations" in snpnodedata.comment:
					snpnodedata.comment["old_SNP_locations"]={}
					
				#snpnodedata.comment["sequence"][alignmentbase].parent_base=parentbase
				
				snpnodedata.comment["old_SNP_locations"][location]=snpnodedata.comment["SNP_locations"][location]
				#nodedata=treeObject.node(recipient).get_data()
				if not snpnodedata.comment["old_SNP_locations"][alignmentsnplocation].recombination:
					
					if node==recipient:
						snpnodedata.comment["old_SNP_locations"][alignmentsnplocation].recombination={}
						snpnodedata.comment["old_SNP_locations"][alignmentsnplocation].recombination["donor"]=donor
						snpnodedata.comment["old_SNP_locations"][alignmentsnplocation].recombination["removal round"]=remround
					snpnodedata.comment["old_SNP_locations"][alignmentsnplocation].homoplasy=False
					snpnodedata.comment["old_SNP_locations"][alignmentsnplocation].homoplasies=[]
					
				#treeObject.node(recipient).set_data(nodedata)
				del snpnodedata.comment["SNP_locations"][location]
				#delcount+=1
				if node==recipient and not alignmentsnplocation in changelocations:
					changelocations.append(alignmentsnplocation)
					
				
			elif "old_SNP_locations" in snpnodedata.comment and location in snpnodedata.comment["old_SNP_locations"]:
				print "Already removed and gone ", location, node, parent_base, daughter_base, donor_parent_base, donor_base
			
#			else:
#				print "EH?:", location, node
				#sys.exit()
			
			treeObject.node(node).set_data(snpnodedata)
			
			
			
		def add_SNP(location, parent, node, parentbase, SNPbase):
			snpnodedata=treeObject.node(node).get_data()
			if "SNP_locations" in snpnodedata.comment and location in snpnodedata.comment["SNP_locations"]:
				print "Already a SNP here!!!", location, node
				sys.exit()
			
			if not "SNP_locations" in snpnodedata.comment:
				snpnodedata.comment["SNP_locations"]={}
			
			
			SNP=Si_SNPs_temp.SNP()
			SNP.position=location
			SNP.parent=parent
			SNP.daughter=node
			SNP.parent_base=parentbase
			SNP.daughter_base=SNPbase
			
			
			snpnodedata.comment["SNP_locations"][location]=SNP
			
			snpnodedata.comment["SNP_locations"][location].addrecombination={}
			snpnodedata.comment["SNP_locations"][location].addrecombination["donor"]=donor
			snpnodedata.comment["SNP_locations"][location].addrecombination["added round"]=remround
			
			treeObject.node(node).set_data(snpnodedata)			
			
			
			
		
		def remove_homoplasies(locations):
			
			def remove_node_homoplasies(node,locations):
			
				nodedata=treeObject.node(node).get_data()
				for location in locations:
					if "SNP_locations" in nodedata.comment and location in nodedata.comment["SNP_locations"] and nodedata.comment["SNP_locations"][location].homoplasy:
						nodedata.comment["SNP_locations"][location].homoplasy=False
						nodedata.comment["SNP_locations"][location].homoplasies=[]
					
						
					
			
				treeObject.node(node).set_data(nodedata)
				daughters=treeObject.node(node).get_succ()
				for daughter in daughters:
					remove_node_homoplasies(daughter,locations)
				
			remove_node_homoplasies(treeObject.root,locations)
			
		
		
		def label_removed_SNPs(node, alignmentbase, remround, parentbase, daughterbase):
			
			nodedata=treeObject.node(node).get_data()
			
			if "SNP_locations" in nodedata.comment and alignmentbase in nodedata.comment["SNP_locations"]:
				#print nodedata.comment["sequence"][alignmentbase], parentbase, daughterbase
				#print nodedata.comment["SNP_locations"][alignmentbase].parent_base, nodedata.comment["SNP_locations"][alignmentbase].daughter_base[0]
				
				if nodedata.comment["SNP_locations"][alignmentbase].parent_base[0]==daughterbase and nodedata.comment["SNP_locations"][alignmentbase].daughter_base[0]==parentbase:
				
					#print alignmentbase, node, remround, nodedata.comment["SNP_locations"][alignmentbase].parent_base, nodedata.comment["SNP_locations"][alignmentbase].daughter_base, "remove this SNP"
					remove_SNP(alignmentbase, node)
					return
				
					#print nodedata.comment["SNP_locations"][alignmentbase].homoplasy
				elif nodedata.comment["SNP_locations"][alignmentbase].parent_base[0]==daughterbase:
					nodedata.comment["SNP_locations"][alignmentbase].parent_base=[parentbase]
					treeObject.node(node).set_data(nodedata)
					return
				
				elif nodedata.comment["SNP_locations"][alignmentbase].parent_base[0]==parentbase and nodedata.comment["SNP_locations"][alignmentbase].daughter_base[0]==daughterbase:
					nodedata.comment["sequence"][alignmentbase]=parentbase
					treeObject.node(node).set_data(nodedata)
				
				
			else:
				#print nodedata.comment["sequence"][alignmentbase], parentbase, daughterbase
				nodedata.comment["sequence"][alignmentbase]=parentbase
				treeObject.node(node).set_data(nodedata)
				#print nodedata.comment["sequence"][alignmentbase], parentbase, daughterbase
			
			daughters=treeObject.node(node).get_succ()
#			if len(daughters)==0:
#				print "terminal"
			for daughter in daughters:
				label_removed_SNPs(daughter, alignmentbase, remround, parentbase, daughterbase)
				
				
		def label_added_SNPs(node, alignmentbase, remround, parentbase, daughterbase, donorbase):
			
			nodedata=treeObject.node(node).get_data()
			
			if "SNP_locations" in nodedata.comment and alignmentbase in nodedata.comment["SNP_locations"]:
				#print node, nodedata.comment["sequence"][alignmentbase], parentbase, daughterbase, donorbase, treeObject.is_terminal(node)
				#print nodedata.comment["SNP_locations"][alignmentbase].parent_base, nodedata.comment["SNP_locations"][alignmentbase].daughter_base[0]
				
				if nodedata.comment["SNP_locations"][alignmentbase].parent_base[0]==daughterbase and nodedata.comment["SNP_locations"][alignmentbase].daughter_base[0]==donorbase:
				
					#print alignmentbase, node, remround, nodedata.comment["SNP_locations"][alignmentbase].parent_base, nodedata.comment["SNP_locations"][alignmentbase].daughter_base, "remove this SNP"
					remove_SNP(alignmentbase, node)
					return
				
					#print nodedata.comment["SNP_locations"][alignmentbase].homoplasy
				elif nodedata.comment["SNP_locations"][alignmentbase].parent_base[0]==daughterbase:
					nodedata.comment["SNP_locations"][alignmentbase].parent_base=[donorbase]
					treeObject.node(node).set_data(nodedata)
					return
				
				elif nodedata.comment["SNP_locations"][alignmentbase].parent_base[0]==parentbase and nodedata.comment["SNP_locations"][alignmentbase].daughter_base[0]==daughterbase:
					nodedata.comment["sequence"][alignmentbase]=donorbase
					treeObject.node(node).set_data(nodedata)
				
				
			else:
				#print nodedata.comment["sequence"][alignmentbase], parentbase, daughterbase
				nodedata.comment["sequence"][alignmentbase]=donorbase
				treeObject.node(node).set_data(nodedata)
				#print nodedata.comment["sequence"][alignmentbase], parentbase, daughterbase
			#print node, nodedata.comment["sequence"][alignmentbase], parentbase, daughterbase, donorbase
			daughters=treeObject.node(node).get_succ()
#			if len(daughters)==0:
#				print "terminal"
			for daughter in daughters:
				label_added_SNPs(daughter, alignmentbase, remround, parentbase, daughterbase, donorbase)
		
		
		
		recipientdata=treeObject.node(recipient).get_data()
		parentdata=treeObject.node(treeObject.node(recipient).get_prev()).get_data()
		if treeObject.common_ancestor(recipient, donor)==donor:
			#print "reversal"
			rtype="c"
			if donor==treeObject.root:
				donorparentdata=treeObject.node(donor).get_data()
			else:
				donorparentdata=treeObject.node(treeObject.node(donor).get_prev()).get_data()
			donordata=treeObject.node(donor).get_data()
		else:
			#print "convergence"
			rtype="c"
			if donor==treeObject.root:
				donorparentdata=treeObject.node(donor).get_data()
			else:
				donorparentdata=treeObject.node(treeObject.node(donor).get_prev()).get_data()
			donordata=treeObject.node(donor).get_data()
				
		
		changelocations=[]
		bases={}
		start=nogapposntoseqposn[recipient][donor][homoplasyposns[recipient][donor][block[3]]]
		end=nogapposntoseqposn[recipient][donor][homoplasyposns[recipient][donor][block[4]]]+1
		#print homoplasyposns[recipient][donor][block[3]:block[4]]
		#print block[4]-block[3]
		
		counta=0
		countb=0
		countc=0
		countd=0
		delcount=0
#		for nogapsnplocation in homoplasyposns[recipient][donor][block[3]:block[4]+1]:
#			
#			alignmentsnplocation=nogapposntoseqposn[recipient][donor][nogapsnplocation]

#		test=get_mutations_between_nodes(recipient, donor, start, end, final=True)
#		print test
		for alignmentsnplocation in xrange(start,end):
			
			if alignmentsnplocation in snplocations:
			
			#if "SNP_locations" in nodedata.comment and alignmentsnplocation in nodedata.comment["SNP_locations"] and nodedata.comment["SNP_locations"][alignmentbase].homoplasy:
			
				parent_base=parentdata.comment["sequence"][alignmentsnplocation]
				daughter_base=recipientdata.comment["sequence"][alignmentsnplocation]
				donor_base=donordata.comment["sequence"][alignmentsnplocation]
				donor_parent_base=donorparentdata.comment["sequence"][alignmentsnplocation]
				bases[alignmentsnplocation]=[parent_base, daughter_base, donor_parent_base, donor_base]
				
				
				if rtype=="r":
					print parent_base, daughter_base, donor_parent_base, donor_base
				
				if parent_base!="-" and daughter_base!="-" and (donor_parent_base!="-" or donor_base!="-"):
				
					if parent_base!=daughter_base and (daughter_base==donor_base or daughter_base==donor_parent_base):
#						if donor_base=="-" or donor_parent_base=="-":
#							print alignmentsnplocation, recipient, donor, "remove this SNP", parent_base, daughter_base, donor_parent_base, donor_base
						label_removed_SNPs(recipient, alignmentsnplocation, remround, parent_base, daughter_base)
						remove_SNP(alignmentsnplocation, recipient)
						
						if not alignmentsnplocation in changelocations:
							counta+=1

#					elif parent_base!=daughter_base and daughter_base!=donor_base and daughter_base!=donor_parent_base:
#						print alignmentsnplocation, "remove this SNP", parent_base, daughter_base, " and add this SNP", parent_base, daughter_base, donor_parent_base, donor_base
#						#label_removed_SNPs(recipient, alignmentsnplocation, remround, parent_base, daughter_base)
#						if not alignmentsnplocation in changelocations:
#							countc+=1
#							#changelocations.append(alignmentsnplocation)
					elif daughter_base!=donor_base and daughter_base!=donor_parent_base and parent_base!=donor_base and parent_base!=donor_parent_base and donor_parent_base!="-" and donor_base!="-":
						if "SNP_locations" in recipientdata.comment and alignmentsnplocation in recipientdata.comment["SNP_locations"]:
							if rtype=="r":
								print alignmentsnplocation, "keep this SNP", parent_base, daughter_base, donor_parent_base, donor_base
						else:
							#print alignmentsnplocation, "add this SNP", parent_base, daughter_base, donor_parent_base, donor_base
							
							label_added_SNPs(recipient, alignmentsnplocation, remround, parent_base, daughter_base, donor_base)
							add_SNP(alignmentsnplocation, treeObject.node(recipient).get_prev(), recipient, parent_base, donor_base)
								
							if not alignmentsnplocation in changelocations:
								countd+=1
								#changelocations.append(alignmentsnplocation)
				parent_base=parentdata.comment["sequence"][alignmentsnplocation]
				daughter_base=recipientdata.comment["sequence"][alignmentsnplocation]
				donor_base=donordata.comment["sequence"][alignmentsnplocation]
				donor_parent_base=donorparentdata.comment["sequence"][alignmentsnplocation]
					
		#print len(changelocations), counta, countb, countc, countd, delcount
		
		if len(changelocations)!=snpcount:
			print len(changelocations), counta, countb, countc, countd, delcount
			print nogapposntoseqposn[recipient][donor]
			print seqposntonogapposn[recipient][donor]
			print homoplasyposns[recipient][donor]
			
#			for nogapsnplocation in homoplasyposns[recipient][donor][block[3]:block[4]+1]:
#			
#				alignmentsnplocation=nogapposntoseqposn[recipient][donor][nogapsnplocation]
			for alignmentsnplocation in xrange(start,end):
				try:
					print alignmentsnplocation, bases[alignmentsnplocation]
				except StandardError:
					print alignmentsnplocation
				parent_base=parentdata.comment["sequence"][alignmentsnplocation]
				daughter_base=recipientdata.comment["sequence"][alignmentsnplocation]
				donor_base=donordata.comment["sequence"][alignmentsnplocation]
				donor_parent_base=donorparentdata.comment["sequence"][alignmentsnplocation]
				if alignmentsnplocation in changelocations:
					inchange=True
				else:
					inchange=False
				print donor, recipient, parent_base, daughter_base, donor_parent_base, donor_base, inchange
#				if alignmentsnplocation in snplocations:
#					numdiff=0
#					types=[]
#					for x in bases[alignmentsnplocation]:
#						if x!="-" and x not in types:
#							numdiff+=1
#							types.append(x)
#					if numdiff>1:		
#						if alignmentsnplocation in nogapposntoseqposn[recipient][donor].values():
#							print "SNPhere"
#						else:
#							print "noSNPhere"
#						print alignmentsnplocation, bases[alignmentsnplocation]
#						parent_base=parentdata.comment["sequence"][alignmentsnplocation]
#						daughter_base=recipientdata.comment["sequence"][alignmentsnplocation]
#						donor_base=donordata.comment["sequence"][alignmentsnplocation]
#						donor_parent_base=donorparentdata.comment["sequence"][alignmentsnplocation]
#						if alignmentsnplocation in changelocations:
#							inchange=True
#						else:
#							inchange=False
#						print parent_base, daughter_base, donor_parent_base, donor_base, inchange
			sys.exit()
		
		remove_homoplasies(snplocations)
		identify_homoplasies(treeObject, snplocations, original=False)
		for node in homoplasyposns.keys():
			del homoplasyposns[node]
		for node in seqposntonogapposn.keys():
			del seqposntonogapposn[node]
		for node in nogapposntoseqposn.keys():
			del nogapposntoseqposn[node]
		node_donor_homoplasies(treeObject.root, locations=snplocations)
		
		#sys.exit()
		
	
	
#	def update_node_donor_homoplasies_old(recipient, donor, snplocations, block, locations, remround):
#	
#		
#		
#		
#		#function to find all nodes downstream of the recipient node which share the same snp as a result of the recombination
#		
#		def get_recipient_nodes(recipientnode,recipientlist, alignmentbase):
#		
#			daughters=treeObject.node(recipientnode).succ
#		
#			for daughter in daughters:
#				if not treeObject.node(daughter).get_data().comment.has_key("SNP_locations") or not alignmentbase in treeObject.node(daughter).get_data().comment["SNP_locations"]:
#					recipientlist=get_recipient_nodes(daughter,recipientlist, alignmentbase)
#			recipientlist.append(recipientnode)
#			return recipientlist
#		
#		
#		
#		def find_ancestral_donor_node(donornode, alignmentbase):
#			
#			node=donornode
#			
#			while node!=treeObject.root and ( not treeObject.node(node).get_data().comment.has_key("SNP_locations") or not alignmentbase in treeObject.node(node).get_data().comment["SNP_locations"]):
#				node=treeObject.node(node).get_prev()
#			#print node, alignmentbase
#			return node
#		
#		
#		
#		#function to find all nodes downstream of the recipient node which share the same snp as a result of the recombination
#		
#		def get_donor_nodes(donornode,donorlist, alignmentbase):
#			node=find_ancestral_donor_node(donornode, alignmentbase)
#		
#			donorlist=get_recipient_nodes(node,[], alignmentbase)
#			
#			return donorlist
#		
#		#function to remove a homoplasy from a node
#		
#		def remove_donor_homoplasy(recipientnode, donornode, alignmentbase):
#			translatedsnplocation=seqposntonogapposn[recipientnode][donornode][alignmentbase]
#			#print "C) removing recipient ", recipientnode, donornode, alignmentbase, translatedsnplocation
#			homoplasyposns[recipientnode][donornode].remove(translatedsnplocation)
#			
#			del seqposntonogapposn[recipientnode][donornode][alignmentbase]
#			del nogapposntoseqposn[recipientnode][donornode][translatedsnplocation]
#			
#			nodedata=treeObject.node(recipientnode).get_data()
#			if "SNP_locations" in nodedata.comment and alignmentbase in nodedata.comment["SNP_locations"] and nodedata.comment["SNP_locations"][alignmentbase].homoplasy:
#				for hom in nodedata.comment["SNP_locations"][alignmentbase].homoplasies:
#					if hom[1]==donornode:
#						nodedata.comment["SNP_locations"][alignmentbase].oldhomoplasies.append(hom)
#						nodedata.comment["SNP_locations"][alignmentbase].homoplasies.remove(hom)
#						#print "A", hom, alignmentbase
#					if len(nodedata.comment["SNP_locations"][alignmentbase].homoplasies)==0:
#						nodedata.comment["SNP_locations"][alignmentbase].homoplasy=False
#						nodedata.comment["SNP_locations"][alignmentbase].oldhomoplasy=True
#			treeObject.node(recipientnode).set_data(nodedata)
#			
#		
#		
#		def remove_recipient_homoplasy(recipientnode, donornode, alignmentbase):
#			translatedsnplocation=seqposntonogapposn[recipientnode][donornode][alignmentbase]
#			#print "C) removing recipient ", recipientnode, donornode, alignmentbase, translatedsnplocation
#			homoplasyposns[recipientnode][donornode].remove(translatedsnplocation)
#			
#			del seqposntonogapposn[recipientnode][donornode][alignmentbase]
#			del nogapposntoseqposn[recipientnode][donornode][translatedsnplocation]
#			
#			nodedata=treeObject.node(recipientnode).get_data()
##			if "SNP_locations" in nodedata.comment:
##				print nodedata.comment["SNP_locations"], alignmentbase
#			
#			if "SNP_locations" in nodedata.comment and alignmentbase in nodedata.comment["SNP_locations"] and nodedata.comment["SNP_locations"][alignmentbase].homoplasy:
#				for hom in nodedata.comment["SNP_locations"][alignmentbase].homoplasies:
#					if hom[1]==donornode:
#						nodedata.comment["SNP_locations"][alignmentbase].oldhomoplasies.append(hom)
#						nodedata.comment["SNP_locations"][alignmentbase].homoplasies.remove(hom)
##						print "A", hom, alignmentbase
#					if len(nodedata.comment["SNP_locations"][alignmentbase].homoplasies)==0:
#						nodedata.comment["SNP_locations"][alignmentbase].homoplasy=False
#						nodedata.comment["SNP_locations"][alignmentbase].oldhomoplasy=True
#			treeObject.node(recipientnode).set_data(nodedata)
#			
#			
#			nodedata=treeObject.node(donornode).get_data()
#			if "SNP_locations" in nodedata.comment and alignmentbase in nodedata.comment["SNP_locations"] and nodedata.comment["SNP_locations"][alignmentbase].homoplasy:
#				for hom in nodedata.comment["SNP_locations"][alignmentbase].homoplasies:
#					if hom[1]==recipientnode:
#						nodedata.comment["SNP_locations"][alignmentbase].oldhomoplasies.append(hom)
#						nodedata.comment["SNP_locations"][alignmentbase].homoplasies.remove(hom)
##						print "B", hom, alignmentbase
#					if len(nodedata.comment["SNP_locations"][alignmentbase].homoplasies)==0:
#						nodedata.comment["SNP_locations"][alignmentbase].homoplasy=False
#						nodedata.comment["SNP_locations"][alignmentbase].oldhomoplasy=True
#			treeObject.node(donornode).set_data(nodedata)
#			
#			
#			
#		
#		#function to remove a homoplasy from a donor
#		
#		def remove_recipient(recipientnode,alignmentbase):
#		
#			for donornode in homoplasyposns[recipientnode]:
#				if alignmentbase in seqposntonogapposn[recipientnode][donornode]:
#					remove_recipient_homoplasy(recipientnode, donornode, alignmentbase)
#					
#		#function to remove a homoplasy from a donor
#		
#		def remove_donor(recipientlist, donorlist, alignmentbase):
#			
#			for donornode in donorlist:
#				for recipientnode in  recipientlist:
##					print recipientnode, donornode, alignmentbase
##					if donornode in homoplasyposns:
##						print "D", seqposntonogapposn[donornode]
##					if recipientnode in homoplasyposns:
##						print "R", seqposntonogapposn[recipientnode]
#					if donornode in homoplasyposns and seqposntonogapposn[donornode].has_key(recipientnode) and seqposntonogapposn[donornode][recipientnode].has_key(alignmentbase):
#						remove_donor_homoplasy(donornode, recipientnode, alignmentbase)
#					
#					if recipientnode in homoplasyposns and seqposntonogapposn[recipientnode].has_key(donornode) and seqposntonogapposn[recipientnode][donornode].has_key(alignmentbase):
#						remove_donor_homoplasy(recipientnode, donornode, alignmentbase)
#		
#		
#
#		
#		recipientdata=treeObject.node(recipient).get_data()
#		recipientseq=recipientdata.comment["sequence"]
#		donorseq=treeObject.node(donor).get_data().comment["sequence"]
#		if donor==treeObject.root:
#			donorparentseq=donorseq
#			donorparent=donor
#		else:
#			donorparentseq=treeObject.node(treeObject.node(donor).get_prev()).get_data().comment["sequence"]
#			donorparent=treeObject.node(donor).get_prev()
#		mutationstoadd=set()
#		
#		count=0
#
#		
#		for alignmentsnplocation in xrange(nogapposntoseqposn[recipient][donor][homoplasyposns[recipient][donor][block[3]]],nogapposntoseqposn[recipient][donor][homoplasyposns[recipient][donor][block[4]]]+1):
#			
#			
#			
#			if alignmentsnplocation in recipientdata.comment["SNP_locations"] and (recipientseq[alignmentsnplocation]==donorseq[alignmentsnplocation] or recipientseq[alignmentsnplocation]==donorparentseq[alignmentsnplocation]) and recipientseq[alignmentsnplocation] not in missing_and_gaps_set:
#				count+=1
#				#print recipientseq[alignmentsnplocation], donorseq[alignmentsnplocation], donorparentseq[alignmentsnplocation]
#
#				#print homoplasyposns.keys(), seqposntonogapposn.keys()
#				
#				recipientlist=get_recipient_nodes(recipient,[], alignmentsnplocation)
#				recipientlist.sort()
#				print recipientlist
#				
#				if recipientseq[alignmentsnplocation]==donorseq[alignmentsnplocation]:
#					donornode=donor
#				elif recipientseq[alignmentsnplocation]==donorparentseq[alignmentsnplocation]:
#					donornode=donorparent
#				
#				donorlist=get_donor_nodes(donornode,[], alignmentsnplocation)
#				donorlist.sort()
#				print alignmentsnplocation, donorlist
#				
#				remove_donor(recipientlist, donorlist, alignmentsnplocation)
#				remove_recipient(recipient,alignmentsnplocation)
#				nodedata=treeObject.node(recipient).get_data()
#				if not nodedata.comment["SNP_locations"][alignmentsnplocation].recombination:
#					nodedata.comment["SNP_locations"][alignmentsnplocation].recombination={}
#					nodedata.comment["SNP_locations"][alignmentsnplocation].recombination["donor"]=donor
#					nodedata.comment["SNP_locations"][alignmentsnplocation].recombination["removal round"]=remround
#				treeObject.node(recipient).set_data(nodedata)

	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	def detect_homoplasy_clusters_using_moving_windows_new(snplocations):
		
	
		def calculate_N():
			#nnode=treeObject.root
			
			def node_N(nnode, N):
				
				ndaughters=treeObject.node(nnode).get_succ()
				for ndaughter in ndaughters:
					N=node_N(ndaughter, N)
				
				comparitors=get_nodes_not_downstream(treeObject, nnode)
				for comparitor in comparitors:
					N+=len(str(treeObject.node(nnode).get_data().comment["sequence"]).replace("N","").replace("-",""))
				#N+=len(str(treeObject.node(nnode).get_data().comment["sequence"]).replace("N","").replace("-",""))
				return N
				
			N=0.0
			N=node_N(treeObject.root, N)
			return N
							


		
		if len(homoplasyposns)!=len(seqposntonogapposn) or len(homoplasyposns)!=len(alignmentlengths):
			return
		
		Final_blocks=[]
		added=True
		
		roundnum=1
		while added:
			
#			if roundnum==6:
#				return Final_blocks
			
			
			print " Round", roundnum, "...",
			sys.stdout.flush(),
			starttime=time.clock()
			if roundnum==1:
				N=calculate_N()
			added=False
			blocks=[0]
			
			
		
			blocks=[]
			C=0.0
			for node in homoplasyposns.keys():
				for donor in homoplasyposns[node]:
			#		if len(homoplasyposns[node][donor])>0:
			#			C+=1
					C+=len(homoplasyposns[node][donor])

			#	C+=len(homoplasyposns[node])
				for donor in homoplasyposns[node].keys():
					if len(homoplasyposns[node][donor])==0:
						del homoplasyposns[node][donor]
						del seqposntonogapposn[node][donor]
						del nogapposntoseqposn[node][donor]
				if len(homoplasyposns[node])==0:
					del homoplasyposns[node]
					del seqposntonogapposn[node]
					del nogapposntoseqposn[node]
			#print homoplasyposns
			
			#sys.exit()
			#print C, "sites remaining in homoplasy database.",
			#print N, C
			#print homoplasyposns
			for node in homoplasyposns:
				
		
				for donornode in homoplasyposns[node]:
			
	
					donor_blocks=[]
					#startbtime= time.clock()
#						if node==47:	
#							print node, donornode, len(homoplasyposns[node][donornode])
					donor_blocks=detect_blocks_using_moving_windows_new(node, donornode, homoplasyposns[node][donornode], alignmentlengths[node][donornode], minsnps, N, C, snplocations)
					
#						if node==47 and len(donor_blocks)>0:
#							print donornode, donor_blocks
					
					#print "b", time.clock()-startbtime
					#blocks are of structure [ll,start posn, end posn, homoplasy list start homoplasy list end]
					for x in xrange(0,len(donor_blocks)):
						donor_blocks[x].append(node)
						donor_blocks[x].append(donornode)
						
						
					blocks=blocks+donor_blocks
				
			
			
			blocks.sort()
			while not added and len(blocks)>0:
				#print blocks[:3]
				#sys.exit()
				if len(blocks)>0:
#					print homoplasyposns[47]
#					for block in blocks:
#						if block[7]==47:
#							print block
					#sys.exit()
					minloglike=blocks[0][0]
					
					recipients={blocks[0][7]:[blocks[0][8]]}
					
					x=1
					while x<len(blocks) and blocks[x][0]==minloglike:
						
						if not blocks[x][7] in recipients:
							recipients[blocks[x][7]]=[blocks[x][8]]
						else:
							recipients[blocks[x][7]].append(blocks[x][8])
					
						x+=1
					
					blocks_with_min_ll=x
					
					recipient = recipients.keys()[0]
					donor=recipients[recipient][0]
					minhomoplasynum=len(homoplasyposns[recipient][donor])
					
#go through all possible recipients and choose the one with the fewest homoplasies in total. In effect this is choosing the most likely recipient when the recombination direction is not possible to estimate from the number of homoplasies in the block
					for r in recipients:
						recipienthomoplasyposns=0
						for d in homoplasyposns[r]:
							recipienthomoplasyposns+=len(homoplasyposns[r][d])
							
						if recipienthomoplasyposns<minhomoplasynum:
						       	minhomoplasynum=recipienthomoplasyposns
						       	recipient=r
						       	donor=d



					
					
					
					common_ancestor=get_common_ancestor_of_taxon_list(treeObject,recipients[recipient])
					alldonors=recipients[recipient]
					alldonors.sort()
					
					
					
					oldtotalsnps=len(homoplasyposns[recipient][donor])
					
					
					
					#first identify which block was the one for the correct donor and recipient
					for x in xrange(0,blocks_with_min_ll):
						if blocks[x][7]==recipient and blocks[x][8]==donor:
							blocktoremove=x
							break
					
					start=nogapposntoseqposn[recipient][donor][blocks[blocktoremove][1]]
					end=nogapposntoseqposn[recipient][donor][blocks[blocktoremove][2]]


					if start>end:
						print blocks[blocktoremove], homoplasyposns[recipient][donor], nogapposntoseqposn[recipient][donor]
					
		
					#Calculate the p-value of the block
					
					
					blocklength=(blocks[blocktoremove][2]+1)-blocks[blocktoremove][1]
					lennogaps=alignmentlengths[recipient][donor]
					snpcount=(blocks[blocktoremove][4]+1)-blocks[blocktoremove][3]
					mutationcount=get_mutations_between_nodes(recipient, donor, start, end+1, snplocations, getimprovement=True)
					
					#mutationcount=blocks[blocktoremove][3]
					diffcount=snpcount-mutationcount
					N=blocks[blocktoremove][5]
					C=blocks[blocktoremove][6]
					
					#print start, end, blocklength, snpcount, mutationcount, diffcount
					
					pvalue=1
					if diffcount>=minsnps:
						x=0
						pvalue=0.0
						while x<diffcount:
							
#							part1=reduce_factorial(blocklength,x)-reduce_factorial(x,x)
#							part2=math.log((float(oldtotalsnps)/lennogaps),10)*x
#							part3=math.log((1.0-(float(oldtotalsnps)/lennogaps)),10)*(blocklength-x)
							part1=reduce_factorial(blocklength,x)-reduce_factorial(x,x)
							part2=math.log((C/N),10)*x
							part3=math.log((1.0-(C/N)),10)*(blocklength-x)
						
							logthing=part1 + part2 + part3
							
							pvalue+=10**logthing
							x+=1
						pvalue=1.0-round(pvalue,10)
						pvaluethreshold=(0.05/float(N))
						#pvaluethreshold=0.05
						#print pvalue, pvaluethreshold, N, C, blocklength, diffcount, minsnps, blocks[blocktoremove], snpcount
						
						if pvalue<pvaluethreshold:
							Final_blocks.append({"start":start, "end":end, "ll":blocks[blocktoremove][0], "common_ancestor":common_ancestor, "donors":alldonors, "snpcount": snpcount, "p":pvalue, "removal_round":roundnum, "recipient":recipient})
							roundnum+=1
							added=True
							print "/".join(map(str,alldonors))+"->"+str(recipient)+", "+str(start)+".."+str(end)+", "+str(snpcount)+" snps in "+str(blocklength)+" bases, mutations="+str(mutationcount)+", improvement="+str(diffcount)+" -ll="+str(blocks[blocktoremove][0])+", pvalue="+str(pvalue)+".",
							#bit to remove homoplasies involved in the recombination identified
							update_node_donor_homoplasies(recipient, donor, snplocations, blocks[blocktoremove], locations, roundnum, snpcount)
					
					blocks=blocks[len(recipients):]
					
			if not added:
				print "No blocks found.",
			print time.clock()-starttime, "secs"
		return 	Final_blocks
			
	
	

	homoplasyposns={}
	seqposntonogapposn={}
	nogapposntoseqposn={}
	alignmentlengths={}
	minsnps=3
	
	print " Creating homoplasy database ...",
	sys.stdout.flush()
	node_donor_homoplasies(treeObject.root, locations)
	print "Done"
	sys.stdout.flush()
	
	Final_blocks=detect_homoplasy_clusters_using_moving_windows_new(locations)
	handle=open(prefix+"_recombination_blocks.tab","w")
	
	print "Found", len(Final_blocks), "homoplasy block(s)"
	sys.stdout.flush()
	print_blocks_to_tab(treeObject, Final_blocks, handle)
	
	handle.close()

	return treeObject
