#!/usr/bin/env python

#################################
# Import some necessary modules #
#################################

import string, re
import os, sys
import random
from math import sqrt, pow, log
from numpy import repeat, convolve, mean, median
from optparse import OptionParser, OptionGroup
from Bio.Nexus import Trees, Nodes
import shlex, subprocess
#on my laptop
#sys.path.extend(map(os.path.abspath, ['/Users/sh16/Documents/scripts/modules/']))
#on pcs4
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_general import *
from Si_SeqIO import *
from Si_nexus import midpoint_root, tree_to_string
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import GenBank
from Bio.Graphics.GenomeDiagram._Colors import ColorTranslator
from math import floor
from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.lib import pagesizes
from reportlab.graphics.shapes import * 
from reportlab.pdfgen.canvas import Canvas
from reportlab.graphics.charts.textlabels import Label
from reportlab.graphics.charts.lineplots import LinePlot
from reportlab.graphics.widgets.markers import makeMarker
from reportlab.graphics.charts.axes import XValueAxis, Label
from reportlab.graphics.charts.legends import Legend, LineLegend
from reportlab.graphics.charts.barcharts import VerticalBarChart
from reportlab.graphics.widgets.grids import ShadedRect
from reportlab.graphics.widgets.signsandsymbols import ArrowOne
from reportlab.graphics import renderPDF
import pysam
#on my laptop
#SAMTOOLS_DIR="/Users/sh16/Applications/samtools-0.1.18/"
#BCFTOOLS_DIR="/Users/sh16/Applications/samtools-0.1.18/bcftools/"
#on pcs4
SAMTOOLS_DIR="/nfs/users/nfs_s/sh16/samtools-0.1.17/"
BCFTOOLS_DIR="/nfs/users/nfs_s/sh16/samtools-0.1.17/bcftools/"

#################################
# Set up some global converters #
#################################

pagesizeconverter={'A0':pagesizes.A0, 'A1':pagesizes.A1, 'A2':pagesizes.A2, 'A3':pagesizes.A3, 'A4':pagesizes.A4, 'A5':pagesizes.A5, 'A6':pagesizes.A6, 'B0':pagesizes.B0, 'B1':pagesizes.B1, 'B2':pagesizes.B2, 'B3':pagesizes.B3, 'B4':pagesizes.B4, 'B5':pagesizes.B5, 'B6':pagesizes.B6, 'LEGAL':pagesizes.LEGAL, 'LETTER':pagesizes.LETTER, 'landscape':pagesizes.landscape, 'legal':pagesizes.legal, 'letter':pagesizes.letter, 'portrait':pagesizes.portrait}


colourconverter={'aliceblue':colors.aliceblue, 'antiquewhite':colors.antiquewhite, 'aqua':colors.aqua, 'aquamarine':colors.aquamarine, 'azure':colors.azure, 'beige':colors.beige, 'bisque':colors.bisque, 'black':colors.black, 'blanchedalmond':colors.blanchedalmond, 'blue':colors.blue, 'blueviolet':colors.blueviolet, 'brown':colors.brown, 'burlywood':colors.burlywood, 'cadetblue':colors.cadetblue, 'chartreuse':colors.chartreuse, 'chocolate':colors.chocolate, 'coral':colors.coral, 'cornflower':colors.cornflower, 'cornflowerblue':colors.cornflowerblue, 'cornsilk':colors.cornsilk, 'crimson':colors.crimson, 'cyan':colors.cyan, 'darkblue':colors.darkblue, 'darkcyan':colors.darkcyan, 'darkgoldenrod':colors.darkgoldenrod, 'darkgray':colors.darkgray, 'darkgreen':colors.darkgreen, 'darkgrey':colors.darkgrey, 'darkkhaki':colors.darkkhaki, 'darkmagenta':colors.darkmagenta, 'darkolivegreen':colors.darkolivegreen, 'darkorange':colors.darkorange, 'darkorchid':colors.darkorchid, 'darkred':colors.darkred, 'darksalmon':colors.darksalmon, 'darkseagreen':colors.darkseagreen, 'darkslateblue':colors.darkslateblue, 'darkslategray':colors.darkslategray, 'darkslategrey':colors.darkslategrey, 'darkturquoise':colors.darkturquoise, 'darkviolet':colors.darkviolet, 'deeppink':colors.deeppink, 'deepskyblue':colors.deepskyblue, 'dimgray':colors.dimgray, 'dimgrey':colors.dimgrey, 'dodgerblue':colors.dodgerblue, 'fidblue':colors.fidblue, 'fidlightblue':colors.fidlightblue, 'fidred':colors.fidred, 'firebrick':colors.floralwhite, 'floralwhite':colors.floralwhite, 'forestgreen':colors.forestgreen, 'fuchsia':colors.fuchsia, 'gainsboro':colors.gainsboro, 'ghostwhite':colors.ghostwhite, 'gold':colors.gold, 'goldenrod':colors.goldenrod, 'gray':colors.gray, 'green':colors.green, 'greenyellow':colors.greenyellow, 'grey':colors.grey, 'honeydew':colors.honeydew, 'hotpink':colors.hotpink, 'indianred':colors.indianred, 'indigo':colors.indigo, 'ivory':colors.ivory, 'khaki':colors.khaki, 'lavender':colors.lavender, 'lavenderblush':colors.lavenderblush, 'lawngreen':colors.lawngreen, 'lemonchiffon':colors.lemonchiffon, 'lightblue':colors.lightblue, 'lightcoral':colors.lightcoral, 'lightcyan':colors.lightcyan, 'lightgoldenrodyellow':colors.lightgoldenrodyellow, 'lightgreen':colors.lightgreen, 'lightgrey':colors.lightgrey, 'lightpink':colors.lightpink, 'lightsalmon':colors.lightsalmon, 'lightseagreen':colors.lightseagreen, 'lightskyblue':colors.lightskyblue, 'lightslategray':colors.lightslategray, 'lightslategrey':colors.lightslategrey, 'lightsteelblue':colors.lightsteelblue, 'lightyellow':colors.lightyellow, 'lime':colors.lime, 'limegreen':colors.limegreen, 'linen':colors.linen, 'magenta':colors.magenta, 'maroon':colors.maroon, 'math':colors.math, 'mediumaquamarine':colors.mediumaquamarine, 'mediumblue':colors.mediumblue, 'mediumorchid':colors.mediumorchid, 'mediumpurple':colors.mediumpurple, 'mediumseagreen':colors.mediumseagreen, 'mediumslateblue':colors.mediumslateblue, 'mediumspringgreen':colors.mediumspringgreen, 'mediumturquoise':colors.mediumturquoise, 'mediumvioletred':colors.mediumvioletred, 'midnightblue':colors.midnightblue, 'mintcream':colors.mintcream, 'mistyrose':colors.mistyrose, 'moccasin':colors.moccasin, 'navajowhite':colors.navajowhite, 'navy':colors.navy , 'oldlace':colors.oldlace, 'olive':colors.olive, 'olivedrab':colors.olivedrab, 'orange':colors.orange, 'orangered':colors.orangered, 'orchid':colors.orchid, 'palegoldenrod':colors.palegoldenrod, 'palegreen':colors.palegreen, 'paleturquoise':colors.paleturquoise, 'palevioletred':colors.palevioletred, 'papayawhip':colors.papayawhip, 'peachpuff':colors.peachpuff, 'peru':colors.peru, 'pink':colors.pink, 'plum':colors.plum, 'powderblue':colors.powderblue, 'purple':colors.purple, 'red':colors.red, 'rosybrown':colors.rosybrown, 'royalblue':colors.royalblue, 'saddlebrown':colors.saddlebrown, 'salmon':colors.salmon, 'sandybrown':colors.sandybrown, 'seagreen':colors.seagreen, 'seashell':colors.seashell, 'sienna':colors.sienna, 'silver':colors.silver, 'skyblue':colors.skyblue, 'slateblue':colors.slateblue, 'slategray':colors.slategray, 'slategrey':colors.slategrey, 'snow':colors.snow, 'springgreen':colors.springgreen, 'steelblue':colors.steelblue, 'tan':colors.tan, 'teal':colors.teal, 'thistle':colors.thistle, 'tomato':colors.tomato, 'turquoise':colors.turquoise, 'violet':colors.violet, 'wheat':colors.wheat, 'white':colors.white, 'whitesmoke':colors.whitesmoke, 'yellow':colors.yellow, 'yellowgreen':colors.yellowgreen}

########################################################
# Up the default recursion limit to allow bigger trees #
########################################################
sys.setrecursionlimit(10000)


################################
# Get the command line options #
################################


def main():
	usage = "usage: %prog [options] args"
	parser = OptionParser(usage=usage)
	
	
	group = OptionGroup(parser, "General Output Options")
	
	group.add_option("-o", "--output", action="store", dest="outputfile", help="output file name [default= %default]", type="string", metavar="FILE", default="test.pdf")
	group.add_option("-O", "--orientation", action="store", choices=['landscape', 'portrait'], dest="orientation", help="page orientation [default= %default]", type="choice", default="landscape")
	group.add_option("-p", "--pagesize", action="store", choices=['A0', 'A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'B0', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'LEGAL', 'LETTER', 'legal', 'letter', 'custom'], dest="page", help="page size [default= %default]", type="choice", default="A4")
	group.add_option("-5", "--custompagesize", action="store", dest="custompage", help="custom page size. Must be in the form of x,y in mm. Only applicable with the custom -p option", default="")
	group.add_option("-P", "--npages", action="store", dest="npages", help="number of pages to  split picture over [default= %default]", type="int", default=1)
	
	parser.add_option_group(group)
	
	
	group = OptionGroup(parser, "Tree Output Options")
	
	group.add_option("-t", "--tree", action="store", dest="tree", help="tree file to align tab files to", default="")
	group.add_option("-2", "--proportion", action="store", dest="treeproportion", help="Proportion of page to take up with the tree", default=0.3, type='float')
	group.add_option("-s", "--support", action="store_true", dest="tree_support", help="Scale tree branch widths by support (if present)", default=False)
	group.add_option("-6", "--brlens", action="store_true", dest="show_branchlengths", help="Label branches with branchlengths", default=False)
	group.add_option("-M", "--midpoint", action="store_true", dest="midpoint", help="Midpoint root tree", default=False)
	group.add_option("-L", "--ladderise", action="store", choices=['right', 'left'], dest="ladderise", help="page size [default= %default]", type="choice", default=None)
	group.add_option("-z", "--names_as_shapes", action="store", choices=['circle', 'square', 'rectangle', 'auto'], dest="names_as_shapes", help="Use shapes rather than taxon names in tree (choose from circle) [default= %default]", type="choice", default="auto")
	group.add_option("-1", "--logbranches", action="store_true", dest="log_branches", help="page size [default= %default]", default=False)
	
	parser.add_option_group(group)
	
	
	group = OptionGroup(parser, "Output Options For Tracks Without Tree")
	
	group.add_option("-x", "--taxon_list", action="store", dest="taxon_list", help="File with ordered taxa", default="")
	group.add_option("-n", "--names", action="store_true", dest="track_names", help="Show track names", default=False)
	
	parser.add_option_group(group)
	
	
	group = OptionGroup(parser, "Names Options")
	
	group.add_option("-a", "--aligntaxa", action="store", dest="aligntaxa", help="Align taxon names (0=no, 1=align left, 2=align right) [default= %default]", default=0, type="int")
	group.add_option("-N", "--taxanames", action="store_false", dest="taxon_names", help="Do not show taxa names on tree. [default=show names]", default=True)
	group.add_option("-m", "--metadata", action="store", dest="metadata", help="metadata file in csv format. Note that this file must have a header row ontaining titles for each column", default="")
	group.add_option("-c", "--columns", action="store", dest="columns", help="column(s) from metadata file to use for track name (comma separated list) [default=%default]", default="1")
	group.add_option("-C", "--colourbycolumns", action="store", dest="colour_columns", help="column(s) from metadata file to use to colour track name and blocks next to name (comma separated list). If names are being shown, the first column will be used to colour names. All following columns will be added as coloured shapes as defined by the -z option. [default=%default]", default=False)
	group.add_option("-k", "--nometadatakey", action="store_false", dest="show_metadata_key", help="Do not show metadata keys. [default=show keys]", default=True)
	group.add_option("-r", "--parsimony_reconstruction", action="store", dest="transformation", help="Reconstruct colours across branches using parsimony. Select from acctran or deltran transformations [default=%default]", default=None, type="choice", choices=['acctran', 'deltran'])
	group.add_option("-i", "--suffix", action="store", dest="suffix", help="suffix to remove from filenames", default="")
	
	parser.add_option_group(group)
	
	
	group = OptionGroup(parser, "Annotation Track Options")
	
	group.add_option("-b", "--beginning", action="store", dest="beginning", help="Start position", default=0, metavar="int", type="int")
	group.add_option("-e", "--end", action="store", dest="end", help="End position", default=-1, metavar="int", type="int")
	group.add_option("-f", "--fragments", action="store", dest="fragments", help="number of fragments (lines on which to split genome) for linear diagram [default= %default]", type="int", default=1, metavar="int")
	group.add_option("-F", "--fragment_separation", action="store", dest="fragment_separation", help="Number of blank tracks between fragments  [default= %default]", default=1, metavar="int", type="int")
	group.add_option("-T", "--tracksize", action="store", dest="tracksize", help="Proportion of space available to each track that should be used in drawing", default=0.8, metavar="float", type="float")
	group.add_option("-E", "--emblheight", action="store", dest="emblheight", help="Relative track height of embl tracks to other tracks [default= %default]", default=2, metavar="int", type="int")
	group.add_option("-l", "--labels", action="store", dest="labels", help="Relative track height of track for CDS labels (0 = do not show labels) [default= %default]", default=0, metavar="int", type="int")
	group.add_option("-g", "--label_angle", action="store", dest="label_angle", help="Angle of CDS labels (anticlockwise from horizontal) [default= %default]", default=45, metavar="int", type="int")
	group.add_option("-B", "--outline_features", action="store_true", dest="outline_features", help="Outline all embl features in black", default=False)
	group.add_option("-G", "--greytracks", action="store_true", dest="greytracks", help="Make tracks behind features grey", default=False)
	group.add_option("-q", "--qualifier", action="store", dest="qualifier", help="qualifier in tab file containing strain names to sort by (must be the same as in the tree provided) [e.g. note]", default="")
	group.add_option("-A", "--arrowheads", action="store", choices=['0', '1', '2'], dest="arrows", help="add arrowheads to features on forward and reverse strands. 0=No arrows, 1=arrow style 1, 2=arrow style 2 [default= %default]", type="choice", default=0)
	group.add_option("-S", "--oneline", action="store_true", dest="oneline", help="Draw forward and reverse embl objects on one line", default=False)
	
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "Bam options")
	
	group.add_option("-u", "--base_qual_filter", action="store", dest="base_qual_filter", help="Base quality filter for bam file mapping plots [default= %default]", default=0, type="float")
	group.add_option("-U", "--mapping_qual_filter", action="store", dest="mapping_qual_filter", help="Mapping quality filter for bam file plots [default= %default]", default=0, type="float")
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "Bcf options")
	group.add_option("-v", "--varianttype", action="store", dest="bcfvariants", choices=["a","A","i","I","s","S", "h", "H"], help="bcf variants to show. Letter code for types of variant to include: a/A=include all i/I:just indels, s/S:just SNPs, h/H just heterozygotous sites. Lower case means exclude non-variant sites, upper case means include them [default= %default]", default="A", type="choice")
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "Plot Options")
	
	group.add_option("-H", "--plotheight", action="store", dest="plotheight", help="Relative track height of plot tracks to other tracks [default= %default]", default=2, metavar="int", type="int")
	group.add_option("-d", "--default_plot_type", choices=["hist", "heat", "bar", "line", "area", "stackedarea"], type="choice", action="store", dest="plottype", help="Set the default plot type (for plots called plot or graph and bam files). Choose from "+", ".join(["hist", "heat", "bar", "line", "area", "stackedarea"])+" [default= %default]", default="line")
	group.add_option("-X", "--scale_plots", action="store_true", dest="scale_plots_same", help="Use the same max and min scale for all plots", default=False)
	group.add_option("-y", "--plot_min", action="store", dest="plot_max", help="Set a maximum value for plots", default="Inf", type="float")
	group.add_option("-Y", "--plot_max", action="store", dest="plot_min", help="Set a minimum value for plots", default="Inf", type="float")
	group.add_option("-Z", "--log_plot", action="store_true", dest="log_plots", help="Show plots on a log scale (doesn't work yet)", default=False)
	group.add_option("-w", "--plot_scales", action="store_true", dest="plot_scales", help="Show legend on heatmaps", default=False)
	group.add_option("-3", "--heatmap_colour", default="bluered", choices=["redblue", "bluered", "blackwhite", "whiteblack"], type="choice", action="store", dest="heat_colour", help="Set the default plot type (for plots called plot or graph and bam files). Choose from "+", ".join(["redblue", "bluered", "blackwhite", "whiteblack"]))
	group.add_option("-4", "--windows", action="store", dest="windows", help="Number of windows per line (be careful, making this value too high will make things very slow and memory intensive) [default= %default]", default=501, type="float")
	
	parser.add_option_group(group)
	
	return parser.parse_args()



###############################################
# Check the command line options are sensible #
###############################################









#######################################################################
# Function to reconstruct a character across the tree using parsimony #
#######################################################################

def parsimony_reconstruction(treeObject, namecolours, colours, transformation="acctran"):
	
	
	def make_transistion_matrix():
		
		for name in namecolours:
			if not namecolours[name] in transition_matrix:
				transition_matrix[colours[namecolours[name]]]={}
		
		if not (0,0,0) in transition_matrix:
			transition_matrix[(0,0,0)]={}
		
		for colour in transition_matrix:
			for colourb in transition_matrix:
				if colour==colourb:
					transition_matrix[colour][colourb]=0
				else:
					transition_matrix[colour][colourb]=1
	
	
	def sankoff(treeObject, node=-1):
		
		if node==-1:
			node=treeObject.root
		
		daughters=treeObject.node(node).get_succ()
	
		if treeObject.is_terminal(node):
			
			node_data=treeObject.node(node).get_data()
			node_data.comment["sankoff"]={}
			
			if "name_colour" in node_data.comment:
				node_data.comment["branch_colour"]=node_data.comment["name_colour"]
			else:
				node_data.comment["branch_colour"]=(0,0,0)
				
			for state in transition_matrix.keys():
				node_data.comment["sankoff"][state]=100
				
			node_data.comment["sankoff"][node_data.comment["branch_colour"]]=0

				
			treeObject.node(node).set_data(node_data)
			
				
			
		else:
			node_data=treeObject.node(node).get_data()
			node_data.comment["sankoff"]={}
			for state in transition_matrix.keys():
				node_data.comment["sankoff"][state]=0
			for daughter in daughters:
				
				treeObject=sankoff(treeObject, daughter)
			 	
			 	daughter_sankoff=treeObject.node(daughter).get_data().comment["sankoff"]
			 	
			 	for state in transition_matrix.keys():
			 		min_cost=100000
					for comparison_state in daughter_sankoff.keys():
						cost=transition_matrix[state][comparison_state]+daughter_sankoff[comparison_state]
						if cost<min_cost:
							min_cost=cost
						
					node_data.comment["sankoff"][state]+=min_cost
				
			treeObject.node(node).set_data(node_data)
		
		return treeObject	
		
		
		
	def sankoff_second_traversal(treeObject, node=-1, transformation="acctran"):
		
		if node==-1:
			node=treeObject.root
			
		daughters=treeObject.node(node).get_succ()	
		
		if node==treeObject.root:
			min_cost_states=set()
		 	min_cost=100000
		 	
		 	
			for state in transition_matrix.keys():
				cost=treeObject.node(node).get_data().comment["sankoff"][state]
		 		if cost<min_cost:
					min_cost=cost
					min_cost_states=set()
					min_cost_states.add(state)
				elif cost==min_cost:
					min_cost_states.add(state)
			
			root_states=set()
			
			
			
			if len(min_cost_states)>1 and (0,0,0) in min_cost_states:
				min_cost_states.remove((0,0,0))
			
			if len(min_cost_states)>1:# and transformation=="deltran":
				
				parent_states_min_cost={}
			 	
			 	for state in min_cost_states:
			 		parent_states_min_cost[state]=[set()]
			 		daughter_min_cost_states=set()
			 		daughter_min_cost=100000
					for comparison_state in transition_matrix.keys():
						cost=0
						for daughter in daughters:
							cost+=transition_matrix[state][comparison_state]+treeObject.node(daughter).get_data().comment["sankoff"][comparison_state]
						if cost<daughter_min_cost:
							daughter_min_cost=cost
							daughter_min_cost_states=set()
							daughter_min_cost_states.add(comparison_state)
							parent_states_min_cost[state]=[set(),cost]
							parent_states_min_cost[state][0].add(comparison_state)
						elif cost==daughter_min_cost:
							daughter_min_cost_states.add(comparison_state)
							parent_states_min_cost[state][0].add(comparison_state)
				

				min_cost=100000
				for state in parent_states_min_cost:
					if parent_states_min_cost[state][1]<min_cost:
						min_cost=parent_states_min_cost[state][1]
						
				root_states=set()
				for state in parent_states_min_cost:
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
				min_cost=100000
				for state in possible_states:
					for comparison_state in transition_matrix.keys():
						cost=transition_matrix[state][comparison_state]+treeObject.node(daughters[0]).get_data().comment["sankoff"][comparison_state]
						if cost<min_cost:
							min_cost=cost
							min_cost_state=state
				possible_states=[min_cost_state]
				#print sitenumber, possible_states
				
			node_data.comment["branch_colour"]=possible_states[0]
				
			treeObject.node(node).set_data(node_data)
			
		node_data=treeObject.node(node).get_data()
		node_state=node_data.comment["branch_colour"]
		node_state_set=set(node_state)
		
		
		for daughter in daughters:
			
		 	daughter_data=treeObject.node(daughter).get_data()
		 	daughter_sankoff=daughter_data.comment["sankoff"]
		 	
		 	
	 		min_cost_states=set()
	 		min_cost=100000
			for comparison_state in daughter_sankoff.keys():
				cost=transition_matrix[node_state][comparison_state]+daughter_sankoff[comparison_state]
				if cost<min_cost:
					min_cost=cost
					min_cost_states=set()
					min_cost_states.add(comparison_state)
				elif cost==min_cost:
					min_cost_states.add(comparison_state)
			
			if len(min_cost_states)>1 and (0,0,0) in min_cost_states:
				min_cost_states.remove((0,0,0))
			
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
			
			
			daughter_data.comment["branch_colour"]=possible_states[0]
				
			treeObject.node(daughter).set_data(daughter_data)
				
			
			sankoff_second_traversal(treeObject, daughter, transformation=transformation)

	
		return treeObject
		
	
#	def print_sankoffs():
#		node=treeObject.root
#		
#		def print_node_sankoff(node):
#			for daughter in treeObject.node(node).get_succ():
#				print_node_sankoff(daughter)
#			if "branch_colour" in treeObject.node(node).get_data().comment:
#				print node, treeObject.node(node).get_data().comment["branch_colour"]
#			else:
#				print node, "unknown"
#		
#		print_node_sankoff(node)
	
		
	
	transition_matrix={}
	make_transistion_matrix()
	
	treeObject=sankoff(treeObject, treeObject.root)
	treeObject=sankoff_second_traversal(treeObject, treeObject.root, transformation=transformation)

			
	return treeObject




####################################################
# Function to round floats to n significant digits #
####################################################

def round_to_n(x, n):
    if n < 1:
        raise ValueError("number of significant digits must be >= 1")
    # Use %e format to get the n most significant digits, as a string.
    format = "%." + str(n-1) + "e"
    as_string = format % x
    if x>=10 or x<=-10:
    	return int(float(as_string))
    else:
	    return float(as_string)



##############################################################################################################
# Function to convert features with subfeatures (e.g. pseudogenes) to a list of locations of the subfeatures #
##############################################################################################################

def iterate_subfeatures(feature, locations):
	if len(feature.sub_features)>0:
		for subfeature in feature.sub_features:
			locations=iterate_subfeatures(subfeature, locations)
	else:
		locations.append((feature.location.start.position+1, feature.location.end.position+1))
	
	
	return locations	


####################################################
# Function to get the pixel width of a text string #
####################################################

def get_text_width(font, size, text):
	c = Canvas(test,pagesize=pagesize)
	length= c.stringWidth(str(text),font,size)
	
	return length

###############################################################################
# Function to find the best font size for a strint to fit in a certain length #
###############################################################################

def set_text_width(font, size, length, text):
	c = Canvas(test,pagesize=pagesize)
	strlen=float("Inf")
	while strlen>length:
		strlen= c.stringWidth(text,font,size)
		size-=0.1
	return size




##############################################
# Function to add a sequence file to a track #
##############################################

def add_sequence_file_to_diagram(fastarecord, name):
	pos=0
	odd=True
	new_track = Track()
	new_track.name=name
	for record in fastarecord:
		seqlen=len(str(record.seq))
		if odd:
			odd=False
			colour=colors.orange
		else:
			odd=True
			colour=colors.brown
		new_track.add_feature(locations=[(pos,pos+seqlen)], fillcolour=colour, strokecolour=colors.black, strokeweight=0, label=record.name, strand=0, arrows=0)
		pos+=seqlen

	return new_track



###########################################
# Function to add an embl file to a track #
###########################################

def add_embl_to_diagram(record, incfeatures=["CDS", "feature", "tRNA", "rRNA", "repeat_region"], emblfile=True, name=""):
	
	########################################
	# Function to get a name for a feature #
	########################################
	
	def get_best_feature_name(feature):
		
		name_types=["gene", "primary_name", "systematic_id", "locus_tag"]
		
		for name in name_types:
			if feature.qualifiers.has_key(name):
				return feature.qualifiers[name][0]
	
		return ""
	
	new_track = Track()
	if name=="":
		new_track.name=record.name
	else:
		new_track.name=name
	
	incfeatures=map(string.lower,incfeatures)
		
	print len(record.features), "features found for", record.name	
	
	if len(record.seq)>500000:
		scale_largetick_interval=int(round((len(record.seq)/10),-5))
		scale_smalltick_interval=int(round((len(record.seq)/10),-5)/5)
	else:
		scale_largetick_interval=len(record.seq)
		scale_smalltick_interval=len(record.seq)/5
	
		
	
#	if options.misc_features:
#		incfeatures.append("misc_feature")
	
	for x, feature in enumerate(record.features):
		if feature.type.lower() not in incfeatures or feature.location.nofuzzy_end<options.beginning or (feature.location.nofuzzy_start>options.end and options.end!=-1):
			#Exclude this feature
			#print "here"
			continue
		
			
		if feature.qualifiers.has_key("colour"):
			colourline=feature.qualifiers["colour"][0]
		elif feature.qualifiers.has_key("color"):
			colourline=feature.qualifiers["color"][0]
		else:
			if feature.type.lower()=="cds":
				colourline = "5"
			elif feature.type.lower()=="rrna":
				colourline = "0"
			elif feature.type.lower()=="trna":
				colourline = "8"
			elif feature.type.lower()=="repeat_region":
				colourline = "9"
			else:
				colourline = "1"

		
		if len(colourline.split())==1:
			try:
				colour=translator.artemis_color(colourline)
			except StandardError:
				colour=translator.artemis_color("5")
		elif len(colourline.split())==3:
			colour=translator.int255_color((int(colourline.split()[0]),int(colourline.split()[1]),int(colourline.split()[2])))
		else:
			print "Can't understand colour code!"
			colour=translator.artemis_color("5")
			
		locations=[]
		#get gene locations (including subfeatures)
		locations=iterate_subfeatures(feature, locations)
		
		if feature.type.lower()=="cds":
			new_track.add_feature(locations, fillcolour=colour, strokecolour=colour, strokeweight=0.5, strand=feature.strand, arrows=int(options.arrows), label=get_best_feature_name(feature))
			#gd_feature_set.add_feature(feature, color=colour, label=0, sigil=sigiltype, arrowhead_length=0.25, locations=locations)
		else:
			new_track.add_feature(locations, fillcolour=colour, strokecolour=colour, label=get_best_feature_name(feature), arrows=int(options.arrows))
			#gd_feature_set.add_feature(feature, color=colour, label=0, strand=0, locations=locations)


	return new_track


#########################################
# Function to add a tab file to a track #
#########################################

def add_tab_to_diagram(filename):
	
	features=[]
	
	if filename.split(".")[-1]=="gz":
		print "need to add gzip functionality to tab reader"
		return
	
	try:
		record=tab_parser(open(filename,"r"))
	except IOError:
		print "Cannot find file", filename
		sys.exit()
	record.name=filename
	track=add_embl_to_diagram(record, incfeatures=["i", "d", "li", "del", "snp", "misc_feature", "core", "CORE", "cds", "insertion", "deletion", "recombination", "feature", "blastn_hit", "fasta_record", "contig", "repeat_region"], emblfile=False)
	
	return track


#####################################################################################
# Function to add an embl file to multiple tracks split using the qualifiers option #
#####################################################################################

def add_ordered_embl_to_diagram(record, incfeatures=["CDS", "feature"], emblfile=True):

	incfeatures=map(string.lower,incfeatures)
	
	new_tracks={}
	
	
	print len(record.features), "features found for", record.name	
	
	if len(record.seq)>500000:
		scale_largetick_interval=int(round((len(record.seq)/10),-5))
		scale_smalltick_interval=int(round((len(record.seq)/10),-5)/5)
	else:
		scale_largetick_interval=len(record.seq)
		scale_smalltick_interval=len(record.seq)/5
	
		
	
#	if options.misc_features:
#		incfeatures.append("misc_feature")
	
	for x, feature in enumerate(record.features):
		if feature.type.lower() not in incfeatures or feature.location.nofuzzy_end<options.beginning or (feature.location.nofuzzy_start>options.end and options.end!=-1):# or options.qualifier not in feature.qualifiers:#,"tRNA","repeat_unit"] :
			#Exclude this feature
			#print "here"
			continue
		
		if options.beginning!=0:
		
			if feature.location.nofuzzy_start<options.beginning and feature.location.nofuzzy_end>options.beginning:
				feature.location=FeatureLocation(options.beginning, feature.location.nofuzzy_end)
		
		if options.end!=-1:
			if feature.location.nofuzzy_start<options.end and feature.location.nofuzzy_end>options.end:
				feature.location=FeatureLocation(feature.location.nofuzzy_start, options.end)
				
		if feature.qualifiers.has_key("colour"):
			colourline=feature.qualifiers["colour"][0]
		elif feature.qualifiers.has_key("color"):
			colourline=feature.qualifiers["color"][0]
		else:
			colourline = "5"
		if len(colourline.split())==1:
			colour=translator.artemis_color(colourline)
		elif len(colourline.split())==3:
			colour=translator.int255_color((int(colourline.split()[0]),int(colourline.split()[1]),int(colourline.split()[2])))
		else:
			print "Can't understand colour code!"
			print colourline
			sys.exit()
			
		locations=[]
		#get gene locations (including subfeatures)
		locations=iterate_subfeatures(feature, locations)
#		if feature.type=="CDS":	
#			gd_feature_set.add_feature(feature, color=colour, label=0, sigil=sigiltype, arrowhead_length=0.25, locations=locations)
#		else:
#			gd_feature_set.add_feature(feature, color=colour, label=0, strand=0, locations=locations)
		
		if options.qualifier in feature.qualifiers:
			qualifiernames=feature.qualifiers[options.qualifier][0].replace(", "," ").split()
			
			for taxonname in qualifiernames:
				taxonname=taxonname.strip()
				if not taxonname in new_tracks:
					newtrack = Track()
					newtrack.name=taxonname
					new_tracks[taxonname]=newtrack
					
				if feature.type.lower()=="cds":
					arrows=int(options.arrows)
				else:
					arrows=0
				new_tracks[taxonname].add_feature(locations, fillcolour=colour, strokecolour=colour, arrows=arrows)
		
		else:
			if not record.name in new_tracks:
				newtrack = Track()
				newtrack.name=record.name
				new_tracks[record.name]=newtrack	
			if feature.type.lower()=="cds":
				arrows=int(options.arrows)
			else:
				arrows=0
			new_tracks[record.name].add_feature(locations, fillcolour=colour, strokecolour=colour, arrows=arrows)
			
		
	if len(new_tracks)>1 and record.name in new_tracks:
		del new_tracks[record.name]
	return new_tracks


###################################################################################
# Function to add a tab file to multiple tracks split using the qualifiers option #
###################################################################################

def add_ordered_tab_to_diagram(filename):
	
	features={"":[]}
	
	featurename=""
	names_to_add_feature_to=[]
	
	if filename.split(".")[-1]=="gz":
		print "need to add gzip functionality to tab reader. Until then, please unzip your input files."
		return
	
	try:
		record=tab_parser(open(filename,"r"))
	except IOError:
		print "Cannot find file", filename
		sys.exit()
	record.name=filename
	new_tracks=add_ordered_embl_to_diagram(record, incfeatures=["i", "d", "li", "del", "snp", "misc_feature", "core", "cds", "insertion", "deletion", "recombination", "feature", "blastn_hit", "fasta_record", "repeat_region"], emblfile=False)
	return new_tracks



###########################################
# Function to add a bcf file to a diagram #
###########################################


def add_bcf_to_diagram(filename):		
	
	new_track = Track()
	
	features=[]
	bcftoolssarg = shlex.split(BCFTOOLS_DIR+"bcftools view "+filename)

	returnval = subprocess.Popen(bcftoolssarg, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

	stdout, stderr  = returnval.communicate()
	
	if len(stderr)>0:
		bcftoolssarg = shlex.split(OLD_BCFTOOLS_DIR+"bcftools view "+filename)
		returnval = subprocess.Popen(bcftoolssarg, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

		stdout, stderr  = returnval.communicate()

		if len(stderr)>0:
			print "Failed to open ", filename, "bcftools error:", stderr
			sys.exit()
	
	
	lines=stdout.split("\n")
	
	
	inmapped=False
	lastbase=0
	for line in lines:
	
		#add filters here?
		words=line.strip().split()
		if len(words)==0:
			continue
		if words[0][0]=="#":
			if words[0][1]!="#":
				headings=words
				headings[0]=headings[0][1:]
			continue
		
		if len(words)!=len(headings):
			print "words not equal to headings"
			print headings
			print words
			sys.exit()
		
		BASEINFO={}
		
		for x, heading in enumerate(headings):
#			if "ALT" in BASEINFO and BASEINFO["ALT"]==".":
#				break
			if heading=="INFO":
				
				BASEINFO[heading]={}
				
				try: info=words[x].split(";")
				except StandardError:
					print "Cannot split info string", words[x]
					sys.exit()
				for i in info:
					
					infotype=i.split("=")[0]
					
					if len(i.split("="))<2:
						if infotype=="INDEL":
							BASEINFO[heading][infotype]=True
					else:
						infodata=i.split("=")[1]
						try: BASEINFO[heading][infotype]=float(infodata)
						except StandardError:
							try: BASEINFO[heading][infotype]=map(float,infodata.split(","))
							except StandardError:
								BASEINFO[heading][infotype]=infodata
				
				
					
			else:
				try: BASEINFO[heading]=float(words[x])
				except StandardError:
					BASEINFO[heading]=words[x]
				
		
		
		#filter the call
		keep=True
		SNP=True
		INDEL=False
		if options.end>-1 and int(BASEINFO["POS"])>options.end:
			break
		if  int(BASEINFO["POS"])<options.beginning:
			continue
		if BASEINFO["ALT"]==".":
			SNP=False
			
#			if options.bcfvariants not in ["A", "I", "S"]:
#				continue
		
		
		#Calculate the ref/alt ratios
#		BASEINFO["INFO"]["DP4ratios"]={}
#		if not "DP4" in BASEINFO["INFO"]:
#			BASEINFO["INFO"]["DP4"]=[0,0,0,0]
#			BASEINFO["INFO"]["DP4ratios"]["fref"]=0.0
#			BASEINFO["INFO"]["DP4ratios"]["rref"]=0.0
#			BASEINFO["INFO"]["DP4ratios"]["falt"]=0.0
#			BASEINFO["INFO"]["DP4ratios"]["ralt"]=0.0
#			BASEINFO["INFO"]["AF1"]=0
#			BASEINFO["INFO"]["MQ"]=0
#		elif "DP4" in BASEINFO["INFO"]:
#			try: BASEINFO["INFO"]["DP4ratios"]["fref"]=float(BASEINFO["INFO"]["DP4"][0])/(BASEINFO["INFO"]["DP4"][0]+BASEINFO["INFO"]["DP4"][2])
#			except ZeroDivisionError:
#				BASEINFO["INFO"]["DP4ratios"]["fref"]=0.0
#			try: BASEINFO["INFO"]["DP4ratios"]["rref"]=float(BASEINFO["INFO"]["DP4"][1])/(BASEINFO["INFO"]["DP4"][1]+BASEINFO["INFO"]["DP4"][3])
#			except ZeroDivisionError:
#				BASEINFO["INFO"]["DP4ratios"]["rref"]=0.0
#			try: BASEINFO["INFO"]["DP4ratios"]["falt"]=float(BASEINFO["INFO"]["DP4"][2])/(BASEINFO["INFO"]["DP4"][0]+BASEINFO["INFO"]["DP4"][2])
#			except ZeroDivisionError:
#				BASEINFO["INFO"]["DP4ratios"]["falt"]=0.0
#			try: BASEINFO["INFO"]["DP4ratios"]["ralt"]=float(BASEINFO["INFO"]["DP4"][3])/(BASEINFO["INFO"]["DP4"][1]+BASEINFO["INFO"]["DP4"][3])
#			except ZeroDivisionError:
#				BASEINFO["INFO"]["DP4ratios"]["ralt"]=0.0
		
	
		
		
		
#		if BASEINFO["QUAL"]<options.QUAL:
#			#print options.QUAL, BASEINFO["QUAL"]
#			keep=False
#		elif  BASEINFO["INFO"]["MQ"]<options.MQUAL:
#			#print options.MQUAL, BASEINFO["INFO"]["MQ"]
#			keep=False
#		elif not SNP and BASEINFO["INFO"]["DP4"][0]+BASEINFO["INFO"]["DP4"][1]<options.depth:
#			keep=False
#		elif not SNP and BASEINFO["INFO"]["DP4"][0]<options.stranddepth:
#			keep=False
#		elif not SNP and BASEINFO["INFO"]["DP4"][1]<options.stranddepth:
#			keep=False
#		elif not SNP and BASEINFO["INFO"]["DP4ratios"]["fref"]<options.ratio:
#			keep=False
#		elif not SNP and BASEINFO["INFO"]["DP4ratios"]["rref"]<options.ratio:
#			keep=False
#		elif SNP and BASEINFO["INFO"]["DP4"][2]+BASEINFO["INFO"]["DP4"][3]<options.depth:
#			keep=False
#		elif SNP and BASEINFO["INFO"]["DP4"][2]<options.stranddepth:
#			keep=False
#		elif SNP and BASEINFO["INFO"]["DP4"][3]<options.stranddepth:
#			keep=False
#		elif SNP and BASEINFO["INFO"]["DP4ratios"]["falt"]<options.ratio:
#			keep=False
#		elif SNP and BASEINFO["INFO"]["DP4ratios"]["ralt"]<options.ratio:
#			keep=False
#		elif BASEINFO["ALT"]=="." and BASEINFO["INFO"]["AF1"]>(1-options.AF1):
#			keep=False
#		elif BASEINFO["ALT"]!="." and BASEINFO["INFO"]["AF1"]<options.AF1:
#			keep=False
#		elif SNP and "PV4" in BASEINFO["INFO"]:
#			if BASEINFO["INFO"]["PV4"][0]<=options.strand_bias:
#				keep=False
#			if BASEINFO["INFO"]["PV4"][1]<=options.baseq_bias:
#				keep=False
#			if BASEINFO["INFO"]["PV4"][2]<=options.mapping_bias:
#				keep=False
#			if BASEINFO["INFO"]["PV4"][3]<=options.tail_bias:
#				keep=False
			
		
		HETERO=False
		#find hetrozygous SNP calls and INDELS
		if len(BASEINFO["ALT"].split(","))>1:
			HETERO=True
#			keep=False
		elif (len(BASEINFO["ALT"].split(",")[0])>1 or len(BASEINFO["REF"].split(",")[0])>1) and "INDEL" in BASEINFO['INFO']:
			INDEL=True
		elif "INDEL" in BASEINFO['INFO']:
			keep=False
		
		if not "DP" in BASEINFO["INFO"] or int(BASEINFO["INFO"]["DP"])<1:
			keep=False
		
		
				
		if inmapped and lastbase+1!=int(BASEINFO["POS"]):
			inmapped=False
			new_track.add_feature([(start,end+1)], fillcolour=colour, strokecolour=colour)
		lastbase=int(BASEINFO["POS"])
		
		if keep:
			if not SNP and not INDEL:
				if inmapped==False:
					colour=colors.Color(230.0/255,230.0/255,230.0/255)
					inmapped=True
					start=int(BASEINFO["POS"])
				end=int(BASEINFO["POS"])
			elif SNP:
				if inmapped:
					inmapped=False
					new_track.add_feature([(start,end+1)], fillcolour=colour, strokecolour=colour)
					
#					try: features.append((feature,colour))
#					except NameError: pass #print "here"
				if INDEL and options.bcfvariants not in ["S", "H", "i"] and len(BASEINFO["ALT"])>len(BASEINFO["REF"]):
					colour=translator.artemis_color("6")
					start=int(BASEINFO["POS"])
					end=int(BASEINFO["POS"])
					#feature = SeqFeature(FeatureLocation(start, end), strand=None)
					new_track.add_feature([(start,end)], fillcolour=colour, strokecolour=colour)
					try: features.append((feature,colour))
					except NameError: pass #print "here"
				elif INDEL and options.bcfvariants not in ["S", "H", "i"] and len(BASEINFO["ALT"])<len(BASEINFO["REF"]):
					#colour='1'
					colour=translator.artemis_color("1")
					start=int(BASEINFO["POS"])
					end=int(BASEINFO["POS"])+(len(BASEINFO["REF"])-len(BASEINFO["ALT"]))
					#feature = SeqFeature(FeatureLocation(start, end), strand=None)
					new_track.add_feature([(start,end)], fillcolour=colour, strokecolour=colour)
#					try: features.append((feature,colour))
#					except NameError: pass #print "here"
				elif SNP and options.bcfvariants not in ["S", "I", "h"] and "AF1" in BASEINFO["INFO"] and BASEINFO["INFO"]["AF1"]<0.8 and BASEINFO["INFO"]["AF1"]>0.2:
#					colour='10'
					colour=translator.artemis_color("10")
					start=int(BASEINFO["POS"])
					end=int(BASEINFO["POS"])
					new_track.add_feature([(start,end)], fillcolour=colour, strokecolour=colour)
#					feature = SeqFeature(FeatureLocation(start, end), strand=None)
#					try: features.append((feature,colour))
#					except NameError: pass #print "here"
				elif SNP and options.bcfvariants not in ["s", "I", "H"] and BASEINFO["ALT"]=="A":
#					colour='3'
					colour=translator.artemis_color("3")
					start=int(BASEINFO["POS"])
					end=int(BASEINFO["POS"])
					new_track.add_feature([(start,end)], fillcolour=colour, strokecolour=colour)
#					feature = SeqFeature(FeatureLocation(start, end), strand=None)
#					try: features.append((feature,colour))
#					except NameError: pass #print "here"
				elif SNP and options.bcfvariants not in ["s", "I", "H"] and BASEINFO["ALT"]=="C":
					colour=translator.artemis_color("2")
					start=int(BASEINFO["POS"])
					end=int(BASEINFO["POS"])
					new_track.add_feature([(start,end)], fillcolour=colour, strokecolour=colour)
#					feature = SeqFeature(FeatureLocation(start, end), strand=None)
#					try: features.append((feature,colour))
#					except NameError: pass #print "here"
				elif SNP and options.bcfvariants not in ["s", "I", "H"] and BASEINFO["ALT"]=="G":
					colour=translator.artemis_color("4")
					start=int(BASEINFO["POS"])
					end=int(BASEINFO["POS"])
					new_track.add_feature([(start,end)], fillcolour=colour, strokecolour=colour)
#					feature = SeqFeature(FeatureLocation(start, end), strand=None)
#					try: features.append((feature,colour))
#					except NameError: pass #print "here"
				elif SNP and options.bcfvariants not in ["s", "I", "H"] and BASEINFO["ALT"]=="T":
					colour=translator.artemis_color("14")
					start=int(BASEINFO["POS"])
					end=int(BASEINFO["POS"])
					new_track.add_feature([(start,end)], fillcolour=colour, strokecolour=colour)
#					feature = SeqFeature(FeatureLocation(start, end), strand=None)
#					try: features.append((feature,colour))
#					except NameError: pass #print "here"
		elif inmapped:
#			feature = SeqFeature(FeatureLocation(start, end), strand=None)
			new_track.add_feature([(start,end+1)], fillcolour=colour, strokecolour=colour)
			inmapped=False
			end=int(BASEINFO["POS"])
#			try: features.append((feature,colour))
#			except NameError: pass #print "here"
	
	if inmapped:
#			feature = SeqFeature(FeatureLocation(start, end), strand=None)
		new_track.add_feature([(start,end+1)], fillcolour=colour, strokecolour=colour)
		inmapped=False
		end=int(BASEINFO["POS"])

	

	print len(new_track.features), "features found in", filename
	
	return new_track




#############################
# Function to draw the tree #
#############################

def drawtree(treeObject, treeheight, treewidth, xoffset, yoffset, name_offset=5):
	
	
	def get_max_branch_depth():
		
		terminals=treeObject.get_terminals()
		maxbrlen=0.0
		for terminal in terminals:
			if treeObject.sum_branchlength(node=terminal)>maxbrlen:
				maxbrlen=treeObject.sum_branchlength(node=terminal)
		
		return maxbrlen
	
	
	def count_downstream_nodes(node,count=0):
		def count_next_node(node, count):
			for daughter in treeObject.node(node).succ:
				count=count_next_node(daughter, count)
			count+=1
			return count
		count=count_next_node(node, count)
		return count
	
	
	
	def draw_scale():
		
		if vertical_scaling_factor<5:
			linewidth=0.5
		else:
			linewidth=1.0
		branchlength=round_to_n(max_branch_depth/10, 2)*horizontal_scaling_factor
		horizontalpos=xoffset+round_to_n(max_branch_depth/10, 2)*horizontal_scaling_factor
		vertpos=treebase-fontsize
		scalestring = str(round_to_n(max_branch_depth/10, 2))
		scalefontsize=fontsize
		if scalefontsize<6:
			scalefontsize=6
		d.add(Line(horizontalpos, vertpos, horizontalpos+branchlength, vertpos, strokeWidth=linewidth))
		d.add(String(horizontalpos+(float(branchlength)/2), vertpos-(scalefontsize+1), scalestring, textAnchor='middle', fontSize=scalefontsize, fontName='Helvetica'))
		
		
	def log_branchlengths():
		
		
		def min_branchlength(node, min_brlen):
			daughters=treeObject.node(node).succ
			brlen=treeObject.node(node).data.branchlength
			
			if brlen!=0 and brlen<min_brlen:
				min_brlen=brlen
			
			for daughter in daughters:
				daughterbrlen=min_branchlength(daughter, min_brlen)
				if daughterbrlen<min_brlen:
					min_brlen=daughterbrlen
			
			return min_brlen
		
		
		def log_branchlength(node, multiplier):
		
			daughters=treeObject.node(node).succ
			if treeObject.node(node).data.branchlength>0:
				treeObject.node(node).data.branchlength=log(treeObject.node(node).data.branchlength*multiplier)
			
			for daughter in daughters:
				log_branchlength(daughter, multiplier)
		
		node=treeObject.root
		min_br_len=min_branchlength(node, float("Inf"))
		
		x=1
		while x*min_br_len<=1:
			x=x*10
		
		print min_br_len
		node=treeObject.root
		log_branchlength(node, x)
	
	
	def get_outgroup_nodes_above_and_below(node):
		
		outgroup_nodes_above=0
		outgroup_nodes_below=0
		
		def get_outgroup_prevnode(node, outgroup_nodes_above, outgroup_nodes_below):
		
			if node==treeObject.root:
				return outgroup_nodes_above, outgroup_nodes_below
			
			parentnode=treeObject.node(node).prev
			outgroup_nodes_above, outgroup_nodes_below=get_outgroup_prevnode(parentnode, outgroup_nodes_above, outgroup_nodes_below)
			sisters=treeObject.node(parentnode).succ
			above=True
			for sister in sisters:
				if sister==node:
					above=False
				elif above:
					if treeObject.is_internal(node=sister):
						outgroup_nodes_above+=treeObject.count_terminals(node=sister)
					else:
						outgroup_nodes_above+=1
				else:
					if treeObject.is_internal(node=sister):
						outgroup_nodes_below+=treeObject.count_terminals(node=sister)
					else:
						outgroup_nodes_below+=1
				
	
			return outgroup_nodes_above, outgroup_nodes_below
		
		
		outgroup_nodes_above, outgroup_nodes_below = get_outgroup_prevnode(node, outgroup_nodes_above, outgroup_nodes_below)
		return outgroup_nodes_above, outgroup_nodes_below
	
	
	
	def get_node_vertical_positions():
		
		def get_node_vertical_position(node):
			
			for daughter in treeObject.node(node).succ:
				get_node_vertical_position(daughter)
			
#			if treeObject.is_terminal(node):
#				outgroup_nodes_above, outgroup_nodes_below=get_outgroup_nodes_above_and_below(node)
#				treeObject.node(node).data.comment={}
#				treeObject.node(node).data.comment["vertpos"]=(outgroup_nodes_below+yoffset)*vertical_scaling_factor
#				treeObject.node(node).data.comment["height"]=vertical_scaling_factor
			
			if not treeObject.is_terminal(node):
				daughters=treeObject.node(node).succ
#				treeObject.node(node).data.comment={}
				treeObject.node(node).data.comment["vertpos"]=float(treeObject.node(daughters[0]).data.comment["vertpos"]+treeObject.node(daughters[-1]).data.comment["vertpos"])/2
				
#				treeObject.node(node).data.comment["height"]=vertical_scaling_factor
				
		
		node=treeObject.root
		get_node_vertical_position(node)
		
	
	def drawbranch(node,horizontalpos):
		
		vertpos=treeObject.node(node).data.comment["vertpos"]+yoffset
		
		horizontalpos+=xoffset
		
		branchlength=treeObject.node(node).data.branchlength*horizontal_scaling_factor
		
		
		if options.tree_support:
			max_width=vertical_scaling_factor*0.8
			
			if treeObject.node(node).data.support==None:
				linewidth=max_width
			else:
				linewidth=(treeObject.node(node).data.support/100)*max_width
		else:
			if vertical_scaling_factor<5:
				linewidth=0.5
			else:
				linewidth=1.0
		
		# if branches have colours, find out now
		if treeObject.node(node).data.comment and "branch_colour" in treeObject.node(node).data.comment:
			r,g,b=treeObject.node(node).data.comment["branch_colour"]
			branch_colour=colors.Color(float(r)/255,float(g)/255,float(b)/255)
		else:
			branch_colour=colors.black
			
#		c.line(horizontalpos,vertpos,horizontalpos+branchlength,vertpos)
		if branchlength<linewidth:
			branchlength=linewidth
		d.add(Line(horizontalpos-(linewidth/2), vertpos, (horizontalpos-(linewidth/2))+branchlength, vertpos, strokeWidth=linewidth, strokeColor=branch_colour))
		
		if options.show_branchlengths and treeObject.node(node).data.branchlength>0:
			d.add(String((horizontalpos-(linewidth/2))+(branchlength/2), vertpos+linewidth, str(treeObject.node(node).data.branchlength).rstrip('0').rstrip('.'), textAnchor='middle', fontSize=fontsize*0.9, fillColor='black', fontName='Helvetica'))
		
		
		if node!=treeObject.root:
	
			parentnode=treeObject.node(node).prev
			sisters=treeObject.node(parentnode).succ
			parentvertpos=treeObject.node(parentnode).data.comment["vertpos"]+yoffset
			d.add(Line(horizontalpos, vertpos, horizontalpos, parentvertpos, strokeWidth=linewidth, strokeColor=branch_colour))
		
		
		if treeObject.is_terminal(node):
			
			if treeObject.node(node).data.comment and "name_colour" in treeObject.node(node).data.comment:
				name_colours=[]
				for x in xrange(0,len(treeObject.node(node).data.comment["name_colour"])):
					r,g,b= treeObject.node(node).data.comment["name_colour"][x]
					name_colours.append(colors.Color(float(r)/255,float(g)/255,float(b)/255))
			else:
				name_colours=[colors.black]
			
			
			# calculate total length of gubbins to add
			
			gubbins_length=0.0
			
			colpos=0
			if options.taxon_names:
				namewidth=get_text_width('Helvetica', fontsize, treeObject.node(node).data.taxon)+name_offset
				gubbins_length += namewidth
				colpos=1
			
			for x in xrange(colpos,len(name_colours)):
				gubbins_length += block_length
				if x!=0:
					gubbins_length += vertical_scaling_factor
			
			#Add the taxon names if present
			if options.taxon_names:
				if options.aligntaxa==2:
					d.add(String(treewidth+xoffset+(max_name_width-gubbins_length)+(fontsize/2), vertpos-(fontsize/3), treeObject.node(node).data.taxon, textAnchor='start', fontSize=fontsize, fillColor=name_colours[0], fontName='Helvetica'))
					block_xpos=treewidth+xoffset+(max_name_width-gubbins_length)+(fontsize/2)+namewidth
				elif options.aligntaxa==1:
					d.add(String(treewidth+(fontsize/2)+xoffset, vertpos-(fontsize/3), treeObject.node(node).data.taxon, textAnchor='start', fillColor=name_colours[0], fontSize=fontsize, fontName='Helvetica'))
					block_xpos=treewidth+(fontsize/2)+xoffset+(fontsize/2)+namewidth
				else:
					d.add(String(horizontalpos+branchlength+(fontsize/2), vertpos-(fontsize/3), treeObject.node(node).data.taxon, textAnchor='start', fontSize=fontsize, fillColor=name_colours[0], fontName='Helvetica'))
					block_xpos=horizontalpos+branchlength+(fontsize/2)+namewidth
			else:
				if options.aligntaxa==2:
					block_xpos=treewidth+xoffset+(max_name_width-gubbins_length)+(fontsize/2)
				elif options.aligntaxa==1:
					block_xpos=treewidth+(fontsize/2)+xoffset+(fontsize/2)
				else:
					block_xpos=horizontalpos+branchlength+(fontsize/2)
			
			
			
			# draw dashed lines
			
			if options.aligntaxa==1:
				d.add(Line(horizontalpos+branchlength, vertpos, treewidth+xoffset, vertpos, strokeDashArray=[1, 2], strokeWidth=linewidth/2, strokeColor=name_colours[0]))
				#d.add(Line(treewidth+xoffset+max_name_width, vertpos, horizontalpos+branchlength+gubbins_length, vertpos, strokeDashArray=[1, 2], strokeWidth=linewidth/2, strokeColor=name_colours[0]))
			elif options.aligntaxa==2:
				#d.add(Line(horizontalpos+branchlength, vertpos, treewidth+xoffset, vertpos, strokeDashArray=[1, 2], strokeWidth=linewidth/2, strokeColor=name_colours[0]))
				d.add(Line(horizontalpos+branchlength, vertpos, treewidth+xoffset+(max_name_width-gubbins_length), vertpos, strokeDashArray=[1, 2], strokeWidth=linewidth/2, strokeColor=name_colours[0]))
				#d.add(Line(horizontalpos+branchlength, vertpos, treewidth+xoffset, vertpos, strokeWidth=linewidth/2, strokeColor=name_colours[0]))
			
			for x, name_colour in enumerate(name_colours[colpos:]):
#				if options.aligntaxa==1:
#					if options.names_as_shapes=="circle":
#						d.add(Circle(cx=horizontalpos+branchlength+vertical_scaling_factor, cy=vertpos, r=(vertical_scaling_factor/2), fillColor=name_colour, strokeColor=name_colour, strokeWidth=0))
#					elif options.names_as_shapes=="square":
#						d.add(Rect(horizontalpos+branchlength+(vertical_scaling_factor/2), vertpos-(vertical_scaling_factor/2), vertical_scaling_factor, vertical_scaling_factor, fillColor=name_colour, strokeColor=name_colour, strokeWidth=0))
#					elif options.names_as_shapes=="rectangle":
#						d.add(Rect(horizontalpos+branchlength+(vertical_scaling_factor/2), vertpos-(vertical_scaling_factor/2), vertical_scaling_factor*2, vertical_scaling_factor, fillColor=name_colour, strokeColor=name_colour, strokeWidth=0))
#					elif options.names_as_shapes=="auto":
#						d.add(Rect(horizontalpos+branchlength+(vertical_scaling_factor/2), vertpos-(vertical_scaling_factor/2), max_name_width, vertical_scaling_factor, fillColor=name_colour, strokeColor=name_colour, strokeWidth=0))
#					elif options.taxon_names:
#						d.add(Line(horizontalpos+branchlength, vertpos, treewidth+xoffset, vertpos, strokeDashArray=[1, 2], strokeWidth=linewidth/2, strokeColor=name_colour))
#						d.add(String(treewidth+(fontsize/2)+xoffset, vertpos-(fontsize/3), treeObject.node(node).data.taxon, textAnchor='start', fillColor=name_colour, fontSize=fontsize, fontName='Helvetica'))
#					else:
#						d.add(Line(horizontalpos+branchlength, vertpos, treewidth+xoffset, vertpos, strokeDashArray=[1, 2], strokeWidth=linewidth/2, strokeColor=name_colour))
#				
#				elif options.aligntaxa==2:





				if options.names_as_shapes=="circle":
					d.add(Circle(cx=block_xpos+block_length, cy=vertpos, r=(block_length/2), fillColor=name_colour, strokeColor=None, stroke=False, strokeWidth=0))
				elif options.names_as_shapes in ["square", "rectangle", "auto"]:
					d.add(Rect(block_xpos, vertpos-(vertical_scaling_factor/2), block_length, vertical_scaling_factor, fillColor=name_colour, strokeColor=None, stroke=False, strokeWidth=0))
						
						
						
						
						
#					elif options.taxon_names:
#						namewidth=get_text_width("Helvetica", fontsize, treeObject.node(node).data.taxon)+name_offset
#						#d.add(Line(horizontalpos+branchlength, vertpos, treewidth+xoffset+(max_name_width-namewidth), vertpos, strokeDashArray=[1, 2], strokeWidth=linewidth/2, strokeColor=name_colour))
#						d.add(String(treewidth+xoffset+(max_name_width-namewidth)+(fontsize/2), vertpos-(fontsize/3), treeObject.node(node).data.taxon, textAnchor='start', fontSize=fontsize, fillColor=name_colour, fontName='Helvetica'))
#					else:
#						d.add(Line(horizontalpos+branchlength, vertpos, treewidth+xoffset, vertpos, strokeDashArray=[1, 2], strokeWidth=linewidth/2, strokeColor=name_colour))
				
#				else:
#					if options.names_as_shapes=="circle":
#						d.add(Circle(cx=horizontalpos+branchlength+vertical_scaling_factor, cy=vertpos, r=(vertical_scaling_factor/2), fillColor=name_colour, strokeColor=name_colour, strokeWidth=0))
#					elif options.taxon_names:
#						d.add(String(horizontalpos+branchlength+(fontsize/2), vertpos-(fontsize/3), treeObject.node(node).data.taxon, textAnchor='start', fontSize=fontsize, fillColor=name_colour, fontName='Helvetica'))
				
				block_xpos+=block_length+vertical_scaling_factor
			
			
		
	def recurse_subtree(node, horizontalpos, treebase=float("Inf"), treetop=float("-Inf")):
		
		daughters=treeObject.node(node).succ
		
		daughterhorizontalpos=horizontalpos+(treeObject.node(node).data.branchlength*horizontal_scaling_factor)
		drawbranch(node,horizontalpos)
		for daughter in daughters:
			treebase, treetop=recurse_subtree(daughter,daughterhorizontalpos, treebase=treebase, treetop=treetop)
		if treeObject.node(node).data.comment["vertpos"]<treebase:
			treebase=treeObject.node(node).data.comment["vertpos"]
		if treeObject.node(node).data.comment["vertpos"]>treetop:
			treetop=treeObject.node(node).data.comment["vertpos"]
		
		return treebase, treetop
		
		
	
	def get_max_name_width(name_offset, fontsize):
		max_width=0.0
		for taxon in treeObject.get_terminals():
			curwidth= get_text_width("Helvetica", fontsize, treeObject.node(taxon).data.taxon)
			if curwidth>max_width:
				max_width=curwidth
		
		return max_width
	
	
	def draw_column_label(xpox, ypos, fontsize, column_label):
		
		mylabel=Label()
		
		mylabel.setText(column_label)
		
		mylabel.angle=45
		mylabel.fontName="Helvetica"
		mylabel.fontSize=fontsize
		mylabel.x=xpox
		mylabel.y=ypos
		mylabel.boxAnchor="sw"
		#Axis.setPosition(self.track_position[0], self.track_position[1]+(self.track_height/2), self.track_length)
		d.add(mylabel)
	
	
	
#	vertical_scaling_factor=float(treeheight)/(treeObject.count_terminals(node=treeObject.root)+2)
	fontsize=vertical_scaling_factor
	if fontsize>12:
		fontsize=12
	name_offset=fontsize
	if options.taxon_names:
		while get_max_name_width(name_offset, fontsize)+name_offset>treewidth/3:
			fontsize-=0.2
		max_name_width=get_max_name_width(name_offset, fontsize)+name_offset
		colblockstart=1
	else:
		max_name_width=0
		colblockstart=0
	
	block_length=0
	if len(colour_dict)>0:
		max_total_name_length=(float(treewidth)/4)*3
		
		max_total_block_length=(max_total_name_length-max_name_width)
		
		
		max_block_length=((max_total_block_length-(vertical_scaling_factor*((len(colour_dict)-1)+colblockstart)))/len(colour_dict))	
		
		if max_block_length<vertical_scaling_factor:
			print ("Not enough space to draw your metadata colour columns")
			sys.exit()
		
		if options.names_as_shapes in ["circle", "square"]:
			block_length=vertical_scaling_factor
		elif options.names_as_shapes=="rectangle":
			if max_block_length>(vertical_scaling_factor*2):
				block_length=vertical_scaling_factor*2
			else:
				block_length=max_block_length
		elif options.names_as_shapes=="auto":
			if (treewidth/20)<max_block_length and (treewidth/20)<20:
				block_length=(treewidth/20)
			elif max_block_length>20:
				block_length=20
			else:
				block_length=max_block_length
		
	gubbins_length=0.0
	for x in range(colblockstart,len(colour_dict)):
		if x>0:
			max_name_width+=vertical_scaling_factor
			gubbins_length+=vertical_scaling_factor
		
		max_name_width+=block_length
		gubbins_length+=block_length
		
	
	treewidth-=(max_name_width+(fontsize/2)+5)
	
#	treewidth-=xoffset
	
	if options.log_branches:
		log_branchlengths()
	
	max_branch_depth=get_max_branch_depth()
	horizontal_scaling_factor=float(treewidth)/max_branch_depth
	
	get_node_vertical_positions()
	
	treebase, treetop=recurse_subtree(treeObject.root, 0)
	treebase+=yoffset
	treetop+=yoffset
	
	#treebase=treeObject.node(treeObject.get_terminals()[-1]).data.comment["vertpos"]+yoffset
	
	if not options.show_branchlengths:
		draw_scale()
	
	if fontsize>6:
		labelfontsize=fontsize
	else:
		labelfontsize=6
	
	
	try:
		if colour_column_names:
			if options.aligntaxa==2:
				column_name_x_pos=treewidth+xoffset+(max_name_width-gubbins_length)+(fontsize/2)
				column_name_y_pos=treetop+(vertical_scaling_factor/2)
				colpos=0
				if options.taxon_names:
					
					draw_column_label(treewidth+xoffset+((max_name_width-gubbins_length)/2), column_name_y_pos, labelfontsize, colour_column_names[0])
				
					colpos=1
				for x in xrange(colpos,len(colour_column_names)):
					
					draw_column_label(column_name_x_pos+(block_length/2), column_name_y_pos, labelfontsize, colour_column_names[x])
					column_name_x_pos += block_length
					column_name_x_pos += vertical_scaling_factor
	except NameError:
		pass
			
		
	
	return



#################
# Drawing class #
#################

class Figure:
	def __init__(self, beginning, end):
		self.begnining=0
		self.end=-1




###############
# Track class #
###############


class Track:
	def __init__(self, track_position=[-1,-1], track_height=0, track_length=0, track_draw_proportion=0.75, scale=False, tick_marks=True, tick_mark_number=5, tick_mark_labels=True, minor_tick_marks=True, minor_tick_mark_number=3, features=[], beginning=0, end=-1):
	
		self.track_position=track_position#horizontal and vertical position of centre of track
		self.track_height=track_height#height of space allocated for track
		self.track_length=track_length
		self.track_draw_proportion=track_draw_proportion#proportion of the track that should be used for drawing
		self.scale=scale
		self.scale_position="middle"
		self.tick_marks=tick_marks
		self.tick_mark_number=tick_mark_number
		self.tick_mark_labels=tick_mark_labels
		self.tick_mark_label_font="Helvetica"
		self.tick_mark_label_size=8
		self.tick_mark_label_angle=45
		self.minor_tick_marks=minor_tick_marks
		self.minor_tick_mark_number=minor_tick_mark_number
		self.feature_label_track_height=0
		self.features=features[:]
		self.scaled_features=features[:]
		self.draw_feature_labels=False
		self.feature_label_size=8
		self.feature_label_angle=0
		self.feature_label_font="Helvetica"
		self.greytrack=False
		self.grey_track_colour=colors.Color(0.25,0.25,0.25)
		self.grey_track_opacity_percent=10
		self.max_feature_length=-1
		self.beginning=0
		self.end=-1
		self.track_number=-1
		self.plots=[]
		self.fragments=1
		self.name=""
		self.show_name=False
		self.name_font="Helvetica"
		self.name_size=12
		self.name_length=0
		self.is_key=False
		self.key_data=[]

	
	def draw_greytrack(self):
		self.grey_track_colour.alpha=float(self.grey_track_opacity_percent)/100
		d.add(Rect(self.track_position[0], self.track_position[1]-((float(self.track_height)/2)*self.track_draw_proportion), self.track_length, float(self.track_height)*self.track_draw_proportion, fillColor=self.grey_track_colour, strokeColor=self.grey_track_colour, strokeWidth = 0))
	
	
	def add_plot(self, filename, plot_type="line", fragments=1):
		
		newplot=Plot()
		
		newplot.number_of_windows=newplot.number_of_windows*(fragments*options.npages)
		newplot.plot_type=plot_type
		if plot_type=="stackedarea":
			newplot.transparency=1.0
		datalines=newplot.read_plot_from_file(filename)
		
		newplot.beginning=self.beginning
		if self.end==-1:
			if len(newplot.raw_data)>0:
				newplot.end=len(newplot.raw_data[0])
		else:
			newplot.end=self.end

		
		newplot.read_data(datalines)
		
		if plot_type=="heat":
			newplot.heat_colour=options.heat_colour
		
		self.plots.append(newplot)
		
	
	
	def add_bam_plot(self, filename, plot_type="line", fragments=1):
		
		print "Calculating coverage for", filename
		sys.stdout.flush()
#		samtoolssarg = shlex.split(SAMTOOLS_DIR+"samtools view -H "+filename)
#		returnval = subprocess.Popen(samtoolssarg, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#	
#		stdout, stderr  = returnval.communicate()
#	
#		if len(stderr)>0:
#			print "Failed to open ", filename, "samtools error:", stderr
#			return
#
#		
#		headerlines=stdout.split("\n")
		
		try:
			samfile = pysam.Samfile( filename, "rb" )
		except StandardError:
			print 'Failed to open '+filename+'. Is it in bam format?'
			return
		
		refs=samfile.references
		lengths=samfile.lengths
		
		contiglocs={}
		totallength=0
		
		for x, ref in enumerate(refs):
			contiglocs[ref]=totallength
			totallength+=lengths[x]
		
#		for line in headerlines:
#			words=line.split()
#			if len(words)==3 and words[0]=="@SQ" and len(words[1].split(":"))==2 and words[1].split(":")[0]=="SN" and len(words[2].split(":"))==2 and words[2].split(":")[0]=="LN":
#				contiglocs[words[1].split(":")[1]]=totallength
#				totallength+=int(words[2].split(":")[1])
				
		
		
		newplot=Plot()
		newplot.number_of_windows=newplot.number_of_windows*fragments*(options.npages)
		newplot.plot_type=plot_type
		if plot_type=="heat":
			newplot.heat_colour=options.heat_colour
		
		poscount=0
		if options.base_qual_filter!=0 or options.mapping_qual_filter!=0:
			depths=[["#BASE", "High_Quality_Coverage", "Low_Quality_Coverage"]]
		else:
			depths=[["#BASE", "Coverage"]]	
		for x, ref in enumerate(refs):
			print ref
			lastcolumn=-1
			zerocount=0
			for pileupcolumn in samfile.pileup(ref):
				
				while pileupcolumn.pos!=lastcolumn+1:
					#print lastcolumn, pileupcolumn.pos
					poscount+=1
					if options.base_qual_filter!=0 or options.mapping_qual_filter!=0:
						depths.append([str(poscount), "0", "0"])
					else:
						depths.append([str(poscount), "0"])
					lastcolumn+=1
					zerocount+=1
				poscount+=1
				if options.base_qual_filter!=0 or options.mapping_qual_filter!=0:
					#print str(pileupcolumn)
					filtered_depth=0
					for pileupread in pileupcolumn.pileups:
						#print pileupread.alignment
						q=ord(pileupread.alignment.qual[pileupread.qpos])-33
						Q=pileupread.alignment.mapq
						#print q, Q, options.base_qual_filter, options.mapping_qual_filter
						if q>=options.base_qual_filter and Q>=options.mapping_qual_filter:
							filtered_depth+=1
					depths.append([str(poscount), str(filtered_depth), str(pileupcolumn.n-filtered_depth)])
					#sys.exit()
				else:
					depths.append([str(poscount), str(pileupcolumn.n)])
				lastcolumn=pileupcolumn.pos
			
			while lastcolumn+1<lengths[x]:
				poscount+=1
				if options.base_qual_filter!=0 or options.mapping_qual_filter!=0:
					depths.append([str(poscount), "0", "0"])
				else:
					depths.append([str(poscount), "0"])
				lastcolumn+=1
				zerocount+=1
			#print len(depths)
		
#		samtoolssarg = shlex.split(SAMTOOLS_DIR+"samtools depth -q "+str(options.base_qual_filter)+" -Q "+str(options.mapping_qual_filter)+" "+filename)
#		returnval = subprocess.Popen(samtoolssarg, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#	
#		stdout, stderr  = returnval.communicate()
#	
#		if len(stderr)>0:
#			print "Failed to open ", filename, "samtools error:", stderr
#			return
#		
#		
		newplot.beginning=self.beginning
		if self.end==-1:
			if len(newplot.raw_data)>0:
				newplot.end=len(newplot.raw_data[0])
		else:
			newplot.end=self.end
#
#		
#		datalines=stdout.split("\n")
		depth_data= map(" ".join,depths)
		
		#newplot.read_data(depths, samtools=True, contiglocs=contiglocs, totallength=totallength)
		newplot.read_data(depth_data)
		
#		newplot.raw_data_to_data()

		self.plots.append(newplot)
	
	
	
	
	def draw_name(self):
		#d.add(String(self.track_position[0]-(self.name_length+vertical_scaling_factor), self.track_position[1]-(self.name_size/3), self.name, textAnchor='start', fontSize=self.name_size, fontName='Helvetica'))
		
		d.add(String(self.track_position[0]-(self.name_length), self.track_position[1]-(self.name_size/3), self.name, textAnchor='start', fontSize=self.name_size, fontName='Helvetica'))
	
	
	def draw_scale(self):
		
		if self.end!=-1:
			length=float(self.end-self.beginning)
		else:
			length=float(self.max_feature_length-self.beginning)
		
		if self.scale_position=="middle":
			yposition=self.track_position[1]
		elif self.scale_position=="top":
			yposition=self.track_position[1]+((float(self.track_height)/2)*self.track_draw_proportion)
		elif self.scale_position=="bottom":
			yposition=self.track_position[1]-((float(self.track_height)/2)*self.track_draw_proportion)
		else:
			yposition=self.track_position[1]
			
		xAxis = XValueAxis() 
		xAxis.setPosition(self.track_position[0], yposition, self.track_length)
		xAxis.configure([(self.beginning,self.end)])
		xAxis.valueMax=self.end
		xAxis.valueMin=self.beginning
		if self.tick_marks:
			xAxis.tickUp=(float(self.track_height)/4)*self.track_draw_proportion
		else:
			xAxis.visibleTicks=0
		if self.minor_tick_marks:
			xAxis.visibleSubTicks=1
		if not self.tick_mark_labels:
			xAxis.visibleLabels=0
		xAxis.tickDown=0
		xAxis.labels.boxAnchor = 'w'
		xAxis.labels.dy=((float(self.track_height)/4)*self.track_draw_proportion)+(inch*0.02)
		xAxis.labels.angle=self.tick_mark_label_angle
		xAxis.labels.fontName=self.tick_mark_label_font
		xAxis.labels.fontSize=self.tick_mark_label_size

		d.add(xAxis) 

				
	def rescale_scale(self, maxvalue):
		if self.scale:
			max_scale_label_length=get_text_width(self.tick_mark_label_font, self.tick_mark_label_size,str(maxvalue))
			max_scale_label_height=sqrt(pow(max_scale_label_length,2)/2)
			height_including_label=(float(self.track_height*3)/4)+(max_scale_label_height)+(inch*0.02)
			if height_including_label<=self.track_height:
				return

			rescaled_label_height=(self.track_height/height_including_label)*max_scale_label_height
			rescaled_label_height=self.track_height-((float(self.track_height*3)/4)+(inch*0.02))
			rescaled_label_length=sqrt(pow(rescaled_label_height,2)*2)
			self.tick_mark_label_size=set_text_width(self.tick_mark_label_font, self.tick_mark_label_size, rescaled_label_length, str(maxvalue))

	
	
	def get_max_feature_length(self):
		max_feature_length=0
		for feature in self.features:
			for location in feature.feature_locations:
				if location[0]>max_feature_length:
					max_feature_length=location[0]
				if location[1]>max_feature_length:
					max_feature_length=location[1]
		return max_feature_length
	
	
	def scale_feature_positions(self):
		
		self.scaled_features=[]
		
		if self.end!=-1:
			length=float(self.end-self.beginning)
		else:
			length=float(self.max_feature_length-self.beginning)
		
		for feature in self.features:
		
			newfeature=Feature()
			newfeature.fillcolour=feature.fillcolour
			newfeature.strokecolour=feature.strokecolour
			newfeature.strokeweight=feature.strokeweight
			newfeature.strand=feature.strand
			newfeature.label=feature.label
			newfeature.arrows=feature.arrows
			scaledlocations=[]
			for location in feature.feature_locations:
				start=location[0]
				finish=location[1]
				if self.beginning!=0:
					if start<self.beginning and finish>self.beginning:
						start=self.beginning
				if self.end!=-1:
					if start<self.end and finish>self.end:
						finish=self.end
				start-=self.beginning
				finish-=self.beginning
				
				scaledlocations.append(((float(start)/length)*self.track_length,(float(finish)/length)*self.track_length))
			
			newfeature.feature_locations=scaledlocations
			self.scaled_features.append(newfeature)
			
	def draw_arrow(self,x, y, length, height, fillcolour, strokecolour, strokeweight, arrows, strand=0, headsize=-1):
		
		
		if headsize==-1:
			#headsize=float(self.track_length)/100
			headsize=float(height)/2
		if headsize>length:
			headsize=length
		
		if arrows==0:
			d.add(Rect(x, y, length, height, strokeColor = strokecolour, strokeWidth = strokeweight, fillColor = fillcolour))
		if arrows==1:
			if strand==1:
				d.add(Polygon(strokeColor = strokecolour, strokeWidth = strokeweight, fillColor = fillcolour, points=[x,y,x+length-headsize,y,x+length,y+(float(height)/2),x+length-headsize,y+height,x,y+height]))
			elif strand==-1:
				d.add(Polygon(strokeColor = strokecolour, strokeWidth = strokeweight, fillColor = fillcolour, points=[x+headsize,y,x+length,y,x+length,y+height,x+headsize,y+height,x,y+(float(height)/2)]))
			else:
				d.add(Rect(x, y, length, height, strokeColor = strokecolour, strokeWidth = strokeweight, fillColor = fillcolour))
				
		elif arrows==2:
			if strand==1:
				d.add(Polygon(strokeColor = strokecolour, strokeWidth = strokeweight, fillColor = fillcolour, points=[x,y+(float(height)/4),x+length-headsize,y+(float(height)/4),x+length-headsize,y,x+length,y+(float(height)/2),x+length-headsize,y+height,x+length-headsize,y+(float(height*3)/4),x,y+(float(height*3)/4)]))
			elif strand==-1:
				d.add(Polygon(strokeColor = strokecolour, strokeWidth = strokeweight, fillColor = fillcolour, points=[x+headsize,y,x+headsize,y+(float(height)/4),x+length,y+(float(height)/4),x+length,y+(float(height*3)/4),x+headsize,y+(float(height*3)/4),x+headsize,y+height,x,y+(float(height)/2)]))
			else:
				d.add(Rect(x, y+(float(height)/4), length, height/2, strokeColor = strokecolour, strokeWidth = strokeweight, fillColor = fillcolour))
		
		
	
	def draw_features(self):
			
		if self.max_feature_length==-1:
			return
		
		else:
			self.scale_feature_positions()
		
		featuresort=[]
		for x, feature in enumerate(self.scaled_features):
			featuresort.append([feature.feature_locations[0][0], x])
		#featuresort.sort()
		joins=[]
		
		
		
		for featurenum in featuresort[::-1]:
			feature=self.scaled_features[featurenum[1]]
			#if the feature is white, outline it in black so we can see it and outline features in black if selected in the options
			if feature.strokecolour==colors.Color(1,1,1,1) or options.outline_features:
				feature.strokecolour=colors.Color(0,0,0,1)
			else:
				feature.strokecolour=None
				feature.strokeweight=0
			
			
			
			subfeaturesort=[]
			for x, subfeature in enumerate(feature.feature_locations):
				subfeaturesort.append([subfeature[0], x])
			subfeaturesort.sort()
			subfeature_locations=[]
			for subfeaturenum in subfeaturesort:
				subfeature_locations.append(feature.feature_locations[subfeaturenum[1]])
			
#			print subfeaturesort
			for x, location in enumerate(subfeature_locations):
				
				
				if (location[0]>0 and location[0]<=self.track_length) or (location[1]>0 and location[1]<=self.track_length):
					
					if options.oneline:# and (feature.strand!=0 or feature.arrows!=2):
						y=self.track_position[1]-((float(self.track_height)/2)*self.track_draw_proportion)
						height=float(self.track_height)*self.track_draw_proportion
						y1=self.track_position[1]
						y2=self.track_position[1]+((float(self.track_height)/4)*self.track_draw_proportion)
					elif feature.strand==0:
						y=self.track_position[1]-((float(self.track_height)/4)*self.track_draw_proportion)
						height=(float(self.track_height)*self.track_draw_proportion)/2
						y1=self.track_position[1]
						y2=self.track_position[1]+((float(self.track_height)/8)*self.track_draw_proportion)
					elif feature.strand==1:
						y=self.track_position[1]
						height=(float(self.track_height)*self.track_draw_proportion)/2
						y1=self.track_position[1]+(height/2)
						y2=self.track_position[1]+((3*height)/4)
					elif feature.strand==-1:
						y=self.track_position[1]-((float(self.track_height)/2)*self.track_draw_proportion)
						height=(float(self.track_height)*self.track_draw_proportion)/2
						y1=self.track_position[1]-(height/2)
						y2=self.track_position[1]-((3*height)/4)
					else:
						y=self.track_position[1]-((float(self.track_height)/2)*self.track_draw_proportion)
						height=float(self.track_height)*self.track_draw_proportion
						y1=self.track_position[1]
						y2=self.track_position[1]+((float(self.track_height)/4)*self.track_draw_proportion)
					
					#print location[1], location[0],
					#seta  minimum width for a feature
					min_feature_width=0.5
					if location[1]-location[0]<min_feature_width:
						location=(location[0],location[0]+min_feature_width)
					#print location[1], location[0]
					
					if feature.arrows==0:
						d.add(Rect(self.track_position[0]+location[0], y, location[1]-location[0], height, fillColor=feature.fillcolour, strokeColor=feature.strokecolour, strokeWidth=feature.strokeweight))
					else:
						self.draw_arrow(self.track_position[0]+location[0], y, location[1]-location[0], height, feature.fillcolour, feature.strokecolour, feature.strokeweight, feature.arrows, strand=feature.strand)
						
					if len(subfeature_locations)>x+1 and subfeature_locations[x+1][0]<=self.track_length:
						if subfeature_locations[x+1][0]<location[1]:
							joinheight=y1
						elif y2>y1:
							if (y2-y1)>(float(subfeature_locations[x+1][0]-location[1])/2):
								joinheight=y1+(float(subfeature_locations[x+1][0]-location[1])/2)
							else:
								joinheight=y2
						else:
							if (y1-y2)>(float(subfeature_locations[x+1][0]-location[1])/2):
								joinheight=y1-(float(subfeature_locations[x+1][0]-location[1])/2)
							else:
								joinheight=y2
						
						joins.append(Line(self.track_position[0]+location[1], y1, self.track_position[0]+location[1]+(float(subfeature_locations[x+1][0]-location[1])/2), joinheight, strokeDashArray=[0.5, 1], strokeWidth=0.5))
						joins.append(Line(self.track_position[0]+((location[1]+subfeature_locations[x+1][0])/2), joinheight, self.track_position[0]+location[1]+(float(subfeature_locations[x+1][0]-location[1])), y1, strokeDashArray=[0.5, 1], strokeWidth=0.5))
						
				
			if self.draw_feature_labels:
				
				fstart=feature.feature_locations[0][0]
				fend=feature.feature_locations[-1][1]
				if fstart<self.track_length and fend>0:
					if fstart<0:
						fstart=0
					if fend>self.track_length:
						fend=self.track_length
					fmid=(float(fend-fstart)/2)+fstart
					
					mylabel=Label()
					mylabel.setText(feature.label)
					
					mylabel.angle=self.feature_label_angle
					while mylabel.angle>=360:
						mylabel.angle-=360
					while mylabel.angle<0:
						mylabel.angle+=360
					if mylabel.angle==0:
						mylabel.boxAnchor='s'
					elif mylabel.angle==90:
						mylabel.boxAnchor='w'
					elif mylabel.angle==180:
						mylabel.boxAnchor='n'
					elif mylabel.angle==270:
						mylabel.boxAnchor='e'
					elif (mylabel.angle>0 and mylabel.angle<90):
						mylabel.boxAnchor='sw'
					elif (mylabel.angle>90 and mylabel.angle<180):
						mylabel.boxAnchor='nw'
					elif (mylabel.angle>180 and mylabel.angle<270):
						mylabel.boxAnchor='ne'
					elif (mylabel.angle>270 and mylabel.angle<=360):
						mylabel.boxAnchor='se'
					mylabel.fontName=self.feature_label_font
					mylabel.fontSize=self.feature_label_size
					mylabel.x=self.track_position[0]+fmid
					mylabel.y=self.track_position[1]+(self.track_height/2)
					#Axis.setPosition(self.track_position[0], self.track_position[1]+(self.track_height/2), self.track_length)
					d.add(mylabel)
				
				
#		if self.draw_feature_labels:
#			labelAxis.valueSteps.sort()
#			print labelAxis.valueSteps
#			d.add(labelAxis)
			
			
		for join in joins:
			d.add(join)
			
		self.scaled_features=[]
			
	
	def draw_key(self):
		
		if (self.track_height/3)<5:
			linewidth=0.5
		else:
			linewidth=1.0
		
		
		def get_total_key_text_length(fontsize):
			totallength=0.0
			count=0
			for datum in self.key_data:
				totallength+=get_text_width("Helvetica", fontsize, str(datum[0]))
				count+=1
			totallength+=(count-1)*(fontsize)
			return totallength+(self.track_height*0.2)
		
		fontsize=self.track_height*0.8
		
		if fontsize>10:
			fontsize=10
		while get_total_key_text_length(fontsize)>self.track_length:
#			print fontsize, get_total_key_text_length(fontsize), self.track_length
			fontsize-=0.1
		
		position=self.track_position[0]+(self.track_height*0.1)
		
		for datum in self.key_data:
			d.add(String(position, self.track_position[1]-(float(fontsize)/3), str(datum[0]), textAnchor='start', fontSize=fontsize, fillColor=datum[1], fontName='Helvetica'))
			position+=get_text_width("Helvetica", fontsize, str(datum[0]))
			position+=(fontsize)
		position-=(fontsize)
		fillcolour=colors.Color(0,0,0,0)
		d.add(Rect(self.track_position[0], (self.track_position[1]-((float(fontsize)/2)))-(self.track_height*0.1), (position+(self.track_height*0.1))-self.track_position[0], float(fontsize)+(self.track_height*0.2), fillColor=colors.Color(0,0,0,0), strokeColor=colors.Color(0,0,0,1), strokeWidth = linewidth))
		
		
		
	
		
	def draw_track(self):
		
		
			
		if self.greytrack and not self.is_key and len(self.plots)==0:
			self.draw_greytrack()
		
		if self.is_key:
			self.draw_key()
			return
		
		self.draw_features()
		#self.scale=False
		if self.scale:
			self.draw_scale()
		
		if self.show_name:
			self.draw_name()
		
		for plot in self.plots:
			plot.beginning=self.beginning
			plot.end=self.end
			
#			plot.plot_position=self.track_position
#			plot.plot_height=self.track_height
#			plot.plot_length=self.track_length
#			plot.plot_draw_proportion=self.track_draw_proportion
			plot.draw_plot(self.track_position[0], self.track_position[1]-((float(self.track_height)/2)*self.track_draw_proportion), self.track_height*self.track_draw_proportion, self.track_length)


			
	def add_feature(self,locations=[(-1,-1)], fillcolour=colors.white, strokecolour=colors.black, strokeweight=0, label="", strand=0, arrows=0):
		
		newfeature=Feature()
		
		feature_locations=[]
		for location in locations:
			if location[0]>location[1]:
				feature_locations.append((location[1],location[0]))
			else:
				feature_locations.append((location[0],location[1]))
		
		newfeature.feature_locations=feature_locations
		newfeature.fillcolour=fillcolour
		newfeature.strokecolour=strokecolour
		newfeature.strokeweight=strokeweight
		newfeature.strand=strand
		newfeature.label=label
		newfeature.arrows=arrows
		
		
		self.features.append(newfeature)
		
		
	
	def sort_features_by_length(self):
		featurelist=[]
		ordered_features=[]
		for x, feature in enumerate(self.features):
			featurelist.append([feature.feature_locations[-1][1]-feature.feature_locations[0][0], x])
		featurelist.sort()
		#featurelist.reverse()
		
		for feature in featurelist:
			ordered_features.append(self.features[feature[1]])
		self.features=ordered_features[:]
		



#################
# Feature class #
#################


class Feature:
	def __init__(self):
	
		self.feature_locations=[(-1,-1)]
		self.strand=0
		self.arrows=0
		self.label=""
		self.fillcolour=colors.blue
		self.strokecolour=colors.black
		self.strokeweight=0


##############
# Plot class #
##############


class Plot:
	def __init__(self, scale=False, tick_marks=True, tick_mark_number=5, tick_mark_labels=True, beginning=0, end=-1):
	
		self.scale=scale
		self.tick_marks=tick_marks
		self.tick_mark_number=tick_mark_number
		self.tick_mark_labels=tick_mark_labels
		self.tick_mark_label_size=8
		self.tick_mark_label_angle=45
		self.label=""
		self.primary_colour=colors.red
		self.alternate_colour=colors.blue
		#self.heat_colour="whiteblack"
		self.heat_colour="bluered"
		self.line_colours=[colors.red, colors.blue, colors.green, colors.black, colors.magenta, colors.cyan, colors.yellow]
		self.strokeweight=0.5
		self.raw_data=[]
		self.data=[]
		self.xdata=[]
		self.max_yaxis=float("-Inf")
		self.min_yaxis=float("Inf")
		self.window_size=1
		self.number_of_windows=options.windows
		self.plot_type="line"
		self.beginning=-1
		self.end=-1
		self.circular=True
		self.labels=[]
		self.legend=True
		self.autolegend=True
		self.legend_font="Helvetica"
		self.legend_font_size=8
		self.reorder_data=True
		self.transparency=0.80
		self.data_order=[]
		self.max_feature_length=0
		
		
		
	def calculate_windowsize(self, data):
		self.window_size=0
		maxdatalength=max(map(len,data))
		
		if self.end!=-1 and self.beginning!=-1:
			self.window_size=int(float((self.end-self.beginning))/self.number_of_windows)
		elif self.end==-1 and self.beginning!=-1:
			self.window_size=int(float(maxdatalength-self.beginning)/self.number_of_windows)
		else:
			self.window_size=int(float(maxdatalength)/self.number_of_windows)
		
		if self.window_size<1:
			self.window_size=1
			
	
	def read_plot_from_file(self, filename):
	
		try: 
			datalines=open(filename)
		except StandardError:
			DoError("Cannot open file "+filename)
		
		return datalines
	
	
	def read_data(self, datalines, samtools=False,  contiglocs=[], totallength=-1):	
		
		data=[[]]
		if int(self.end)==-1:
			endpos=float("Inf")
		else:
			endpos=self.end
			self.circular=False
		
		if endpos<totallength:
			totallength=endpos
		elif totallength>0 and endpos>totallength:	
			endpos=totallength
		
		newplot=False
		currpos=0
		
		for line in datalines:
			if len(line.strip())<1:
				continue
			if line.strip()[0]=='#' and not newplot:
				newplot=True
				
				words=line.strip().split()
				for label in words[1:]:
					self.labels.append(label)
				for x in words[2:]:
					data.append([])
				continue
			elif line.strip()[0]=='#':
				continue
			else:
				currpos+=1
				if currpos>endpos:
					break
			if newplot:
				words=line.strip().split()
				while float(words[0])>currpos and float(words[0])<=endpos:# and float(words[0])>=self.beginning:
					for x in xrange(len(data)):
						data[x].append(0.0)
					currpos+=1
				if float(words[0])>=endpos:
					break
				sumtot=0.0
				for x in xrange(len(data)):
					data[x].append(float(words[x+1])+sumtot)
					if self.plot_type=="stackedarea":
						sumtot+=float(words[x+1])
			elif samtools:
				words=line.strip().split()
				while (float(words[1])+contiglocs[words[0]])>=currpos and (float(words[1])+contiglocs[words[0]])<=endpos:# and (float(words[1])+contiglocs[words[0]])>=self.beginning:
					data[0].append(0.0)
					currpos+=1
				if (float(words[1])+contiglocs[words[0]])>endpos:
					break
				data[0].append(float(words[2]))
			else:
				words=line.strip().split()
				if float(words[0])>=endpos:
					break
				data[0].append(float(words[0]))
		
		while currpos<totallength:
			for x in xrange(len(data)):
				data[x].append(0)
			currpos+=1
		
		if self.circular:
			for datum in data:
				if len(datum)>0:
					datum.insert(0,datum[-1])
				else:
					datum.append(0)
		
		self.calculate_windowsize(data)
		
		
		if len(data)>1 and self.reorder_data:# and self.plot_type!="stackedarea" :
			data_order=[]
			datameans=[]
			for x, datum in enumerate(data):
				datameans.append([mean(datum),x])
			datameans.sort()
			for datamean in datameans[::-1]:
				data_order.append(datamean[1])
			#remember to reorder the labels as well
			labels=self.labels[:]
			colours=self.line_colours[:]
			for x, y in enumerate(data_order):
				self.labels[x]=labels[y]
				colours[x].alpha=self.transparency
				self.line_colours[x]=colours[x]
				
		else:
			data_order=xrange(len(data))
		
		self.data=[]
		self.xdata=[]
		
		for x in data_order:
			windowdata=[]
			xdata=[]
#			mymax=0
			
			if self.end>len(data[x]):
				end=len(data[x])
			else:
				end=self.end
#			self.window_size=3
#			print self.end
			if self.plot_type in ["line", "area", "stackedarea"]:
				for y in xrange(0,len(data[x]),self.window_size):
					
					windowstart=int(y-floor((self.window_size-1)/2))
					windowend=int(y+floor((self.window_size-1)/2))+1
					
					if windowstart>=0 and windowend<len(data[x]):
						region=data[x][windowstart:windowend]
					elif windowend>=len(data[x]):
						region=data[x][windowstart:]+data[x][:windowend-len(data[x])]
					elif windowstart<0:
						region=data[x][len(data[x])-windowstart:]+data[x][:windowend]
					
					#print windowstart, windowend, y, len(data[x]), data[x][windowstart:windowend], region
					windowdata.append(mean(region))
					xdata.append(y)
#					
#				windowdata.insert(0,(0, mean(region)))
					
			elif self.plot_type in ["bar", "heat"]:
				bardata=[]
				barxdata=[]
				for y in xrange(0, len(data[x]),self.window_size):
					windowstart=int(y-floor((self.window_size-1)/2))
					windowend=int(y+floor((self.window_size-1)/2))+1
					
					if windowstart>=0 and windowend<len(data[x]):
						region=data[x][windowstart:windowend]
					elif windowend>=len(data[x]):
						region=data[x][windowstart:]+data[x][:windowend-len(data[x])]
					elif windowstart<0:
						region=data[x][len(data[x])-windowstart:]+data[x][:windowend]
					
#					print windowstart, windowend
					bardata.append(mean(region))
					barxdata.append(y)
#				bardata.insert(0,mean(region))
#				print bardata
				for y in xrange(len(bardata)-1):
					if y+1<len(bardata):
						windowdata.append(mean([bardata[y], bardata[y+1]]))
						xdata.append(barxdata[y+1])
					else:
						windowdata.append(mean([bardata[y], bardata[0]]))
			
			self.data.append(windowdata)
			self.xdata.append(xdata)
	
	
	
	def clear_data(self):
		self.raw_data=[]



	def get_data_to_print(self):
		
		
		if options.plot_min!=float("Inf"):
			valueMin = options.plot_min
		elif self.plot_type=="bar":
			valueMin = min(self.data[0])
		else:
			valueMin = min(map(min,self.data))
		if options.plot_max!=float("Inf"):
			valueMax = options.plot_max 
		elif self.plot_type=="bar":
			valueMax = max(self.data[0])
		else:
			valueMax = max(map(max,self.data))
		
		printdata=[]
		if self.plot_type in ["line", "area", "stackedarea"]:
			for x, data in enumerate(self.data):
				currdata=[]
				for y, datum in enumerate(data):
					#print datum, self.xdata[x][y], self.beginning, self.end
					if self.xdata[x][y]>=self.beginning and self.xdata[x][y]<=self.end:
						if datum>self.max_yaxis:
							datum=self.max_yaxis
						elif datum<self.min_yaxis:
							datum=self.min_yaxis
						currdata.append((self.xdata[x][y],datum))
				#print currdata
				printdata.append(currdata)
		elif self.plot_type in ["bar", "heat"]:
			for x, data in enumerate(self.data):
				currdata=[]
				for y, datum in enumerate(data):
					if self.xdata[x][y]>=self.beginning and self.xdata[x][y]<=self.end:
						if datum>self.max_yaxis:
							datum=self.max_yaxis
						elif datum<self.min_yaxis:
							datum=self.min_yaxis
						currdata.append(datum)
						end=self.xdata[x][y]
				printdata.append(currdata)
				
		if self.plot_type in ["line", "area", "stackedarea"]:
			return printdata, valueMin, valueMax
		elif self.plot_type in ["bar", "heat"]:
			return printdata, end, valueMin, valueMax
			
	
	
	def draw_heatmap(self, x, y, height, length):
	
		data, end, valueMin, valueMax=self.get_data_to_print()
				
		datalength=length*(float(end-self.beginning)/(self.end-self.beginning))
		
		feature_width=float(datalength)/len(data[0])
		
		
		
#		self.legend=True
#		if self.legend and len(self.labels)>0:
		if self.legend and options.plot_scales:
			
			if self.autolegend:
				self.legend_font_size=float(height)/6
				if self.legend_font_size>10:
					self.legend_font_size=10
			
			if height>self.legend_font_size:
				height-=self.legend_font_size
				y+=self.legend_font_size

#			legend = Legend()
			my_legend=ShadedRect()
			
#			sys.exit()
			my_legend.fillColorStart=colors.Color(0,0,1)
			my_legend.fillColorEnd=colors.Color(1,0,0)
			my_legend.height=self.legend_font_size
			
			
			my_legend.y=y-(self.legend_font_size+(self.legend_font_size/2))
			my_legend.strokeWidth=self.strokeweight
			my_legend.numShades=20
			my_legend.orientation='vertical'
			
			minlabel="%.2f" % valueMin
			maxlabel="%.2f" % valueMax
			minlabelwidth=get_text_width(self.legend_font, self.legend_font_size, minlabel)
			maxlabelwidth=get_text_width(self.legend_font, self.legend_font_size, maxlabel)
			legendwidth=float(length-(minlabelwidth+maxlabelwidth))/5
			if legendwidth<0:
				legendwidth=10
			my_legend.width=legendwidth
			my_legend.x=x+4+minlabelwidth
			

			d.add(String(x, my_legend.y, minlabel, textAnchor='start', fontSize=self.legend_font_size, fontName=self.legend_font))
			d.add(String(x+8+minlabelwidth+my_legend.width, my_legend.y, maxlabel, textAnchor='start', fontSize=self.legend_font_size, fontName=self.legend_font))
			d.add(my_legend)

		lastvalue=data[0][0]
		draw_width=0.0
		draw_start=0
		for i,datum in enumerate(data[0]):
			if valueMax-valueMin>0:
				value=round(float(datum-valueMin)/(valueMax-valueMin),2)
			else:
				value=0.0
			if value>1:
				value=1.0
			if value!=lastvalue:
				if self.heat_colour=="blackwhite":
					colour=colors.Color(lastvalue,lastvalue,lastvalue)
				elif self.heat_colour=="whiteblack":
					colour=colors.Color(1.0-lastvalue,1.0-lastvalue,1.0-lastvalue)
				elif self.heat_colour=="bluered":
					colour=colors.Color(lastvalue,0,1.0-lastvalue)
				elif self.heat_colour=="redblue":
					colour=colors.Color(1.0-lastvalue,0,lastvalue)
				if not (self.heat_colour=="whiteblack" and lastvalue==0) and not (self.heat_colour=="blackwhite" and lastvalue==1):
					d.add(Rect(x+(draw_start*feature_width), y, draw_width, height, fillColor=colour, strokeColor=None, stroke=False, strokeWidth=0))

#				myrect=Rect(x+(draw_start*feature_width), y, draw_width+1, height, fillColor=colour, strokeColor=None, stroke=False, strokeWidth=0)
				
#				print dir(myrect)
#				sys.exit()
				draw_width=0.0
				draw_start=i
				lastvalue=value
			
			draw_width+=feature_width
			
		
		if self.heat_colour=="blackwhite":
			colour=colors.Color(value,value,value)
		elif self.heat_colour=="whiteblack":
			colour=colors.Color(1.0-value,1.0-value,1.0-value)
		elif self.heat_colour=="bluered":
			colour=colors.Color(value,0,1.0-value)
		elif self.heat_colour=="redblue":
			colour=colors.Color(1.0-value,0,value)
		if not (self.heat_colour=="whiteblack" and value==0) and not (self.heat_colour=="blackwhite" and value==1):
			d.add(Rect(x+(draw_start*feature_width), y, draw_width, height, fillColor=colour, strokeColor=None, stroke=False, strokeWidth=0))
		
			
#		if (((i+1)*feature_width)+x)<self.max_feature_length:
#			d.add(Rect(x+((i+1)*feature_width), y, self.max_feature_length-(((i+1)*feature_width)+x), height, fillColor=colors.Color(0,0,1), strokeColor=colour, strokeWidth=0))
	
	
	def draw_line_plot(self, x, y, height, length):
		
		data, valueMin, valueMax=self.get_data_to_print()
		if height<30:
			self.legend=False
			
		
		if self.legend and len(self.labels)>0:
			
			if self.autolegend:
				self.legend_font_size=float(height)/6
				if self.legend_font_size>10:
					self.legend_font_size=10
			
			if height>self.legend_font_size:
				height-=self.legend_font_size
				y+=self.legend_font_size
				
			
			legend = LineLegend()
			
			legend.x = x
			if self.plot_type=="line":
				legend.y = y-(self.legend_font_size)
			elif self.plot_type in ["area", "stackedarea"]:
				legend.y = y-(self.legend_font_size/2)
			legend.strokeWidth=self.strokeweight
			legend.alignment='right'
			legend.fontName=self.legend_font
			legend.fontSize=self.legend_font_size
			legend.boxAnchor="nw"
			if self.plot_type=="stackedarea":
				legend.colorNamePairs  = []
				for i in xrange(len(data)):
					legend.colorNamePairs.append((self.line_colours[i], self.labels[len(data)-(i+1)]))
					#lp.lines[i].strokeColor=self.line_colours[len(lp.data)-(j+1)]
			else:
				legend.colorNamePairs  = [(self.line_colours[i], self.labels[i]) for i in xrange(len(data))]
			
			maxlabelwidth=0
			for label in self.labels:
				labelwidth=get_text_width(self.legend_font, self.legend_font_size, label)
				if labelwidth>maxlabelwidth:
					maxlabelwidth=labelwidth
			
			legend.columnMaximum=1  
			legend.yGap=0
			
			if self.plot_type=="line":
				legend.dy=self.strokeweight
				legend.dx=10
			elif self.plot_type in ["area", "stackedarea"]:
				legend.dy=self.legend_font_size
				legend.dx=self.legend_font_size
			legend.swdy=legend.dy/2
#			if self.legend_font_size>=10:
#				
			
				
			legend.dxTextSpace=4
			legend.deltax=maxlabelwidth+4  
			d.add(legend)
		
		
		lp = LinePlot()
		lp.x = x
		lp.y = y 
		lp.height = height 
		lp.width = length 
		lp.data = data
		
		lp.joinedLines = 1
		lp.xValueAxis.visibleLabels=0
		if not options.plot_scales:
			lp.yValueAxis.visibleLabels=0
		lp.xValueAxis.visibleTicks=0
		lp.xValueAxis.valueMin = self.beginning 
		if self.end!=-1:
			lp.xValueAxis.valueMax = self.end
		
		if options.plot_min!=float("Inf"):
			lp.yValueAxis.valueMin = options.plot_min
		else:
			lp.yValueAxis.valueMin = float("Inf")
			for datum in data:
				datum_min= min(map(lambda x: x[1], datum))
				if datum_min<lp.yValueAxis.valueMin:
					lp.yValueAxis.valueMin=datum_min

		if options.plot_max!=float("Inf"):
			lp.yValueAxis.valueMax = options.plot_max
		else:
			lp.yValueAxis.valueMax = float("-Inf")
			for datum in data:
				datum_max= max(map(lambda x: x[1], datum))
				if datum_max>lp.yValueAxis.valueMax:
					lp.yValueAxis.valueMax=datum_max
				
			
		
		lp.yValueAxis.tickLeft=0
		lp.yValueAxis.tickRight=5
		lp.yValueAxis.labels.fontSize=height/2
		if lp.yValueAxis.labels.fontSize<1:
			lp.yValueAxis.visibleLabels=0
		elif lp.yValueAxis.labels.fontSize>12:
			lp.yValueAxis.labels.fontSize=12
#		print dir(lp.yValueAxis.scale)
#		print lp.yValueAxis.scale(0)
		
#		lp.yValueAxis.valueSteps=[lp.yValueAxis.valueMin, lp.yValueAxis.valueMax]
		if self.plot_type=="stackedarea":
			for i in xrange(len(self.line_colours)):
				self.line_colours[i].alpha=1.0
		for i in xrange(len(lp.data)):
			j=i
			while j>=len(self.line_colours):
				j-=len(self.line_colours)
			if self.plot_type=="stackedarea":
				if j>len(self.line_colours):
					lp.lines[i].strokeColor=self.line_colours[len(self.line_colours)-(j+1)]
				else:
					lp.lines[i].strokeColor=self.line_colours[len(lp.data)-(j+1)]
			else:
				lp.lines[i].strokeColor=self.line_colours[j]
			lp.lines[i].strokeWidth=self.strokeweight
		if self.plot_type in ["area", "stackedarea"]:
			for i in xrange(len(lp.data)):
				lp.data[i].append((lp.data[i][-1][0], lp.yValueAxis.valueMin))
			lp._inFill=True
		d.add(lp)
	
	
	
	def draw_bar_plot(self, x, y, height, length):
		
		data, end, valueMin, valueMax=self.get_data_to_print()
		
		if self.legend and len(self.labels)>0:
			
			if self.autolegend:
				self.legend_font_size=float(height)/6
				if self.legend_font_size>10:
					self.legend_font_size=10
			
			if height>self.legend_font_size:
				height-=self.legend_font_size
				y+=self.legend_font_size

			legend = Legend()
			
			legend.x = x
			legend.y = y-(float(self.legend_font_size)/2)
			legend.strokeWidth=self.strokeweight
			legend.alignment='right'
			legend.fontName=self.legend_font
			legend.fontSize=self.legend_font_size
			legend.boxAnchor="nw"
			legend.colorNamePairs  = [(self.line_colours[i], self.labels[i]) for i in xrange(len(data))]
			
			maxlabelwidth=0
			for label in self.labels:
				labelwidth=get_text_width(self.legend_font, self.legend_font_size, label)
				if labelwidth>maxlabelwidth:
					maxlabelwidth=labelwidth
			
			legend.columnMaximum=1  
			legend.yGap=0
			
			legend.dy=self.legend_font_size
			legend.dx=self.legend_font_size
#			legend.swdy=legend.dy/2
		
			legend.dxTextSpace=4
			legend.deltax=maxlabelwidth+4  
			d.add(legend)		
		
		bc = VerticalBarChart()
		bc.x = x
		bc.y = y 
		bc.height = height 
		bc.width = length
		bc.data = data
		
		datalength=length*(float(end-self.beginning)/(self.end-self.beginning))
		
		bc.categoryAxis.visibleLabels=0
		bc.categoryAxis.visibleTicks=0
		if not options.plot_scales:
			bc.valueAxis.visibleLabels=0
		bc.categoryAxis.style = 'stacked'
		bc.groupSpacing = 0
		bc.barSpacing = 0

		
		bc.barWidth=datalength/len(bc.data[0])
		
		bc.useAbsolute=1
		if options.plot_min!=float("Inf"):
			bc.valueAxis.valueMin = options.plot_min
#		if self.min_yaxis!=float("Inf"):
#			bc.valueAxis.valueMin = round_to_n(self.min_yaxis, 2)
		else:
			bc.valueAxis.valueMin = round_to_n(min(map(min,data)), 2)
		if options.plot_max!=float("Inf"):
			bc.valueAxis.valueMax = options.plot_max
#		elif len(self.data)==1:
#			bc.valueAxis.valueMax = round_to_n(self.max_yaxis, 2)
		else:
			maxsum=0.0
			for datum in data:
				maxsum+=max(datum)
			bc.valueAxis.valueMax = round_to_n(maxsum, 2)
			
		
		bc.valueAxis.tickLeft=0
		bc.valueAxis.tickRight=5
		
#		print bc.valueAxis.valueMax, max(map(max,self.data))
#		bc.valueAxis.valueSteps=[bc.valueAxis.valueMin, bc.valueAxis.valueMax]
		
		
		for i in xrange(len(bc.data)):
			j=i
			while j>=len(self.line_colours):
				j-=len(self.line_colours)
			bc.bars[i].fillColor=self.line_colours[j]
			bc.bars[i].strokeWidth=self.strokeweight
		
		

		
		d.add(bc)

		

	
	def draw_plot(self, x, y, height, length):
#		self.raw_data_to_data(self)
		if self.plot_type in ["line", "area", "stackedarea"]:
			self.draw_line_plot(x, y, height, length)
		elif self.plot_type=="heat":
			self.draw_heatmap(x, y, height, length)
		elif self.plot_type in ["bar", "stack"]:
			self.draw_bar_plot(x, y, height, length)
		



def get_tree_colour_comments(treeobject):
	
	def get_node_colour_object(node):
	
		daughters=treeobject.node(node).succ
		
		for daughter in daughters:
			get_node_colour_object(daughter)
	
		if tree.node(node).data.comment!=None:
			comment=tree.node(node).data.comment.replace("[&", "").replace("]", "").split('=')
		else:
			comment=[]
		tree.node(node).data.comment={}
		if len(comment)>0 and comment[0] in ["color", "colour"]:
			rgb=comment[1].split()
			if len(rgb)==3:
				try:
					tree.node(node).data.comment["branch_colour"]=colors.Color(r, g, b)
				except StandardError:
					return
			

	get_node_colour_object(treeobject.root)
	return tree


################
# Main program #
################		

if __name__ == "__main__":

	#starttime=time.clock()

	#Get command line arguments

	(options, args) = main()
	
	if options.fragment_separation<0:
		options.fragment_separation=0
		options.fragment_separation=int(options.fragment_separation)
	
	if options.page=="custom":
		try:
			xy=map(float,options.custompage.split(","))
			if len(xy)!=2:
				print "Custom page size option format must be x,y"
				sys.exit()
			pagesize=(float(xy[0]), float(xy[1]))
		except StandardError:
			print "Invalid custom page size option (-5)"
			sys.exit()
	else:
		pagesize=pagesizeconverter[options.page]
	
#	print pagesize
#	sys.exit()
	
	#options.plottype="line"
	
	if options.orientation=="landscape":
		height, width = pagesize
	else:
		width, height = pagesize
	
	
	if len(args)==0 and options.tree!="":
		options.treeproportion=1
	elif len(args)==0:
		print "Found nothing to draw"
		sys.exit()
	
	
	
	
	margin=0.5*inch
	topmargin=margin
		
	
	
	
	metadatanames={}
	namecolours={}
	colour_dict=[]
	

		
	
	
	my_tracks={}
	
	
	#create translator object for translating artemis colours to GenomeDiagram colours
	
	translator = ColorTranslator()
	
	
	track_count=0
	
	tree_positions=[]

	track_names={}

	input_order=[]
	
	metadata_keylist=[]
	colour_columns=[0]
	max_length_col_name=0
	
	#Parse metadata file and create track for key
	
	if options.metadata!="":
		if options.columns:
			try: 
				columns=[]
				cols=options.columns.split(",")
				for col in cols:
					if len(col.split(".."))==1:
						try:
							columns.append(int(col))
						except StandardError:
							print "Could not understand columns selected:", options.columns
							columns=[1]
					elif len(col.split(".."))==2:
						colbits=col.split("..")
						for x in xrange(int(colbits[0]),int(colbits[1])+1):
							columns.append(x)
					else:
						print "Could not understand columns selected:", options.columns
						columns=[1]
#			
#				colour_columns=map(int,options.colour_columns.split(","))
			except StandardError:
				print "Could not understand columns selected:", options.columns
				columns=[1]
		else:
			columns=[1]
#		try: columns=map(int,options.columns.split(","))
#		except StandardError:
#			print "Could not understand columns selected:", options.columns
#			columns=[1]
		
		if options.colour_columns:
			try: 
				colour_columns=[]
				cols=options.colour_columns.split(",")
				for col in cols:
					if len(col.split(".."))==1:
						try:
							colour_columns.append(int(col))
						except StandardError:
							print "Could not understand columns selected:", options.colour_columns
							colour_columns=[]
					elif len(col.split(".."))==2:
						colbits=col.split("..")
						for x in xrange(int(colbits[0]),int(colbits[1])+1):
							colour_columns.append(x)
					else:
						print "Could not understand columns selected:", options.colour_columns
						colour_columns=[]
#			
#				colour_columns=map(int,options.colour_columns.split(","))
			except StandardError:
				print "Could not understand columns selected:", options.colour_columns
				colour_columns=[]
		else:
			colour_columns=[]
		
		print columns, colour_columns
		try: lines=open(options.metadata,"rU").readlines()
		except StandardError:
			print "Could not open metadatafile:", options.metadata
			sys.exit()
		colourslist=[]
		allcolourslist=[]
		
		if len(colour_columns)>0:
			colour_column_names=[]
			for column in colour_columns:
				if column<1:
					print "column numbers must be positive and greater than zero"
					sys.exit()
				colourslist.append([])
				if len(lines[0].strip().split(","))>=column:
					colour_column_names.append(lines[0].strip().split(",")[column-1])
				else:
					colour_column_names.append("")
		
		
		for line in lines[1:]:
			words=line.strip().split(",")
			
			if len(words)>0 and len(words[0].strip())>0:
				newname=[]
				for column in columns:
					if len(words)>=column and len(words[column-1].strip())>0:
						newname.append(words[column-1].strip())
				
				metadatanames[words[0].strip()]="/".join(newname)
				namecolours[words[0].strip()]={}
				if len(colour_columns)>0:
					for x, column in enumerate(colour_columns):
						try:
							colour_column_entry=int(words[column-1].strip())
						except ValueError:
							try:
								colour_column_entry=float(words[column-1].strip())
							except ValueError:
								colour_column_entry=words[column-1].strip()
						
						if colour_column_entry in ["", "-"]:
							continue
						
						if colour_column_entry not in colourslist[x]:
							colourslist[x].append(colour_column_entry)
						if colour_column_entry not in allcolourslist:
							allcolourslist.append(colour_column_entry)
						namecolours[words[0].strip()][x]=colour_column_entry
				
		found_keys=[]
		if len(colour_columns)>0:
			
			for x, colour_column in enumerate(colour_columns):
				colour_dict.append({})
				newtrack=Track()
				newtrack.track_height=1
				if options.end!=-1:
					newtrack.end=options.end
				newtrack.beginning=options.beginning
				
		#		newtrack.track_number=track_number
				newtrack.track_height=5
				newtrack.scale=False
				newtrack.name="metadata_key"+str(x)
				metadata_keylist.append("metadata_key"+str(x))
				newtrack.is_key=True
				if len(colour_column_names[x].split(":")[0])>0:
					newtrack.key_data=[["Key ("+colour_column_names[x].split(":")[0]+"):", colors.Color(0, 0, 0)]]
				else:
					newtrack.key_data=[["Key:", colors.Color(0, 0, 0)]]
				words=colour_column_names[x].split(":")
				if get_text_width("Helvetica", 10, words[0])>max_length_col_name:
					max_length_col_name=get_text_width("Helvetica", 10, words[0])
				
				if len(words)>1:
					if words[1] in ["C", "c"]:
						newtrack.datatype="continuous"
						if len(words)==4:
							try:
								newtrack.datamin=float(words[2])
							except StandardError:
								try:
									newtrack.datamin=min(map(float,colourslist[x]))
								except StandardError:
									newtrack.datatype="discrete"
							try:
								newtrack.datamax=float(words[3])
							except StandardError:
								try:
									newtrack.datamax=max(map(float,colourslist[x]))
								except StandardError:
									newtrack.datatype="discrete"
						else:
							try:
								newtrack.datamin=min(map(float,colourslist[x]))
							except StandardError:
								newtrack.datatype="discrete"
							try:
								newtrack.datamax=max(map(float,colourslist[x]))
							except StandardError:
								newtrack.datatype="discrete"
					elif words[1] in ["D", "d"]:
						newtrack.datatype="discrete"		
					else:
						newtrack.datatype="discrete"
				else:
					newtrack.datatype="discrete"
				
#				try:
#					colourslist[x]=map(float,colourslist[x])
#					colourslist[x].sort()
#				except StandardError:
#					try:
#						colourslist[x]=map(int,colourslist[x])
#						colourslist[x].sort()
#					except StandardError:
#						colourslist[x].sort()
				
				colourslist[x].sort()
				
				if "" in colourslist[x]:
					colour_dict[x][""]=(0,0,0)
					colourslist[x].remove("")
				
				
				if len(colourslist[x])==1 and newtrack.datatype=="discrete":
					newtrack.key_data.append([colourslist[x][0], colors.Color(1, 0, 0)])
					colour_dict[x][colourslist[x][0]]=(255, 0, 0)
				elif len(colourslist[x])==2 and newtrack.datatype=="discrete":
					newtrack.key_data.append([colourslist[x][0], colors.Color(0, 0, 1)])
					newtrack.key_data.append([colourslist[x][1], colors.Color(1, 0, 0)])
					colour_dict[x][colourslist[x][0]]=(0, 0, 255)
					colour_dict[x][colourslist[x][1]]=(255, 0, 0)
				elif newtrack.datatype=="continuous":
					
					for y, name in enumerate(colourslist[x]):
						value=name
						if value<newtrack.datamin:
							value=newtrack.datamin
						elif value>newtrack.datamax:
							value=newtrack.datamax
					
						proportion=((float(value)-newtrack.datamin)/((newtrack.datamax-newtrack.datamin)))*255
						
						
						red=proportion
						blue=255-proportion
						green=0
						colour_dict[x][name]=(red, green, blue)
					
					newtrack.key_data.append([newtrack.datamin, colors.Color(0, 0, 1)])
					newtrack.key_data.append(["===>", colors.Color(0, 0, 0)])
#					newtrack.key_data.append([newtrack.datamin, colors.Color(float(0)/255, float(green)/255, float(blue)/255)])
#					newtrack.key_data.append([newtrack.datamin, colors.Color(float(0)/255, float(green)/255, float(blue)/255)])
#					newtrack.key_data.append([newtrack.datamin, colors.Color(float(0)/255, float(green)/255, float(blue)/255)])
					newtrack.key_data.append([newtrack.datamax, colors.Color(1, 0, 0)])
				elif len(colourslist[x])>2:
					if newtrack.datatype=="discrete":
						for y, name in enumerate(colourslist[x]):
							#proportion=(float(x)/(len(colourslist[x])-1))*1275
							proportion=(float(y)/(len(colourslist[x])-1))*1175
							#proportion=1175-proportion
							
							red=510-proportion
							blue=proportion-510
							green=proportion
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
				#			if green>255 and green<765:
				#				green=255
							if green>155 and green<765:
								green=155
							elif green>=765:
								greendiff=765-green
								green=155+greendiff
								if green<0:
									green=0
	#						if not name in colour_dict[x]:
	#							colour_dict[x][name]=[]
							colour_dict[x][name]=(red, green, blue)
							newtrack.key_data.append([name, colors.Color(float(red)/255, float(green)/255, float(blue)/255)])
					
					
				if options.show_metadata_key and not colour_column_names[x].split(":")[0] in found_keys:
					track_count+=5	
					my_tracks["metadata_key"+str(x)]=newtrack
					found_keys.append(colour_column_names[x].split(":")[0])
	
	
	topmargin=topmargin+(float(max_length_col_name)/2)
	
	for arg in args[::-1]:
		if arg.lower() in ["tree", "list"]:
			input_order.append(arg.lower())
			continue
		if arg.split('.')[-1].lower() in ["plot", "hist", "heat", "bar", "line", "graph", "area", "stackedarea","embl", "gb", "gbk", "tab", "bam", "bcf", "fas", "fasta", "mfa", "dna", "fst", "phylip", "phy", "nexus", "nxs"]:
			
			if arg.split('.')[-1].lower() in ["plot", "hist", "heat", "bar", "line", "graph", "area", "stackedarea", "bam"] or options.qualifier=="":
				newtrack = Track()
				if options.suffix != "":
					namelen=len(arg.split("/")[-1])
					name='.'.join(arg.split("/")[-1].split(options.suffix)[:-1])
					if len(name)==namelen:
						name='.'.join(arg.split("/")[-1].split('.')[:-1])
				else:
					name='.'.join(arg.split("/")[-1].split('.')[:-1])
				if options.end!=-1:
					newtrack.end=options.end
				newtrack.beginning=options.beginning
			
			if arg.split('.')[-1].lower() in ["plot","graph"]:
				track_count+=options.plotheight
				plot_type=options.plottype
				newtrack.track_height=options.plotheight
				newtrack.scale=False
				newtrack.add_plot(arg, plot_type, options.fragments)
			elif arg.split('.')[-1].lower() in ["bam"]:
				track_count+=options.plotheight
				plot_type=options.plottype
				newtrack.track_height=options.plotheight
				newtrack.scale=False
				newtrack.add_bam_plot(arg, plot_type, options.fragments)
			elif arg.split('.')[-1].lower() in ["bcf"]:
				newtrack=add_bcf_to_diagram(arg)
				newtrack.track_height=1
				if options.end!=-1:
					newtrack.end=options.end
				newtrack.beginning=options.beginning
				track_count+=1
				
				newtrack.scale=False
				if options.suffix != "":
					namelen=len(arg.split("/")[-1])
					name='.'.join(arg.split("/")[-1].split(options.suffix)[:-1])
					if len(name)==namelen:
						name='.'.join(arg.split("/")[-1].split('.')[:-1])
				else:
					name='.'.join(arg.split("/")[-1].split('.')[:-1])
				newtrack.name=name
				x=1
				while name in my_tracks:
					name='.'.join(arg.split('.')[:-1])+"_"+str(x)
					x+=1
				if not newtrack.name in track_names:
					track_names[newtrack.name]=[]
				input_order.append(name)
				my_tracks[name]=newtrack
			elif arg.split('.')[-1].lower() in ["hist", "heat", "bar", "line", "area", "stackedarea"]:
				track_count+=options.plotheight
				newtrack.track_height=options.plotheight
				plot_type=arg.split('.')[-1].lower()
				
				newtrack.scale=False
				newtrack.add_plot(arg, plot_type, options.fragments)
			elif arg.split('.')[-1].lower() in ["embl", "gb", "gbk"]:
				track_count+=options.emblheight
				
				newtrack = Track()
				
				try:
					emblrecord=open_annotation(arg)
				except (StandardError, SimonError) as e:
					#print e.strerror
					DoError("Cannot open annotation file "+arg+" please check the format")
				name=emblrecord.name
				while name in my_tracks:
					name=emblrecord.name+"_"+str(x)
					x+=1
				if options.end!=-1:
					newtrack.end=options.end
				newtrack.beginning=options.beginning
				newtrack=add_embl_to_diagram(emblrecord, incfeatures=["CDS", "feature", "tRNA", "rRNA"], emblfile=True, name=name)
				
				newtrack.scale=True
				newtrack.name=name
				newtrack.track_height=options.emblheight
				
				if options.labels>0:
					track_count+=options.labels
					newtrack.feature_label_angle=options.label_angle
					newtrack.feature_label_track_height=options.labels
					newtrack.draw_feature_labels=True
				
				if not newtrack.name in track_names:
					track_names[newtrack.name]=[]
				input_order.append(name)
				track_names[newtrack.name].append(name)
				
				my_tracks[name]=newtrack
			
			elif arg.split('.')[-1].lower() in ["fas", "fasta", "mfa", "dna", "fst", "phylip", "phy", "nexus", "nxs"]:
				track_count+=options.emblheight
				
				newtrack = Track()
				try:
					fastarecord=read_seq_file(arg)
				except (StandardError, SimonError):
					DoError("Cannot open sequence file "+arg+" please check the format")
				shortname='.'.join(arg.split("/")[-1].split('.')[:-1])
				name=shortname
				while name in my_tracks:
					name=shortname+"_"+str(x)
					x+=1
				if options.end!=-1:
					newtrack.end=options.end
				newtrack.beginning=options.beginning
				newtrack=add_sequence_file_to_diagram(fastarecord, name)
				newtrack.scale=True
				
				newtrack.scale_position="middle"
				newtrack.name=name
				newtrack.track_height=options.emblheight
				if options.labels>0:
					track_count+=options.labels
					newtrack.feature_label_angle=options.label_angle
					newtrack.feature_label_track_height=options.labels
					newtrack.draw_feature_labels=True
				if not newtrack.name in track_names:
					track_names[newtrack.name]=[]
				input_order.append(name)
				track_names[newtrack.name].append(name)
				
				my_tracks[name]=newtrack
				
				
			elif arg.split('.')[-1].lower() =="tab":
			
				if options.qualifier!="":#(options.tree!="" or options.taxon_list!="") and options.qualifier!="":
					
					new_tracks=add_ordered_tab_to_diagram(arg)
					
					for track in new_tracks:
						newtrack=new_tracks[track]
						if options.end!=-1:
							newtrack.end=options.end
						newtrack.beginning=options.beginning
						newtrack.name=new_tracks[track].name
						name=newtrack.name
						x=1
						while name in my_tracks:
							name=newtrack.name+"_"+str(x)
							x+=1
						if not newtrack.name in track_names:
							track_names[newtrack.name]=[]
						input_order.append(name)
						track_names[newtrack.name].append(name)
						
						track_count+=1
						newtrack.track_height=1
						my_tracks[name]=newtrack
						
				
				else:
					newtrack=add_tab_to_diagram(arg)
					newtrack.track_height=1
					if options.end!=-1:
						newtrack.end=options.end
					newtrack.beginning=options.beginning
					track_count+=1
					
					newtrack.scale=False
					if options.suffix != "":
						namelen=len(arg.split("/")[-1])
						name='.'.join(arg.split("/")[-1].split(options.suffix)[:-1])
						if len(name)==namelen:
							name='.'.join(arg.split("/")[-1].split('.')[:-1])
					else:
						name='.'.join(arg.split("/")[-1].split('.')[:-1])
					newtrack.name=name
					x=1
					while name in my_tracks:
						name='.'.join(arg.split('.')[:-1])+"_"+str(x)
						x+=1
					if not newtrack.name in track_names:
						track_names[newtrack.name]=[]
					input_order.append(name)
					my_tracks[name]=newtrack
				
			
			if arg.split('.')[-1].lower() in ["plot", "hist", "heat", "bar", "line", "graph", "area", "stackedarea", "bam"]:# or (arg.split('.')[-1].lower()=="tab" and options.qualifier==""):
				newtrack.name=name
				x=1
				while name in my_tracks:
					name=newtrack.name+"_"+str(x)
					x+=1
				if not newtrack.name in track_names:
					track_names[newtrack.name]=[]
				input_order.append(name)
				track_names[newtrack.name].append(name)
				
				
				
				my_tracks[name]=newtrack
			

	treenames=[]
	tree_name_to_node={}
	listnames=[]
	
	if options.tree!="":
		if not os.path.isfile(options.tree):
			print "Cannot find file:", options.tree
			options.tree=""
		else:
			treestring=open(options.tree,"rU").read().strip()
			tree=Trees.Tree(treestring, rooted=True)
			if options.midpoint:
				tree=midpoint_root(tree)
#				treestring=tree_to_string(tree,plain=False,ladderize=options.ladderise)
#				tree=Trees.Tree(treestring, rooted=True)
				
			tree.root
			
			
			def get_ordered_terminals(treeObject):
				
				def measure_downstream_branches(basenode,count=0, maxlen=0.0):
					def measure_next_node(node,count, maxlen):
						for daughter in treeObject.node(node).succ:
							count, maxlen=measure_next_node(daughter,count, maxlen)
						count+=1
						
						if treeObject.is_terminal(node):
							nodedist=treeObject.distance(basenode,node)
							if nodedist>maxlen:
								maxlen=nodedist
						
						return count, maxlen
					count, maxlen=measure_next_node(basenode,count, maxlen)
					
					return count, maxlen+treeObject.node(basenode).data.branchlength
				
				
				def get_next_terminal(node, order):
				
					succs=treeObject.node(node).succ
										
					succ_counts=[]
					
					for succ in succs:
						brcount,brlen=measure_downstream_branches(succ)
						succ_counts.append([brcount,brlen,succ])
					
					succ_counts.sort()
					
					for succ_count in succ_counts:
						succ=succ_count[2]
						
						if not treeObject.is_terminal(succ):
							order=get_next_terminal(succ, order)
						else:
							order.append(succ)
							#print node, succ_count
					
					return order
					

				order=[]
				order=get_next_terminal(treeObject.root, order)
				return order
			
			
			def get_terminals(treeObject):
				
				def get_next_taxon(node, order):
				
					succs=treeObject.node(node).succ
										
					
					for succ in succs:
						
						if not treeObject.is_terminal(succ):
							order=get_next_taxon(succ, order)
						else:
							order.append(succ)
							#print node, succ_count
					
					return order
					

				order=[]
				order=get_next_taxon(treeObject.root, order)
				return order
				
			
			if options.ladderise in ["left", "right"]:
				order=get_ordered_terminals(tree)
				if options.ladderise=="left":
					order.reverse()
			else:
				#order=tree.get_terminals()
				order=get_terminals(tree)
			
			totalbr=0.0
			
			tree=get_tree_colour_comments(tree)
			
			for terminal_node in order:
				terminal=tree.node(terminal_node).data.taxon
				terminal=terminal.strip("'")
				terminal=terminal.strip('"')
				treenames.append(terminal)
				if not terminal in track_names:
					track_count+=1
				tree_name_to_node[terminal]=terminal_node
				
				if terminal in namecolours and len(namecolours[terminal])>0:
					tree.node(terminal_node).data.comment["name_colour"]=[]
					if len(colour_columns)>0:
						for x in xrange(len(colour_columns)):
							if x in namecolours[terminal]:
								namecolour=namecolours[terminal][x]
								tree.node(terminal_node).data.comment["name_colour"].append(colour_dict[x][namecolour])
							else:
								tree.node(terminal_node).data.comment["name_colour"].append((0,0,0))
								namecolours[terminal][x]="-"
					else:
						tree.node(terminal_node).data.comment["name_colour"].append((0,0,0))
						namecolours[terminal][0]="-"
				else:
					namecolours[terminal]={}
					tree.node(terminal_node).data.comment["name_colour"]=[]
					if len(colour_columns)>0:
						for x in xrange(len(colour_columns)):
							tree.node(terminal_node).data.comment["name_colour"].append((0,0,0))
							namecolours[terminal][x]="-"
					else:
						tree.node(terminal_node).data.comment["name_colour"].append((0,0,0))
						namecolours[terminal][0]="-"	
				
				if terminal in metadatanames:
					terminal=metadatanames[terminal]
					tree.node(terminal_node).data.taxon=terminal
				totalbr+=tree.sum_branchlength(root=tree.root, node=terminal_node)
			if totalbr==0:
				
				def make_branches_equal(node):
					for daughter in tree.node(node).succ:
						make_branches_equal(daughter)
						tree.node(daughter).data.branchlength=1
				make_branches_equal(tree.root)
		
			if options.transformation in ["acctran", "deltran"]:
				parsimony_reconstruction(tree, namecolours, colour_dict[0], transformation=options.transformation)
			
				
				
				
				
				
				
				
		
		
	elif options.taxon_list!="":
		if not os.path.isfile(options.taxon_list):
			print "Cannot find file:", options.taxon_list
			options.taxon_list=""
		else:
			for line in open(options.taxon_list,"rU"):
				taxonname=line.strip().split()[0]
				if taxonname!="":
					listnames.append(taxonname)
#					track_names.append(taxonname)
			for terminal in listnames:
				if not terminal in track_names:
					track_count+=1
					

	



	#Find the maximum and minimum plot values across all plots (will be used to calculate plot scale size for name offset)
	
	maxplotheight=float("-Inf")
	minplotheight=float("Inf")
	for track in my_tracks:
		for plot in my_tracks[track].plots:
			for data in plot.data:
				if len(data)==0:
					continue
				plotmax=max(data)#, key=lambda x: x[1])[1]
				if plotmax>maxplotheight:
					maxplotheight=plotmax
				plotmin=min(data)#, key=lambda x: x[1])[1]
				if plotmin<minplotheight:
					minplotheight=plotmin
	
	#if the user has asked all plots to be scaled the same, change the max and min for all plots
	#in all cases, work out the length of the max and min scale point text
	
	maxplot_scale_text=0
	for track in my_tracks:
		for plot in my_tracks[track].plots:
			
			if options.scale_plots_same:
				plot.max_yaxis=maxplotheight
				plot.min_yaxis=minplotheight
			elif plot.plot_type in ["bar", "stack"]:
				
				plot.min_yaxis = round_to_n(min(map(min,plot.data)), 2)
				
				maxsum=0.0
				for datum in plot.data:
					maxsum+=max(datum)
				plot.max_yaxis = round_to_n(maxsum, 2)
			else:
				if len(plot.data)==0:
					plot.data=[minplotheight]
				plot.max_yaxis=max(map(max,plot.data))
				plot.min_yaxis=min(map(min,plot.data))
			
			if options.plot_min!=float("Inf"):
				plot.min_yaxis=options.plot_min
			if options.plot_max!=float("Inf"):
				plot.max_yaxis=options.plot_max
				
					
			scalemax=round_to_n(plot.max_yaxis, 2)
			scalemin=round_to_n(plot.min_yaxis, 2)
			if get_text_width("Helvetica", 10, str(scalemax))>maxplot_scale_text:
				maxplot_scale_text=get_text_width("Helvetica", 10, str(scalemax))
			if get_text_width("Helvetica", 10, str(scalemin))>maxplot_scale_text:
				maxplot_scale_text=get_text_width("Helvetica", 10, str(scalemin))
			


	if track_count==0:
		print "Error: No data found to print"
		sys.exit()

	#Calculate the total number of tracks we need, taking into account the number of fragments requested

	if options.fragments>1:
		track_count+=options.fragment_separation
		
	track_fragment_count=track_count*options.fragments
	
	
	#from this we can work out a constant for the height of a track which takes into account the height of the page and margin sizes
	
	vertical_scaling_factor=float(height-(margin+topmargin))/(track_fragment_count)
	
	#to make sure names can be printed in the space of a track, we can scale the name to the same size as the vertical scaling factor, but limit it to 12pt so it doesn't get crazily big
	
	name_font_size=vertical_scaling_factor
	if name_font_size>12:
		name_font_size=12

	
	#We can now calculate the width of the longest track name, so we can allocate space to print the names
	
	max_name_width=0
	
	for track in my_tracks:
		my_tracks[track].name_size=name_font_size*my_tracks[track].track_height
		if my_tracks[track].name_size>12:
			my_tracks[track].name_size=12
		
		if my_tracks[track].name in metadatanames:	
			my_tracks[track].name = metadatanames[my_tracks[track].name]
		
		my_tracks[track].name_length=get_text_width("Helvetica", my_tracks[track].name_size, my_tracks[track].name)+maxplot_scale_text+3
		if (my_tracks[track].name_length+my_tracks[track].name_size)>max_name_width:
			max_name_width=(my_tracks[track].name_length+my_tracks[track].name_size)
	
	
	#from the name width we can calculate the proportion of the printable image that will be take up by the name, or set it to zero if the user doesn't want to print names
	
	if options.track_names:
		name_proportion=float(max_name_width)/(width-(margin*2))
	else:
		name_proportion=0
	
	#if we have a tree to put in, we compare the name and tree widths and use the largest as the proportion of the page to leave for names/trees
	
	if options.treeproportion<name_proportion or options.tree=="":
		left_proportion=name_proportion
	else:
		left_proportion=options.treeproportion
	
	
	#Set the orders of the tracks	
	
	output_order=[]
	if options.show_metadata_key:
		if len(metadata_keylist)>0:
			metadata_keylist.reverse()
			output_order=metadata_keylist
	treetrack=0
	if not "tree" in input_order:
		output_order=output_order+treenames[::-1]
	if not "list" in input_order:
		output_order=output_order+listnames[::-1]
	
	
	added_tree=False
	added_list=False
	for name in input_order[::-1]:
		if name=="tree" and not added_tree:
			treetrack=len(output_order)
			added_tree=True
			output_order=output_order+treenames[::-1]
		elif name=="list" and not added_list:
			added_list=True
			output_order=output_order+listnames[::-1]
		elif not name in treenames+listnames:
			output_order.append(name)
		
	
	
	
	track_number=0
	
	for track in output_order:
	
		if not track in my_tracks:
			newtrack=Track()
			newtrack.track_height=1
			if options.end!=-1:
				newtrack.end=options.end
			newtrack.beginning=options.beginning
			#track_count+=1
			newtrack.track_number=track_number
			newtrack.track_height=1
			newtrack.scale=False
			name=track
			newtrack.name=name
			x=1
			while name in my_tracks:
				name=track+"_"+str(x)
				x+=1
			if not newtrack.name in track_names:
				track_names[newtrack.name]=[]
			my_tracks[name]=newtrack
#			if track in treenames:
#				tree.node(tree_name_to_node[track]).data.comment["vertpos"]=margin+((track_number)*vertical_scaling_factor)+float(vertical_scaling_factor)/2
#				my_tracks[track]
#				if track in namecolours:
#					r,g,b=colour_dict[namecolours[track]]
#					my_tracks[track].grey_track_colour=colors.Color(float(r)/255,float(g)/255,float(b)/255)
#			track_number+=1
#			continue
		track_height=my_tracks[track].track_height
		label_height=my_tracks[track].feature_label_track_height
		
		my_tracks[track].track_draw_proportion=options.tracksize
		my_tracks[track].track_height=track_height*vertical_scaling_factor
		if left_proportion==1:
			my_tracks[track].track_length=(width-margin)-((width-(margin*2))*0.2+margin)
			my_tracks[track].track_position=[(width-(margin*2))*0.2+margin, margin+((track_number)*vertical_scaling_factor)+float((my_tracks[track].track_height)/2)]

		else:
			my_tracks[track].track_length=(width-margin)-((width-(margin*2))*left_proportion+margin)
			my_tracks[track].track_position=[(width-(margin*2))*left_proportion+margin, margin+((track_number)*vertical_scaling_factor)+float((my_tracks[track].track_height)/2)]

		my_tracks[track].track_number=track_number
		if options.track_names and not track in treenames:
			my_tracks[track].show_name=True
		if track in treenames:
			tree.node(tree_name_to_node[track]).data.comment["vertpos"]=margin+((track_number)*vertical_scaling_factor)+float((my_tracks[track].track_height)/2)
			if track in namecolours and '0' in namecolours[track] and namecolours[track][0]!="-":
				r,g,b=colour_dict[0][namecolours[track][max(namecolours[track].keys())]]
				my_tracks[track].grey_track_colour=colors.Color(float(r)/255,float(g)/255,float(b)/255)
			else:
				r,g,b=(0,0,0)
				my_tracks[track].grey_track_colour=colors.Color(float(r)/255,float(g)/255,float(b)/255)
		
		track_number+=track_height+label_height
		


	#find the maximum feature endpoint to scale by
	max_feature_length=0
	
	for track in my_tracks:
		max_track_feature_length=my_tracks[track].get_max_feature_length()
		if max_track_feature_length>max_feature_length:
			max_feature_length=max_track_feature_length
		for plot in my_tracks[track].plots:
			for data in plot.xdata:
				if data[-1]>max_feature_length:
					max_feature_length=data[-1]
		
			
	for track in my_tracks:
		if my_tracks[track].max_feature_length<max_feature_length:
				my_tracks[track].max_feature_length=max_feature_length
#		for plot in my_tracks[track].plots:
#			for data in plot.xdata:
#				print len(data), data[-10:]
#			for data in plot.data:
#				print len(data), data[-10:]
	

	if options.end!=-1:
		length=options.end-options.beginning
	else:
		length=max_feature_length-options.beginning
	
	
	#We can now start to print the tracks
	c = Canvas(options.outputfile, pagesize=(width, height))
	
	pagelen=float(length)/options.npages
	for page in xrange(options.npages):
		if options.npages>1:
			print "Printing page", page+1
		else:
			print "Printing figure"
		d = Drawing(width, height)
		for fragment in xrange(1, options.fragments+1):
			
			beginning=int(((float(pagelen)/options.fragments)*(fragment-1))+options.beginning)+(page*pagelen)
			end=int(((float(pagelen)/options.fragments)*fragment)+options.beginning)+(page*pagelen)
			#print beginning, end
			for track in output_order:
				
				if not track in my_tracks or (my_tracks[track].is_key and fragment!=options.fragments):
					continue
				
				my_tracks[track].beginning=beginning
				my_tracks[track].end=end
				
				
				
				my_tracks[track].track_position[1]=margin+(((((options.fragments-fragment)*track_count)+my_tracks[track].track_number)*vertical_scaling_factor)+(my_tracks[track].track_height)/2)
				
				if options.greytracks:
					my_tracks[track].greytrack=True
				my_tracks[track].sort_features_by_length()
				my_tracks[track].draw_track()
			
			if options.tree!="":
				drawtree(tree, height-(margin*2), (width-(margin*2))*left_proportion, margin, ((((options.fragments-fragment)*track_count)*vertical_scaling_factor)), maxplot_scale_text+3)
				
		
		
		 
		renderPDF.draw(d, c, 0, 0)
		c.showPage()
	c.save()
	print "Done"
	
	
	
