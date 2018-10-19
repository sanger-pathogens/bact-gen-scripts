#!/usr/bin/env python
##!/software/python-2.7.6/bin/python
#################################
# Import some necessary modules #
#################################

import os, sys
sys.path.insert(1,'/software/python-2.7.6/lib/python2.7/site-packages/reportlab-3.1.8-py2.7-linux-x86_64.egg')
sys.path.insert(1,'/software/python-2.7.6/lib/python2.7/site-packages/')
sys.path.insert(1,'/software/python-2.7.6/lib/python2.7/')
sys.path.insert(1, '/nfs/users/nfs_s/sh16/scripts/modules/')
import dendropy
import string, re
import random
from math import sqrt, pow, log, floor, sin, log10, ceil
from numpy import repeat, convolve, mean, median
from optparse import OptionParser, OptionGroup
from Bio.Nexus import Trees, Nodes
import shlex, subprocess
from colorsys import hsv_to_rgb
#on my laptop
#sys.path.extend(map(os.path.abspath, ['/Users/sh16/Documents/scripts/modules/']))
#on pcs4
#sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_general import *
from Si_SeqIO import *
#from Si_nexus import midpoint_root, tree_to_string
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import GenBank
from Bio.Graphics.GenomeDiagram._Colors import ColorTranslator
from math import floor
import imp
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
from itertools import izip
import datetime
#on my laptop
#SAMTOOLS_DIR="/Users/sh16/Applications/samtools-0.1.18/"
#BCFTOOLS_DIR="/Users/sh16/Applications/samtools-0.1.18/bcftools/"
#on pcs4
OLD_SAMTOOLS_DIR="/nfs/users/nfs_s/sh16/samtools-0.1.17/"
OLD_BCFTOOLS_DIR="/nfs/users/nfs_s/sh16/samtools-0.1.17/bcftools/"
SAMTOOLS_DIR=""
BCFTOOLS_DIR=""

#################################
# Set up some global converters #
#################################

pagesizeconverter={'A0':pagesizes.A0, 'A1':pagesizes.A1, 'A2':pagesizes.A2, 'A3':pagesizes.A3, 'A4':pagesizes.A4, 'A5':pagesizes.A5, 'A6':pagesizes.A6, 'B0':pagesizes.B0, 'B1':pagesizes.B1, 'B2':pagesizes.B2, 'B3':pagesizes.B3, 'B4':pagesizes.B4, 'B5':pagesizes.B5, 'B6':pagesizes.B6, 'LEGAL':pagesizes.LEGAL, 'LETTER':pagesizes.LETTER, 'landscape':pagesizes.landscape, 'legal':pagesizes.legal, 'letter':pagesizes.letter, 'portrait':pagesizes.portrait}


colourconverter={'aliceblue':colors.aliceblue, 'antiquewhite':colors.antiquewhite, 'aqua':colors.aqua, 'aquamarine':colors.aquamarine, 'azure':colors.azure, 'beige':colors.beige, 'bisque':colors.bisque, 'black':colors.black, 'blanchedalmond':colors.blanchedalmond, 'blue':colors.blue, 'blueviolet':colors.blueviolet, 'brown':colors.brown, 'burlywood':colors.burlywood, 'cadetblue':colors.cadetblue, 'chartreuse':colors.chartreuse, 'chocolate':colors.chocolate, 'coral':colors.coral, 'cornflower':colors.cornflower, 'cornflowerblue':colors.cornflowerblue, 'cornsilk':colors.cornsilk, 'crimson':colors.crimson, 'cyan':colors.cyan, 'darkblue':colors.darkblue, 'darkcyan':colors.darkcyan, 'darkgoldenrod':colors.darkgoldenrod, 'darkgray':colors.darkgray, 'darkgreen':colors.darkgreen, 'darkgrey':colors.darkgrey, 'darkkhaki':colors.darkkhaki, 'darkmagenta':colors.darkmagenta, 'darkolivegreen':colors.darkolivegreen, 'darkorange':colors.darkorange, 'darkorchid':colors.darkorchid, 'darkred':colors.darkred, 'darksalmon':colors.darksalmon, 'darkseagreen':colors.darkseagreen, 'darkslateblue':colors.darkslateblue, 'darkslategray':colors.darkslategray, 'darkslategrey':colors.darkslategrey, 'darkturquoise':colors.darkturquoise, 'darkviolet':colors.darkviolet, 'deeppink':colors.deeppink, 'deepskyblue':colors.deepskyblue, 'dimgray':colors.dimgray, 'dimgrey':colors.dimgrey, 'dodgerblue':colors.dodgerblue, 'fidblue':colors.fidblue, 'fidlightblue':colors.fidlightblue, 'fidred':colors.fidred, 'firebrick':colors.floralwhite, 'floralwhite':colors.floralwhite, 'forestgreen':colors.forestgreen, 'fuchsia':colors.fuchsia, 'gainsboro':colors.gainsboro, 'ghostwhite':colors.ghostwhite, 'gold':colors.gold, 'goldenrod':colors.goldenrod, 'gray':colors.gray, 'green':colors.green, 'greenyellow':colors.greenyellow, 'grey':colors.grey, 'honeydew':colors.honeydew, 'hotpink':colors.hotpink, 'indianred':colors.indianred, 'indigo':colors.indigo, 'ivory':colors.ivory, 'khaki':colors.khaki, 'lavender':colors.lavender, 'lavenderblush':colors.lavenderblush, 'lawngreen':colors.lawngreen, 'lemonchiffon':colors.lemonchiffon, 'lightblue':colors.lightblue, 'lightcoral':colors.lightcoral, 'lightcyan':colors.lightcyan, 'lightgoldenrodyellow':colors.lightgoldenrodyellow, 'lightgreen':colors.lightgreen, 'lightgrey':colors.lightgrey, 'lightpink':colors.lightpink, 'lightsalmon':colors.lightsalmon, 'lightseagreen':colors.lightseagreen, 'lightskyblue':colors.lightskyblue, 'lightslategray':colors.lightslategray, 'lightslategrey':colors.lightslategrey, 'lightsteelblue':colors.lightsteelblue, 'lightyellow':colors.lightyellow, 'lime':colors.lime, 'limegreen':colors.limegreen, 'linen':colors.linen, 'magenta':colors.magenta, 'maroon':colors.maroon, 'mediumaquamarine':colors.mediumaquamarine, 'mediumblue':colors.mediumblue, 'mediumorchid':colors.mediumorchid, 'mediumpurple':colors.mediumpurple, 'mediumseagreen':colors.mediumseagreen, 'mediumslateblue':colors.mediumslateblue, 'mediumspringgreen':colors.mediumspringgreen, 'mediumturquoise':colors.mediumturquoise, 'mediumvioletred':colors.mediumvioletred, 'midnightblue':colors.midnightblue, 'mintcream':colors.mintcream, 'mistyrose':colors.mistyrose, 'moccasin':colors.moccasin, 'navajowhite':colors.navajowhite, 'navy':colors.navy , 'oldlace':colors.oldlace, 'olive':colors.olive, 'olivedrab':colors.olivedrab, 'orange':colors.orange, 'orangered':colors.orangered, 'orchid':colors.orchid, 'palegoldenrod':colors.palegoldenrod, 'palegreen':colors.palegreen, 'paleturquoise':colors.paleturquoise, 'palevioletred':colors.palevioletred, 'papayawhip':colors.papayawhip, 'peachpuff':colors.peachpuff, 'peru':colors.peru, 'pink':colors.pink, 'plum':colors.plum, 'powderblue':colors.powderblue, 'purple':colors.purple, 'red':colors.red, 'rosybrown':colors.rosybrown, 'royalblue':colors.royalblue, 'saddlebrown':colors.saddlebrown, 'salmon':colors.salmon, 'sandybrown':colors.sandybrown, 'seagreen':colors.seagreen, 'seashell':colors.seashell, 'sienna':colors.sienna, 'silver':colors.silver, 'skyblue':colors.skyblue, 'slateblue':colors.slateblue, 'slategray':colors.slategray, 'slategrey':colors.slategrey, 'snow':colors.snow, 'springgreen':colors.springgreen, 'steelblue':colors.steelblue, 'tan':colors.tan, 'teal':colors.teal, 'thistle':colors.thistle, 'tomato':colors.tomato, 'turquoise':colors.turquoise, 'violet':colors.violet, 'wheat':colors.wheat, 'white':colors.white, 'whitesmoke':colors.whitesmoke, 'yellow':colors.yellow, 'yellowgreen':colors.yellowgreen}


default_rgb_colours=["blue", "red",  "limegreen", "aqua", "fuchsia", "orange", "green", "gray", "purple", "olive", "teal", "silver", "navy", "white", "black", "maroon", "yellow"]
default_plot_colours=["red", "blue",  "limegreen", "aqua", "fuchsia", "orange", "green", "gray", "purple", "olive", "teal", "silver", "navy", "black", "maroon", "yellow"]




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
	group.add_option("-s", "--support", action="store", dest="tree_support", help="Scale tree branch widths by value. For newick trees this can be any value stored in the tree. Otherwise, use 'support' to scale by branch support values (if present)", default="")
	group.add_option("-7", "--height_HPD", action="store_true", dest="show_height_HPD", help="show branch 95% HPD heights (if present in tree) [default= %default]", default=False)
	group.add_option("-6", "--brlabels", action="store", dest="branch_labels", help="Label branches with value. For newick trees this can be any value stored in the tree. Otherwise, use 'support' to label with branch support values (if present), or 'brlen' to label with branch lengths.", default="")
	group.add_option("-M", "--midpoint", action="store_true", dest="midpoint", help="Midpoint root tree", default=False)
	group.add_option("-L", "--ladderise", action="store", choices=['right', 'left'], dest="ladderise", help="ladderise tree (choose from right or left) [default= %default]", type="choice", default=None)
	group.add_option("-z", "--names_as_shapes", action="store", choices=['circle', 'square', 'rectangle', 'auto'], dest="names_as_shapes", help="Use shapes rather than taxon names in tree (choose from circle) [default= %default]", type="choice", default="auto")
	group.add_option("-1", "--logbranches", action="store_true", dest="log_branches", help="log branch lengths [default= %default]", default=False)
	group.add_option("-D", "--date_separator", action="store", dest="date_separator", help="For trees with dating information, the script will read these from the taxon names if they are a suffix separated from the rest of the taxon name by a particular character. To do this you need to choose the character used to separate the dates with this option", default="")
	group.add_option("-8", "--most_recent_isolation_date", action="store", dest="most_recent", help="Most recent isolation date for samples in tree. This must be in the same units used in the BEAST tree, and as defined by the -I option.", default=0, type="float")
	group.add_option("-I", "--time_type", action="store", choices=['days', 'years'], dest="time_type", help="Unit of time used in trees with dating information (choose from days or years) [default= %default]", type="choice", default="years")
	group.add_option("-j", "--Figtree_colours", action="store_true", dest="figtree", help="Use colours in tree saved from Figtree [default= %default]", default=False)
	group.add_option("-J", "--colour_by_annotation", action="store", dest="colour_by_annotation", help="Colour tree based on annotations (e.g. BEAST trait states). You must specify the annotation to use.", default="")

	
	parser.add_option_group(group)
	
	
	group = OptionGroup(parser, "Output Options For Tracks Without Tree")
	
	group.add_option("-x", "--taxon_list", action="store", dest="taxon_list", help="File with ordered taxa", default="")
	group.add_option("-n", "--names", action="store_true", dest="track_names", help="Show track names", default=False)
	
	parser.add_option_group(group)
	
	
	group = OptionGroup(parser, "Names Options")
	
	group.add_option("-a", "--aligntaxa", action="store", dest="aligntaxa", help="Align taxon names (0=no, 1=align left, 2=align right, 3=name on tree but metadata aligned to right) [default= %default]", default=0, type="int")
	group.add_option("-N", "--taxanames", action="store_false", dest="taxon_names", help="Do not show taxa names on tree. [default=show names]", default=True)
	group.add_option("-m", "--metadata", action="store", dest="metadata", help="metadata file in csv format. Note that this file must have a header row ontaining titles for each column", default="")
	group.add_option("-c", "--columns", action="store", dest="columns", help="column(s) from metadata file to use for track name (comma separated list) [default=%default]", default="1")
	group.add_option("-C", "--colourbycolumns", action="store", dest="colour_columns", help="column(s) from metadata file to use to colour track name and blocks next to name (comma separated list). If names are being shown, the first column will be used to colour names. All following columns will be added as coloured shapes as defined by the -z option. [default=%default]", default=False)
	group.add_option("-k", "--nometadatakey", action="store_false", dest="show_metadata_key", help="Do not show metadata keys. [default=show keys]", default=True)
	group.add_option("-r", "--parsimony_reconstruction", action="store", dest="transformation", help="Reconstruct colours across branches using parsimony. Select from acctran or deltran transformations [default=%default]", default=None, type="choice", choices=['acctran', 'deltran'])
	group.add_option("-i", "--suffix", action="store", dest="suffix", help="suffix to remove from filenames", default="")
	group.add_option("-W", "--listcolours", action="store_true", dest="printcolours", help="Print list of available colours for columns and exit", default=False)
	
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
	group.add_option("-y", "--plot_max", action="store", dest="plot_max", help="Set a maximum value for plots. Can be a number or coverage relative to the median (e.g. 2x)", default="Inf", type="string")
	group.add_option("-Y", "--plot_min", action="store", dest="plot_min", help="Set a minimum value for plots. Can be a number or coverage relative to the median (e.g. 1x)", default="Inf", type="string")
	group.add_option("-Z", "--log_plot", action="store_true", dest="log_plots", help="Show plots on a log scale (doesn't work yet)", default=False)
	group.add_option("-w", "--plot_scales", action="store_true", dest="plot_scales", help="Show legend on heatmaps", default=False)
	group.add_option("-3", "--heatmap_colour", default="bluered", choices=["redblue", "bluered", "blackwhite", "whiteblack", "redandblue"], type="choice", action="store", dest="heat_colour", help="Set the default plot type (for plots called plot or graph and bam files). Choose from "+", ".join(["redblue", "bluered", "blackwhite", "whiteblack"]))
	group.add_option("-4", "--windows", action="store", dest="windows", help="Number of windows per line (be careful, making this value too high will make things very slow and memory intensive) [default= %default]", default=501, type="float")
	
	parser.add_option_group(group)
	
	return parser.parse_args()



###############################################
# Check the command line options are sensible #
###############################################


#############################
# Remove values from a list #
#############################

def remove_zeros_from_list(the_list):
   return [value for value in the_list if value != 0]



##########################################
# Function to read in Excel date formats #
##########################################

def xldate_as_datetime(xldate, datemode=0):
	
	if datemode not in (0, 1):
		print "Error with Excel date mode", xldate
		sys.exit()
		
	xldays = int(xldate)+1

	if xldays >= 1000000:
		print "Error with Excel date", xldate
		print "Cannot handle dates above 25/11/4637 (equals 999999 days post 1900)"
		sys.exit()
	
	if (xldays + 693594 + (1462 * datemode))<=0:
		print "Error with Excel date", xldate
		print "Cannot handle dates below 1/1/1 (equals 693594 days pre 1900)"
		sys.exit()
	
	return datetime.datetime.fromordinal(xldays + 693594 + (1462 * datemode))



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
		
	
	def print_count_of_changes_to_state(treeObject):
		node=treeObject.root
		changes={}
		print "here"
		def get_node_change(node):
			for daughter in treeObject.node(node).get_succ():
				if treeObject.node(node).get_data().comment["branch_colour"]!= treeObject.node(daughter).get_data().comment["branch_colour"]:
					if not treeObject.node(daughter).get_data().comment["branch_colour"] in changes:
						changes[treeObject.node(daughter).get_data().comment["branch_colour"]]=0
					changes[treeObject.node(daughter).get_data().comment["branch_colour"]]+=1
				get_node_change(daughter)
		
		print_node_sankoff(node)
		for change in changes:
			print change, changes[change]



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
	print_count_of_changes_to_state(treeObject)
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
	if hasattr(feature, "sub_features") and len(feature.sub_features)>0:
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
# Function to find the best font size for a string to fit in a certain length #
###############################################################################

def set_text_width(font, size, length, text):
	c = Canvas(test,pagesize=pagesize)
	strlen=float("Inf")
	while strlen>length:
		strlen= c.stringWidth(text,font,size)
		size-=0.1
	return size



##########################################
# Function to produce colours for a list #
##########################################

def get_colours_for_list(data_list, data_type="discrete", data_max=float("-Inf"), data_min=float("Inf"), direction="anticlockwise", start_angle=0.0, end_angle=240, start_saturation=0.8, end_saturation=0.8, start_value=0.9, end_value=0.9):
	s_start=control.metadata_colour_start_saturation
	s_end=control.metadata_colour_end_saturation
	v_start=control.metadata_colour_start_value
	v_end=control.metadata_colour_end_value
#	start_angle=control.metadata_colour_start_angle
#	end_angle=control.metadata_colour_end_angle
#	direction=control.metadata_colour_direction
	
	print start_angle, end_angle, direction
	
	colour_db={}
	if data_type=="continuous" and data_max==float("-Inf"):
		data_max=max(data_list)
	if data_type=="continuous" and data_min==float("Inf"):
		data_min=min(data_list)
	
	if direction=="clockwise" and start_angle>end_angle:
		rotation_degrees=start_angle-end_angle
		direction_multiplier=-1
	elif direction=="clockwise":
		rotation_degrees=start_angle+(360-end_angle)
		direction_multiplier=-1
	elif direction in ["anticlockwise", "counterclockwise"] and start_angle>end_angle:
		rotation_degrees=(360-start_angle)+end_angle
		direction_multiplier=1
	else:
		rotation_degrees=end_angle-start_angle
		direction_multiplier=1

	
	if len(data_list)==1 and data_type=="discrete":
		colour_db[data_list[0]]=colors.Color(1, 0, 0)
	elif len(data_list)==2 and data_type=="discrete":
		colour_db[data_list[0]]=colors.Color(0, 0, 1)
		colour_db[data_list[1]]=colors.Color(1, 0, 0)
	elif data_type=="continuous":
		
		for y, name in enumerate(data_list):
			value=name
			if value<data_min:
				value=data_min
			elif value>data_max:
				value=data_max
		
			proportion=((float(value)-data_min)/((data_max-data_min)))
			
			h=(start_angle/360)+(direction_multiplier*(((proportion/360)*rotation_degrees)))
			
			v=v_start+(proportion*(v_end-v_start))
			s=s_start+(proportion*(s_end-s_start))
			red, green, blue = hsv_to_rgb(h,s,v)
			colour_db[name]=colors.Color(float(red), float(green), float(blue))
		
	elif len(data_list)>2:
		if data_type=="discrete":
			for y, name in enumerate(data_list):
				proportion=(float(y)/(len(data_list)-1))
	
				h=(start_angle/360)+(direction_multiplier*(((proportion/360)*rotation_degrees)))
				v=v_start+(proportion*(v_end-v_start))
				s=s_start+(proportion*(s_end-s_start))
				red, green, blue = hsv_to_rgb(h,s,v)
				colour_db[name]=colors.Color(float(red), float(green), float(blue))
	
	return colour_db




##############################################
# Function to convert RGB ints to RGB tuples #
##############################################

def rgbint2rgbtuple(RGBint):
	blue =  RGBint & 255
	green = (RGBint >> 8) & 255
	red =   (RGBint >> 16) & 255
	return (red, green, blue)



def hex_to_rgb(value):
	value = value.lstrip('#')
	lv = len(value)
	return tuple(int(value[i:i+lv/3], 16) for i in range(0, lv, lv/3))


##############################################
# Function toread a tree file using dendropy #
##############################################

def read_dendropy_tree(treefile):
		
		def log_branchlengths(t):
		
			print "Log transforming branch lengths"
			
			for node in t.postorder_node_iter():
				if node.edge_length!=None:
					node.edge_length=log(node.edge_length+1)
		
		
		#Some fixes to dendropy to make it read annotations properly
		dendropy.dataio.nexustokenizer.NHX_COMMENT_FIELD_PATTERN = re.compile(r'(.+?)=({.+?}|.+?)(,|$)')
		dendropy.dataio.nexustokenizer.FIGTREE_COMMENT_FIELD_PATTERN = re.compile(r'(.+?)=({.+?}|.+?)(,|$)')
		
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
			
				if str(cur_node.edge.length)==str(plen):
					break_on_node = cur_node
					#FIX
					break         #when find the  midpoint, it should break the loop
				elif cur_node.edge.length > plen:
					target_edge = cur_node.edge
					head_node_edge_len = plen #cur_node.edge.length - plen
					plen = 0
					break
				elif cur_node.edge.length < plen:
					plen -= cur_node.edge.length
					cur_node = cur_node.parent_node
				else:
					print "Error midpoint rooting tree"
					sys.exit()
			
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
		
		dendropy.dataobject.Tree.reroot_at_midpoint=reroot_at_midpoint
		
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
		
		
		
#		for node in t.postorder_node_iter():
#			print node.label, node.edge.length
		
		
		#log the branch lengths if the option has been chosen
		if options.log_branches:
			log_branchlengths(t)
		
		
		#some code here to make the branch lengths equal - not actually possible to choose this yet
		t.draw_scale=True
		if options.branch_labels=="brlen":
			t.draw_scale=False
			
		
		equal_branches=False
		if equal_branches:
			for node in t.postorder_node_iter():
				if node.edge_length!=None:
					node.edge_length=0.1
			t.draw_scale=False
			
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
		
		#Make sure the tree is rooted on an edge rather than a node
		if t.is_unrooted:
			print "Tree is unrooted. Rooting it now."
			t.is_rooted=True
			t.update_splits(delete_outdegree_one=False)
		root = t.seed_node
		root_children = root.child_nodes()
		
		if len(root_children) != 2:
			print "Tree rooted at node. Rerooting on first edge from that node."
			t.reroot_at_edge(root_children[0].edge, update_splits=True)
		#print the tree in ascii as a cladogram
		#print(t.as_ascii_plot())
		#print the tree in ascii including branch lengths
		#print(t.as_ascii_plot(plot_metric='length'))
		
		
		for node in t.postorder_node_iter():
			if node.edge_length==None:
				node.edge_length=0
		
		
		for node in t.postorder_node_iter():
			for x, a in enumerate(node.annotations):
				if isinstance(a.value, str):
					a.value=a.value.replace('"','')
					try:
						node.annotations[x].value=float(a.value)
					except:
						node.annotations[x].value=a.value
#					print node.annotations[x]
				elif isinstance(a.value, list):
					for y in xrange(len(a.value)):
						if isinstance(a.value[y], str):
							a.value[y]=a.value[y].replace('"','')
							node.annotations[x].value[y]=a.value[y]
					
					try:
						node.annotations[x].value=map(float,node.annotations[x].value)
					except:		
						if a.name=="!hilight" and len(a.value)==3 and len(a.value[2])>1 and a.value[2][0]=="#":
							try:
								rgbint=int(a.value[2][1:])
								r,g,b=rgbint2rgbtuple(rgbint)
							except:
								try:
									r,g,b=hex_to_rgb(a.value[2])
								except:
									print a.value[1:]
									break
							node.annotations[x].name="Figtree_hilight"
							node.annotations[x].value=colors.Color(float(r)/255,float(g)/255,float(b)/255, 0.5)
				
				
				if isinstance(a.value, str):
					if a.name=="!color" and len (a.value)>1 and a.value[0]=="#":
						try:
							rgbint=int(a.value[1:])
							r,g,b=rgbint2rgbtuple(rgbint)
							if g!=0:
								print r,g,b, a.value[1:]
						except:
							try:
								r,g,b=hex_to_rgb(a.value)
							except:
								print "Figtree hilight value not a valid colour:", a.value[1:]
								break
						node.annotations[x].name="Figtree_colour"
						node.annotations[x].value=colors.Color(float(r)/255,float(g)/255,float(b)/255)
					
			
		#parse leaf annotation information
		for leaf in t.leaf_iter():
			#print leaf.taxon
			if hasattr(leaf, "taxon"):
				if hasattr(leaf.taxon, "annotations"):
					for x, a in enumerate(leaf.taxon.annotations):
						if isinstance(a.value, str):
							a.value=a.value.replace('"','')
							try:
								leaf.taxon.annotations[x].value=float(a.value)
							except:
								leaf.taxon.annotations[x].value=a.value
						elif isinstance(a.value, list):
							for y in xrange(len(a.value)):
								if isinstance(a.value[y], str):
									a.value[y]=a.value[y].replace('"','')
									leaf.taxon.annotations[x].value[y]=a.value[y]
							
							try:
								leaf.taxon.annotations[x].value=map(float,leaf.taxon.annotations[x].value)
							except:
								continue
						
#						print leaf.taxon.label, a.name, a.value
						if isinstance(a.value, str) and a.name=="!color" and len (a.value)>1 and a.value[0]=="#":
							try:
								rgbint=int(a.value[1:])
								r,g,b=rgbint2rgbtuple(rgbint)
							except:
								try:
									r,g,b=hex_to_rgb(a.value)
								except:
									print a.value[1:]
									break
							leaf.taxon.annotations[x].name="Figtree_colour"
							leaf.taxon.annotations[x].value=colors.Color(float(r)/255,float(g)/255,float(b)/255)
							
			
		
		#Check if tree is from BEAST (not sure best way to check this, but will check from height on root node)
		root=t.seed_node
		t.is_BEAST=False
		if hasattr(root, "annotations"):
			for a in root.annotations:
				if a.name=="height" or a.name=="CI":
					t.is_BEAST=True
					break
		
		
		if t.is_BEAST and options.date_separator!="":
			min_leaf_sampling_date=float("Inf")
			max_leaf_sampling_date=float("-Inf")
			for leaf in t.leaf_iter():
				try:
					last_part_of_name=str(leaf.taxon).split(options.date_separator)[-1]
				except:
					"Found name that failed to split on", options.date_separator+":", leaf.taxon
				
				try:
					leaf.sampling_date=int(float(last_part_of_name))
					if leaf.sampling_date>max_leaf_sampling_date:
						max_leaf_sampling_date=leaf.sampling_date
					if leaf.sampling_date<min_leaf_sampling_date:
						min_leaf_sampling_date=leaf.sampling_date
				except:
					if last_part_of_name=="NA":
						leaf.sampling_date=None
					else:
						"Failed to find date when splitting on", options.date_separator+":", leaf.taxon
			if min_leaf_sampling_date!=float("Inf") and max_leaf_sampling_date!=float("-Inf"):
				t.min_leaf_sampling_date=min_leaf_sampling_date
				t.max_leaf_sampling_date=max_leaf_sampling_date
			else:
				print "Failed to find any dates for taxa"
		
		if t.is_BEAST and options.most_recent>0:
			t.max_leaf_sampling_date=options.most_recent
			min_height=float("Inf")
			max_height=float("-Inf")
			for leaf in t.leaf_iter():
				if hasattr(leaf, "annotations"):
					for a in leaf.annotations:
						if a.name=="height":
							if a.value>max_height:
								max_height=a.value
							elif a.value<min_height:
								min_height=a.value
			t.min_leaf_sampling_date=options.most_recent-max_height
		
		
		
		#colour branches by annotation
		if options.colour_by_annotation!="":
			
			annotation_list=[]
			
			for node in t.postorder_node_iter():
				if hasattr(node, "annotations"):
					for a in node.annotations:
						if a.name==options.colour_by_annotation and a.value not in annotation_list:
							annotation_list.append(a.value)	
			annotation_list.sort()
			colour_dictionary=get_colours_for_list(annotation_list, data_type="continuous", start_angle=240.0, end_angle=0.0)
			
			for node in t.postorder_node_iter():
				if hasattr(node, "annotations"):
					for a in node.annotations:
						if a.name==options.colour_by_annotation and a.value in colour_dictionary:
							node.edge_colour=colour_dictionary[a.value]
			
		

		
		return t


def get_leaf_nodes(t):
	leaves=[]
	for leaf in t.leaf_iter():
		leaves.append(leaf)
		
	return leaves


def get_vertical_positions_of_leaves(t):
		#Set the vertical position of the leaves
		vert_scaling_factor=float(height-20)/leaf_count
		count=1.0
		for leaf in t.leaf_iter():
			leaf.vertpos=vert_scaling_factor*count
			count+=1
			
		return
		mycolours=[(1,0,0),(0,1,0),(0,0,1)]
		
		for leaf in t.leaf_iter():
			leaf.edge_colours=[set([mycolours[randrange(0, 3)]])]
		postorder_node_list=[]
		for node in t.postorder_node_iter():
			postorder_node_list.append(node)
		preorder_node_list=[]
		for node in t.preorder_node_iter():
			preorder_node_list.append(node)
		
		dendropy.treecalc.fitch_down_pass(postorder_node_list, attr_name='edge_colours', weight_list=None, taxa_to_state_set_map=None)
		dendropy.treecalc.fitch_up_pass(preorder_node_list, attr_name='edge_colours', taxa_to_state_set_map=None)
		
#		for node in t.postorder_node_iter():
#			print node.edge_colour
#			print len(n.postorder_node_iter())
			
		draw_dendropy_tree(t, height-20, width-20, 10, 10, name_offset=5)
		
		return 
		
		sys.exit()


def deltran_parsimony_reconstruction(t, transformation="deltran"):
	
	print "Reconstructing first colour across tree using parsimony"
	
	for leaf in t.leaf_iter():
		if hasattr(leaf, 'name_colour') and len(leaf.name_colour)>0:
			r,g,b=leaf.name_colour[0]
			leaf.edge_colours=[set([(r/255,g/255,b/255)])]
		else:
			leaf.edge_colours=[set([(0.0,0.0,0.0)])]
		#print leaf.edge_colours
	postorder_node_list=[]
	for node in t.postorder_node_iter():
		postorder_node_list.append(node)
	preorder_node_list=[]
	for node in t.preorder_node_iter():
		preorder_node_list.append(node)
	
	
	reclen=dendropy.treecalc.fitch_down_pass(postorder_node_list, attr_name='edge_colours', weight_list=None, taxa_to_state_set_map=None)
	print "Reconstruction tree length =", reclen

	if t.seed_node.parent_node!= None and not hasattr(t.seed_node.parent_node,'edge_colours'):
		t.seed_node.parent_node.edge_colours=t.seed_node.edge_colours
	dendropy.treecalc.fitch_up_pass(preorder_node_list, attr_name='edge_colours', taxa_to_state_set_map=None)
	



#############################
# Function to draw the tree #
#############################

def draw_dendropy_tree(treeObject, treeheight, treewidth, xoffset, yoffset, name_offset=5):
	
#	vertical_scaling_factor=5
	def make_edge_lengths_equal():#Updated for dendropy
		for edge in treeObject.preorder_edge_iter():
			edge.length=1
	
	def get_max_and_min_branch_depth():#Updated for dendropy
			
		max_height=float("-Inf")
		min_height=float("Inf")
		for leaf in treeObject.leaf_iter():
			if leaf.distance_from_root()>max_height:
				max_height=leaf.distance_from_root()
		
		for node in treeObject.preorder_node_iter():
			try:
				if node.distance_from_root()<min_height:
					min_height=node.distance_from_root()
			except:
				continue
		
		
		if min_height>0:
			min_height=0
		
		return max_height, min_height
	
	
	def node_heights_from_edge_lengths():#Updated for dendropy
		
		max_height, min_height=get_max_and_min_branch_depth()
		
		for node in treeObject.preorder_node_iter():
			node.height=max_height-node.distance_from_root()
	
	
	def count_downstream_nodes(node):#Updated for dendropy
		count=-1
		for n in node.postorder_iter():
			count+=1

		return count
	
	
	def draw_column_label(xpox, ypos, fontsize, column_label):
		
		mylabel=Label()
		
		mylabel.setText(column_label.split(":")[0])
		
		mylabel.angle=control.metadata_column_label_angle
		mylabel.fontName=control.metadata_column_label_font#"Helvetica"
		mylabel.fontSize=control.metadata_column_label_size#fontsize
		mylabel.x=xpox
		mylabel.y=ypos
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
		#Axis.setPosition(self.track_position[0], self.track_position[1]+(self.track_height/2), self.track_length)
		d.add(mylabel)

	
	
	def draw_scale():
		if vertical_scaling_factor<5:
			linewidth=0.5
		else:
			linewidth=1.0
		
		if treeObject.is_BEAST:
			max_depth=max_branch_depth+(-1*min_branch_depth)
			branchlength=(max_depth)*horizontal_scaling_factor
			horizontalpos=xoffset
			scalestring = str(round_to_n(max_depth/10, 2))
			scalefontsize=fontsize
			if scalefontsize<10:
				scalefontsize=10
			if fontsize<10:
				gap=10
			else:
				gap=fontsize
			vertposbottom=treebase-gap
			d.add(Line(horizontalpos, vertposbottom, horizontalpos+branchlength, vertposbottom, strokeDashArray=[1, 2], strokeWidth=linewidth))
			vertpostop=treetop+gap
			d.add(Line(horizontalpos, vertpostop, horizontalpos+branchlength, vertpostop, strokeDashArray=[1, 2], strokeWidth=linewidth))
			
			if hasattr(treeObject, "max_leaf_sampling_date"):
				max_date=float(treeObject.max_leaf_sampling_date)
				min_date=max_date-max_depth
			else:
				max_date=0
				min_date=max_depth
			
			
				
			date_range=max_depth
			date_chunks=floor(date_range)
			if date_chunks>10:
				while date_chunks>10:
					date_chunks=round(date_chunks/2)
			if date_chunks<5:
				while date_chunks<5:
					date_chunks=round(date_chunks*2)
			
			
			if options.time_type=="years":
				datestring="Year"
			if options.time_type=="days":
				datestring="Day"
				if max_date>min_date:
					datestring="Date"
			
			if max_date<min_date:
				
				datestring=datestring+"s Before Last Sample Date"
				
			d.add(String(horizontalpos+(float(branchlength)/2), vertposbottom-(scalefontsize*2+2), datestring, textAnchor='middle', fontSize=scalefontsize, fontName='Helvetica'))
			
			if date_chunks<3:
				date_chunks=date_chunks*2
			for chunk in xrange(0,int(date_chunks)):
				chunklength=(chunk*(date_range/date_chunks))*horizontal_scaling_factor
				d.add(Line(horizontalpos+branchlength-chunklength, vertposbottom, horizontalpos+branchlength-chunklength, vertpostop, strokeDashArray=[1, 2], strokeWidth=linewidth/2))
				if max_date>min_date:
					scale_date=str(int(max_date-(chunk*(date_range/date_chunks))))
					if options.time_type=="days":
						datetime=xldate_as_datetime(float(scale_date))
						scale_date=str(datetime.day)+"/"+str(datetime.month)+"/"+str(datetime.year)
				else:
					scale_date=str(int(max_date+(chunk*(date_range/date_chunks))))
				d.add(String(horizontalpos+branchlength-chunklength, vertposbottom-(scalefontsize+1), scale_date, textAnchor='middle', fontSize=scalefontsize, fontName='Helvetica'))
				d.add(String(horizontalpos+branchlength-chunklength, vertpostop+1, scale_date, textAnchor='middle', fontSize=scalefontsize, fontName='Helvetica'))
			
			if max_date>min_date:
				scale_date=str(int(min_date))
				if options.time_type=="days":
					datetime=xldate_as_datetime(float(scale_date))
					scale_date=str(datetime.day)+"/"+str(datetime.month)+"/"+str(datetime.year)
			else:
				scale_date=str(int(max_date))
			
			if get_text_width('Helvetica', scalefontsize, scale_date)/1.9<chunklength/2 and get_text_width('Helvetica', scalefontsize, scale_date)/1.9<horizontalpos:
				d.add(String(horizontalpos, vertposbottom-(scalefontsize+1), scale_date, textAnchor='middle', fontSize=scalefontsize, fontName='Helvetica'))
				d.add(String(horizontalpos, vertpostop+1, scale_date, textAnchor='middle', fontSize=scalefontsize, fontName='Helvetica'))
				
#				if (chunk % 2)==0:
#					if chunk<int(date_chunks)-1:
#						start_of_box=horizontalpos+branchlength-(((chunk+1)*(date_range/date_chunks))*horizontal_scaling_factor)
#					else:
#						start_of_box=horizontalpos+branchlength-(max_depth*horizontal_scaling_factor)
#					end_of_box=horizontalpos+branchlength-chunklength
#					d.add(Rect(start_of_box, vertposbottom, end_of_box-start_of_box, vertpostop-vertposbottom, fillColor=colors.antiquewhite, strokeColor=None, strokeWidth=0))
#				else:
#					if chunk<int(date_chunks)-1:
#						start_of_box=horizontalpos+branchlength-(((chunk+1)*(date_range/date_chunks))*horizontal_scaling_factor)
#					else:
#						start_of_box=horizontalpos+branchlength-(max_depth*horizontal_scaling_factor)
#					end_of_box=horizontalpos+branchlength-chunklength
#					d.add(Rect(start_of_box, vertposbottom, end_of_box-start_of_box, vertpostop-vertposbottom, fillColor=colors.lightgrey, strokeColor=None, strokeWidth=0))
				
			if (chunk*(date_range/date_chunks))!=max_depth:
				d.add(Line(horizontalpos+branchlength-(max_depth*horizontal_scaling_factor), vertposbottom, horizontalpos+branchlength-(max_depth*horizontal_scaling_factor), vertpostop, strokeDashArray=[1, 2], strokeWidth=linewidth/2))
			
		else:
			vertpos=treebase-fontsize
			branchlength=round_to_n(max_branch_depth/10, 2)*horizontal_scaling_factor
			horizontalpos=xoffset+round_to_n(max_branch_depth/10, 2)*horizontal_scaling_factor
			scalestring = str(round_to_n(max_branch_depth/10, 2))
			scalefontsize=fontsize
			if scalefontsize<6:
				scalefontsize=6
			d.add(Line(horizontalpos, vertpos, horizontalpos+branchlength, vertpos, strokeWidth=linewidth))
			d.add(String(horizontalpos+(float(branchlength)/2), vertpos-(scalefontsize+1), scalestring, textAnchor='middle', fontSize=scalefontsize, fontName='Helvetica'))
	
		
	
	
	def set_node_vertical_positions():
	
		for node in treeObject.postorder_internal_node_iter():
			child_max=float("-Inf")
			child_min=float("Inf")
			for child in node.child_nodes():
				if child.vertpos>child_max:	
					child_max=child.vertpos
				if child.vertpos<child_min:
					child_min=child.vertpos
			
			node.vertpos=(child_max+child_min)/2
		

	def get_min_95HPD_value(hsf, mbd):
		#print hsf
		max_HPD_height=0
		for node in treeObject.preorder_node_iter():
			if node.distance_from_root()==0:
				for a  in node.annotations:
					if a.name=="height":
						root_height=a.value
			for a  in node.annotations:
				if a.name=="height_95%_HPD":
					if a.value[1]>max_HPD_height:
						max_HPD_height=a.value[1]
					# HPDmax=((max_branch_depth-a.value[0])*horizontal_scaling_factor)+xoffset-branchlength+(-1*min_branch_depth*horizontal_scaling_factor)
		#print root_height, max_HPD_height
		hsf=(hsf*root_height)/max_HPD_height
		HPD_offset=(max_HPD_height-root_height)*(root_height/mbd)



		#*horizontal_scaling_factor
		#print hsf, HPD_offset
		return hsf, HPD_offset



	
	def draw_branch(node):
		
		vertpos=node.vertpos+yoffset
		
		# if options.show_height_HPD:
		# 	horizontal_scaling_factor, HPD_extra_bit=resize_for_HPDs()
		# else:
		# 	HPD_extra_bit=0

		branchlength=node.edge_length*horizontal_scaling_factor
		
		horizontalpos=(node.distance_from_root()*horizontal_scaling_factor)+HPD_offset+xoffset-branchlength+(-1*min_branch_depth*horizontal_scaling_factor)
		
		brlabel=''
		
		if options.branch_labels!="":
			if options.branch_labels=="brlen" and node.edge_length>0:
				brlabel=str(node.edge_length).rstrip('0').rstrip('.')
			elif options.branch_labels=="support" and node.edge_length>0.0001:
				if (treeObject.schema=="nexus" and hasattr(node, "posterior") and node.posterior!=None):
					brlabel=str(float(node.posterior)).rstrip('0').rstrip('.')
				elif hasattr(node, "label") and node.label!=None:
					brlabel=str(float(node.label)).rstrip('0').rstrip('.')
			else:
				if hasattr(node, "annotations"):
					for a in node.annotations:
						if a.name==options.branch_labels:
							try:
								brlabel=str(float(a.value)).rstrip('0').rstrip('.')
							except:
								brlabel=''
		
		
		if options.tree_support!="":
			max_width=vertical_scaling_factor*0.8
			
			if options.tree_support=="support":
				
				if (treeObject.schema=="nexus" and hasattr(node, "posterior") and node.posterior!=None):
					linewidth=float(node.posterior)*max_width
				elif hasattr(node, "label") and node.label!=None:
					linewidth=(float(node.label)/100)*max_width
				else:
					linewidth=max_width
			
			else:
				if hasattr(node, "annotations"):
					linewidth=max_width
					for a in node.annotations:
						if a.name==options.tree_support:
							try:
								linewidth=float(a.value)*max_width
							except:
								linewidth=max_width
				else:
					linewidth=max_width
			
		else:
			if vertical_scaling_factor<5:
				linewidth=0.5
			else:
				linewidth=1.0
		
		if vertical_scaling_factor<5:
			vlinewidth=0.5
		else:
			vlinewidth=1.0
		
		
		
		if branchlength>0 and branchlength<vlinewidth:
			branchlength=linewidth
		
		if branchlength<0 and (-1*branchlength)<vlinewidth:
			branchlength=(-1*vlinewidth)
			
		
		# if branches have colours or hilights, find out now
		if hasattr(node, 'edge_colours') and len(node.edge_colours)==1 and len(node.edge_colours[0])==1:
			r,g,b=list(node.edge_colours[0])[0]
			branch_colour=colors.Color(float(r),float(g),float(b))
			clade_hilight=None
		elif hasattr(node, 'edge_colour'):
			branch_colour=node.edge_colour
			clade_hilight=None
		elif hasattr(node, "annotations"):
			branch_colour=colors.black
			clade_hilight=None
			for a  in node.annotations:
				if a.name=="Figtree_colour" and options.figtree:
					branch_colour=a.value
				if a.name=="Figtree_hilight" and options.figtree:
					clade_hilight=a.value
				elif a.name=="branch_colour":
					branch_colour=a.value
		else:
			branch_colour=colors.black
			clade_hilight=None
		
		#If leaves have colours, find out now
		if node.is_leaf():
			#print dir(node)
			if hasattr(node, 'name_colour') and len(node.name_colour)>0:
				name_colours=[]
				for x in xrange(0,len(node.name_colour)):
					r,g,b= node.name_colour[x]
					name_colours.append(colors.Color(float(r)/255,float(g)/255,float(b)/255))
			elif hasattr(node, "taxon"):
				name_colours=[colors.black]
				if hasattr(node.taxon, "annotations"):
					for a  in node.taxon.annotations:
						if a.name=="Figtree_colour" and options.figtree:
							name_colours=[a.value]
						if a.name=="leaf_colour":
							name_colours=[a.value]
			else:
				name_colours=[colors.black]
		
		
		#If a node is hilighted, draw the hilight box now
		if clade_hilight!=None:
			max_down_leaf_vertpos=float("-Inf")
			min_down_leaf_vertpos=float("Inf")
			#we need to find the top and bottom leaves downstream of the node to define the top and bottom of the box the box
			for downsteam_leaf in node.leaf_iter():
				if downsteam_leaf.vertpos+yoffset>max_down_leaf_vertpos:
					max_down_leaf_vertpos=downsteam_leaf.vertpos+yoffset
				if downsteam_leaf.vertpos+yoffset<min_down_leaf_vertpos:
					min_down_leaf_vertpos=downsteam_leaf.vertpos+yoffset
			d.add(Rect(horizontalpos+(branchlength/2), min_down_leaf_vertpos-(vertical_scaling_factor/2), treewidth+HPD_offset+xoffset+max_name_width+(fontsize/2)-(horizontalpos+(branchlength/2)), (max_down_leaf_vertpos-min_down_leaf_vertpos)+vertical_scaling_factor, fillColor=clade_hilight, strokeColor=None, strokeWidth=0))
		
		#if the user has chosen to show 95%HPD, draw them now:
		if options.show_height_HPD:
			if hasattr(node, "annotations"):
				HPDmax=float("-Inf")
				HPDmin=float("Inf")
				median=0.0
				for a  in node.annotations:
					if a.name=="height_95%_HPD":
						HPDmax=((max_branch_depth-a.value[0])*horizontal_scaling_factor)+HPD_offset+xoffset-branchlength+(-1*min_branch_depth*horizontal_scaling_factor)
						HPDmin=((max_branch_depth-a.value[1])*horizontal_scaling_factor)+HPD_offset+xoffset-branchlength+(-1*min_branch_depth*horizontal_scaling_factor)
					elif a.name=="height_median":
						median=a.value
				if HPDmax-HPDmin>0:
					d.add(Rect(HPDmin+branchlength, vertpos-(vertical_scaling_factor*0.4), HPDmax-HPDmin, vertical_scaling_factor*0.8, fillColor=colors.lightblue, strokeColor=None, strokeWidth=0))
		
		#draw the horizontal part of the branch
		d.add(Line(horizontalpos-(vlinewidth/2), vertpos, (horizontalpos-(vlinewidth/2))+branchlength, vertpos, strokeWidth=linewidth, strokeColor=branch_colour))
		
		#If the user has chosen to show branchlengths on branches, draw them now
		if brlabel!="":
			d.add(String((horizontalpos-(linewidth/2))+(branchlength/2), vertpos+linewidth, brlabel, textAnchor='middle', fontSize=fontsize*0.9, fillColor='black', fontName='Helvetica'))


		#draw the vertical part of the branch that goes to the parent node
		d.add(Line(horizontalpos, node.vertpos+yoffset, horizontalpos, node.parent_node.vertpos+yoffset, strokeWidth=vlinewidth, strokeColor=branch_colour))
		
			
		if node.is_leaf():
			# calculate total length of gubbins to add
			
			gubbins_length=0.0
			
			colpos=0
			
			if options.taxon_names:
				namewidth=get_text_width('Helvetica', fontsize, str(node.taxon))+name_offset
				gubbins_length += namewidth
				colpos=1
			
			for x in xrange(colpos,len(name_colours)):
				gubbins_length += block_length
				if vertical_scaling_factor>2:
					spacer=2
				else:
					spacer=vertical_scaling_factor
				if x!=0:
					gubbins_length += spacer
			
			#Add the taxon names if present
			if options.taxon_names:
				if options.aligntaxa==3:
					d.add(String(horizontalpos+branchlength+(fontsize/2), vertpos-(fontsize/3), str(node.taxon), textAnchor='start', fontSize=fontsize, fillColor=name_colours[0], fontName='Helvetica'))
					block_xpos=treewidth+xoffset+(max_name_width-gubbins_length)+(fontsize/2)+namewidth
				elif options.aligntaxa==2:
					d.add(String(treewidth+xoffset+(max_name_width-gubbins_length)+(fontsize/2), vertpos-(fontsize/3), str(node.taxon), textAnchor='start', fontSize=fontsize, fillColor=name_colours[0], fontName='Helvetica'))
					block_xpos=treewidth+xoffset+(max_name_width-gubbins_length)+(fontsize/2)+namewidth
				elif options.aligntaxa==1:
					d.add(String(treewidth+(fontsize/2)+xoffset, vertpos-(fontsize/3), str(node.taxon), textAnchor='start', fillColor=name_colours[0], fontSize=fontsize, fontName='Helvetica'))
					block_xpos=treewidth+(fontsize/2)+xoffset+(fontsize/2)+namewidth
				else:
					d.add(String(horizontalpos+branchlength+(fontsize/2), vertpos-(fontsize/3), str(node.taxon), textAnchor='start', fontSize=fontsize, fillColor=name_colours[0], fontName='Helvetica'))
					block_xpos=horizontalpos+branchlength+(fontsize/2)+namewidth
			else:
				if options.aligntaxa==2 or options.aligntaxa==3:
					block_xpos=treewidth+xoffset+(max_name_width-gubbins_length)+(fontsize/2)
				elif options.aligntaxa==1:
					block_xpos=treewidth+(fontsize/2)+xoffset+(fontsize/2)
				else:
					block_xpos=horizontalpos+branchlength+(fontsize/2)
				
			
			# draw dashed lines
			
			if options.aligntaxa==1:
				d.add(Line(horizontalpos+branchlength, vertpos, treewidth+xoffset, vertpos, strokeDashArray=[1, 2], strokeWidth=linewidth/2, strokeColor=name_colours[0]))
			elif options.aligntaxa==2:
				d.add(Line(horizontalpos+branchlength, vertpos, treewidth+xoffset+(max_name_width-gubbins_length), vertpos, strokeDashArray=[1, 2], strokeWidth=linewidth/2, strokeColor=name_colours[0]))
			elif options.aligntaxa==3:
				if options.taxon_names:
					d.add(Line(horizontalpos+branchlength+(fontsize/2)+namewidth, vertpos, block_xpos, vertpos, strokeDashArray=[1, 2], strokeWidth=linewidth/2, strokeColor=name_colours[0]))
				else:
					d.add(Line(horizontalpos+branchlength, vertpos, block_xpos, vertpos, strokeDashArray=[1, 2], strokeWidth=linewidth/2, strokeColor=name_colours[0]))
			
			
			if vertical_scaling_factor>2:
				spacer=2
			else:
				spacer=vertical_scaling_factor
			#print metadata blocks
			
			for x, name_colour in enumerate(name_colours[colpos:]):
				if block_length==0:
					break

				if options.names_as_shapes=="circle":
					d.add(Circle(cx=block_xpos+block_length, cy=vertpos, r=(block_length/2), fillColor=name_colour, strokeColor=name_colour, strokeWidth=0))
				elif options.names_as_shapes in ["square", "rectangle", "auto"]:
					d.add(Rect(block_xpos, vertpos-(vertical_scaling_factor/2), block_length, vertical_scaling_factor, fillColor=name_colour, strokeColor=name_colour, strokeWidth=0))
						
				
				block_xpos+=block_length+spacer
			
			
		
	def draw_branches():
		for node in treeObject.preorder_node_iter():
			if node.edge_length==None:
				node.edge_length=0.0
			if node!=treeObject.seed_node:
				draw_branch(node)
			
		
	
	def get_max_name_width(name_offset, fontsize):
		max_width=0.0
		for leaf in treeObject.leaf_iter():
			curwidth= get_text_width("Helvetica", fontsize, str(leaf.taxon))
			if curwidth>max_width:
				max_width=curwidth
		
		return max_width
	
	
	
#	vertical_scaling_factor=float(treeheight)/(treeObject.count_terminals(node=treeObject.root)+2)
	fontsize=vertical_scaling_factor
	if fontsize>12:
		fontsize=12
	
	if fontsize<1:
		options.taxon_names=False
	
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
		max_total_name_length=(3*float(treewidth)/4)
		
		max_total_block_length=(max_total_name_length-max_name_width)
		
		if vertical_scaling_factor>2:
			spacer=2
		else:
			spacer=vertical_scaling_factor
		max_block_length=((max_total_block_length-(spacer*((len(colour_dict)-1)+colblockstart)))/len(colour_dict))	
		
		
		if max_block_length<1:
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
		if vertical_scaling_factor>2:
			spacer=2
		else:
			spacer=vertical_scaling_factor
		if x>0:
			max_name_width+=spacer
			gubbins_length+=spacer
		max_name_width+=block_length
		gubbins_length+=block_length

	treewidth-=(max_name_width+(fontsize/2))
	

	
	max_branch_depth, min_branch_depth=get_max_and_min_branch_depth()
#	max_branch_depth+=(-1*min_branch_depth)
	
	horizontal_scaling_factor=float(treewidth)/(max_branch_depth-min_branch_depth)
	print max_branch_depth, min_branch_depth
	if options.show_height_HPD:
		horizontal_scaling_factor, HPD_offset=get_min_95HPD_value(horizontal_scaling_factor, max_branch_depth)
#	xoffset+=horizontal_scaling_factor*(-1*min_branch_depth)
	else: HPD_offset=0
	
	set_node_vertical_positions()
	
	def get_tree_base_and_top(t):
		base=float("Inf")
		top=float("-Inf")
		for leaf in t.leaf_iter():
			if leaf.vertpos<base:
				base=leaf.vertpos
			if leaf.vertpos>top:
				top=leaf.vertpos
		return base, top
	
	treebase, treetop=get_tree_base_and_top(tree)
	treebase+=yoffset
	treetop+=yoffset
	if tree.draw_scale:
		draw_scale()
	
	draw_branches()
	
	if fontsize>6:
		labelfontsize=fontsize
	else:
		labelfontsize=6
	
	if control.metadata_column_labels:
		try:
			if colour_column_names:
				if vertical_scaling_factor>2:
					spacer=2
				else:
					spacer=vertical_scaling_factor
				if options.aligntaxa==2:
					column_name_x_pos=treewidth+xoffset+(max_name_width-gubbins_length)+(fontsize/2)
					column_name_y_pos=treetop+(vertical_scaling_factor/2)
					colpos=0
					if options.taxon_names and (len(colour_column_names)==1 or colour_column_names[1]!=colour_column_names[0]):
						
						draw_column_label(treewidth+xoffset+((max_name_width-gubbins_length)/2), column_name_y_pos, labelfontsize, colour_column_names[0])
						colpos=1
					elif len(colour_column_names)>1 and colour_column_names[1]==colour_column_names[0]:
						colpos+=1
					
					for x in xrange(colpos,len(colour_column_names)):
					
						draw_column_label(column_name_x_pos+(block_length/2), column_name_y_pos, labelfontsize, colour_column_names[x])
						column_name_x_pos += block_length
						column_name_x_pos += spacer
						
		except NameError:
			pass
	
	return
	

















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
	
	def get_best_feature_name(feature, name_types=["gene", "primary_name", "systematic_id", "locus_tag", "label"]):
		
		
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
		
		
		
		
		if feature.qualifiers.has_key("border"):
			borderline=feature.qualifiers["border"][0]
		elif feature.qualifiers.has_key("border"):
			borderline=feature.qualifiers["border"][0]
		else:
			borderline=colourline

		
		if len(borderline.split())==1:
			try:
				border=translator.artemis_color(borderline)
			except StandardError:
				border=translator.artemis_color("5")
		elif len(borderline.split())==3:
			border=translator.int255_color((int(borderline.split()[0]),int(borderline.split()[1]),int(borderline.split()[2])))
		else:
			print "Can't understand colour code!"
			border=translator.artemis_color("5")
		
		
			
		locations=[]
		#get gene locations (including subfeatures)
		locations=iterate_subfeatures(feature, locations)
		
		if feature.type.lower()=="cds" and emblfile:
			new_track.add_feature(locations, fillcolour=colour, strokecolour=border, strokeweight=0.5, strand=feature.strand, arrows=int(options.arrows), label=get_best_feature_name(feature))
			#gd_feature_set.add_feature(feature, color=colour, label=0, sigil=sigiltype, arrowhead_length=0.25, locations=locations)
		elif feature.type.lower()=="cds":
			new_track.add_feature(locations, fillcolour=colour, strokecolour=border, strokeweight=0.5, strand=feature.strand, arrows=int(options.arrows), label=get_best_feature_name(feature, name_types=["label"]))
		else:
			new_track.add_feature(locations, fillcolour=colour, strokecolour=border, arrows=int(options.arrows), label=get_best_feature_name(feature, name_types=["label"]))
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
	track=add_embl_to_diagram(record, incfeatures=["i", "d", "li", "del", "snp", "misc_feature", "core", "CORE", "cds", "insertion", "deletion", "recombination", "feature", "blastn_hit", "fasta_record", "contig", "repeat_region", "variation"], emblfile=False)
	
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
		
		
		if feature.qualifiers.has_key("border"):
			borderline=feature.qualifiers["border"][0]
		elif feature.qualifiers.has_key("border"):
			borderline=feature.qualifiers["border"][0]
		else:
			borderline=colourline

		
		if len(borderline.split())==1:
			try:
				border=translator.artemis_color(borderline)
			except StandardError:
				border=translator.artemis_color("5")
		elif len(borderline.split())==3:
			border=translator.int255_color((int(borderline.split()[0]),int(borderline.split()[1]),int(borderline.split()[2])))
		else:
			print "Can't understand colour code!"
			border=translator.artemis_color("5")
		
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
				new_tracks[taxonname].add_feature(locations, fillcolour=colour, strokecolour=border, arrows=arrows)
		
		else:
			if not record.name in new_tracks:
				newtrack = Track()
				newtrack.name=record.name
				new_tracks[record.name]=newtrack	
			if feature.type.lower()=="cds":
				arrows=int(options.arrows)
			else:
				arrows=0
			new_tracks[record.name].add_feature(locations, fillcolour=colour, strokecolour=border, arrows=arrows)
			
		
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
	new_tracks=add_ordered_embl_to_diagram(record, incfeatures=["i", "d", "li", "del", "snp", "misc_feature", "core", "cds", "insertion", "deletion", "recombination", "feature", "blastn_hit", "fasta_record", "repeat_region", "variation"], emblfile=False)
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
				elif SNP and options.bcfvariants not in ["S", "I", "h"] and BASEINFO["INFO"]["AF1"]<0.8 and BASEINFO["INFO"]["AF1"]>0.2:
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
	
	
	
	def set_node_vertical_positions():
		
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
		
		if options.branch_labels!="" and treeObject.node(node).data.branchlength>0:
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
		
		mylabel.setText(column_label.split(":")[0])
		
		mylabel.angle=control.metadata_column_label_angle
		mylabel.fontName=control.metadata_column_label_font#"Helvetica"
		mylabel.fontSize=control.metadata_column_label_size#fontsize
		mylabel.x=xpox
		mylabel.y=ypos
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
	
	set_node_vertical_positions()
	
	
	treebase, treetop=recurse_subtree(treeObject.root, 0)
	treebase+=yoffset
	treetop+=yoffset
	
	if tree.draw_scale:
		draw_scale()
	
	if fontsize>6:
		labelfontsize=fontsize
	else:
		labelfontsize=6
	
	
	if control.metadata_column_labels:
		try:
			if colour_column_names:
				#print options.taxon_names, colour_column_names
				if vertical_scaling_factor>2:
					spacer=2
				else:
					spacer=vertical_scaling_factor
				if options.aligntaxa==2:
					column_name_x_pos=treewidth+xoffset+(max_name_width-gubbins_length)+(fontsize/2)
					column_name_y_pos=treetop+(vertical_scaling_factor/2)
					colpos=0
					if options.taxon_names and (len(colour_column_names)==1 or colour_column_names[1]!=colour_column_names[0]):
						
						draw_column_label(treewidth+xoffset+((max_name_width-gubbins_length)/2), column_name_y_pos, labelfontsize, colour_column_names[0])
					
						colpos=1
					for x in xrange(colpos,len(colour_column_names)):
						#print spacer, column_name_x_pos
						draw_column_label(column_name_x_pos+(block_length/2), column_name_y_pos, labelfontsize, colour_column_names[x])
						column_name_x_pos += block_length
						column_name_x_pos += spacer
		except NameError:
			pass
			
		
	
	return




class blast_result:
	def __init__(self):
		self.blast_matches=[]
		self.min_id=90.0
		self.min_length=100
		self.max_e=0.0001
		self.min_bitscore=2000
	
	
	def add_m8_match(self, m8blastline):
		
		blast_match={}
		blastwords=m8blastline.strip().split()
		blast_match['query']=blastwords[0]
		blast_match['subject']=blastwords[1]
		blast_match['percent_id']=float(blastwords[2])
		blast_match['length']=float(blastwords[3])
		blast_match['query_start']=int(blastwords[6])
		blast_match['query_end']=int(blastwords[7])
		blast_match['subject_start']=int(blastwords[8])
		blast_match['subject_end']=int(blastwords[9])
		blast_match['e']=float(blastwords[10])
		blast_match['bitscore']=float(blastwords[11])
		#print blast_match
		if blast_match['e']<self.max_e and blast_match['percent_id']>self.min_id and blast_match['length']>self.min_length and blast_match['bitscore']>self.min_bitscore:
			self.blast_matches.append(blast_match)
	
	
	def parse_m8_blast_output(self, blastfile):
		for line in open(blastfile, "rU"):
			self.add_m8_match(line)
		print 'Found', len(self.blast_matches), "blast matches matching cutoffs in file", blastfile
	
	
	def print_blast_match(self, query_x1, query_x2, subject_x1, subject_x2, query_y, subject_y):
		if query_x2>query_x1 and subject_x2<subject_x1:
			fill=colors.blue
		else:
			fill=colors.red
		stroke=None
		d.add(Polygon([query_x1, query_y, query_x2, query_y, subject_x2, subject_y, subject_x1, subject_y], fillColor=fill, strokeColor=stroke))
	
	
	def print_blast_matches(self, x0, y0, trackheight, tracklength, start_loc, end_loc):
		for blast_match in self.blast_matches:
			
			length_of_printed_genomic_region=end_loc-start_loc
			horizontal_scaling_factor=float(tracklength)/length_of_printed_genomic_region
			
			if (blast_match['query_start']>start_loc and blast_match['query_start']<end_loc) or (blast_match['query_end']>start_loc and blast_match['query_end']<end_loc) or (blast_match['query_start']<start_loc and blast_match['query_end']>end_loc) or (blast_match['query_end']<start_loc and blast_match['query_start']>end_loc):
				if (blast_match['subject_start']>start_loc and blast_match['subject_start']<end_loc) or (blast_match['subject_end']>start_loc and blast_match['subject_end']<end_loc) or (blast_match['subject_start']<start_loc and blast_match['subject_end']>end_loc) or (blast_match['subject_end']<start_loc and blast_match['subject_start']>end_loc):
					
					
					if blast_match['query_start']<start_loc:
						query_x1=x0
					if blast_match['query_start']>end_loc:
						query_x1=((end_loc-start_loc)*horizontal_scaling_factor)+x0
					elif blast_match['query_start']>start_loc and blast_match['query_start']<end_loc:
						query_x1=((blast_match['query_start']-start_loc)*horizontal_scaling_factor)+x0
					
					if blast_match['query_end']<start_loc:
						query_x2=x0
					if blast_match['query_end']>end_loc:
						query_x2=((end_loc-start_loc)*horizontal_scaling_factor)+x0
					elif blast_match['query_end']>start_loc and blast_match['query_end']<end_loc:
						query_x2=((blast_match['query_end']-start_loc)*horizontal_scaling_factor)+x0
					
					if blast_match['subject_start']<start_loc:
						subject_x1=x0
					if blast_match['subject_start']>end_loc:
						subject_x1=((end_loc-start_loc)*horizontal_scaling_factor)+x0
					elif blast_match['subject_start']>start_loc and blast_match['subject_start']<end_loc:
						subject_x1=((blast_match['subject_start']-start_loc)*horizontal_scaling_factor)+x0
					
					if blast_match['subject_end']<start_loc:
						subject_x2=x0
					if blast_match['subject_end']>end_loc:
						subject_x2=((end_loc-start_loc)*horizontal_scaling_factor)+x0
					elif blast_match['subject_end']>start_loc and blast_match['subject_end']<end_loc:
						subject_x2=((blast_match['subject_end']-start_loc)*horizontal_scaling_factor)+x0
					
#					print query_x1, query_x2, subject_x1, subject_x2, y0, y0+trackheight
#					print trackheight, horizontal_scaling_factor, length_of_printed_genomic_region, blast_match['query_start'], blast_match['query_end'], blast_match['subject_start'], blast_match['subject_end']
					
					self.print_blast_match(query_x1, query_x2, subject_x1, subject_x2, y0, y0+trackheight)





def add_blast_track(filename):
	 newtrack= Track()
	 newtrack.blast_output=blast_result()
	 newtrack.blast_output.parse_m8_blast_output(filename)
	 #newtrack.blast_output.print_blast_matches()
	 return newtrack





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
	def __init__(self, track_position=[-1,-1], track_height=1, track_length=0, track_draw_proportion=0.75, scale=False, tick_marks=True, tick_mark_number=5, tick_mark_labels=True, minor_tick_marks=True, minor_tick_mark_number=3, features=[], beginning=0, end=-1, minimum_feature_length=0.5):
	
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
		self.feature_scale=False
		self.feature_maxticks=7
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
		self.minimum_feature_length=minimum_feature_length

	
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
				if finish<start+1:
					finish=start+1
				
				if self.beginning!=0:
					if start<self.beginning and finish>self.beginning:
						start=self.beginning
				if self.end!=-1:
					if start<self.end and finish>self.end:
						finish=self.end
				start-=self.beginning
				finish-=self.beginning
				
				scaled_start=(float(start)/length)*self.track_length
				scaled_finish=(float(finish)/length)*self.track_length
				scaledlocations.append((scaled_start,scaled_finish))
				feature.scaled_feature_locations=[(start,finish)]
				
			
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
			unscaled_feature=self.features[featurenum[1]]
				
#			#if the feature is white, outline it in black so we can see it and outline features in black if selected in the options
#			if feature.strokecolour==colors.Color(1,1,1,1) or options.outline_features:
#				feature.strokecolour=colors.Color(0,0,0,1)
#			el
			if feature.strokecolour==feature.fillcolour:
				feature.strokecolour=None
				feature.strokeweight=0
			else:
				feature.strokeweight=2
			
			
			
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
					if location[1]-location[0]<self.minimum_feature_length:
						location=(location[0],location[0]+self.minimum_feature_length)
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
						
			
			
					if self.feature_scale:
						
						
						def round_to_1sf(x):
						
							return round(x, -int(floor(log10(x))))
						
						def round_up_to_1sf(x):
						
							multiplier=1
							y=x
							while y>10:
								multiplier=multiplier*10
								y=y/10
							
							return ceil(y)*multiplier
						
						
						feature_scaling_factor=(unscaled_feature.scaled_feature_locations[0][1]-unscaled_feature.scaled_feature_locations[0][0])/(feature.feature_locations[0][1]-feature.feature_locations[0][0])
#						if feature_scaling_factor*location[0]>self.end:
#							print feature_scaling_factor*location[0], feature_scaling_factor*self.track_position[0], self.beginning, self.end
						
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
						xAxis.setPosition(self.track_position[0]+location[0], yposition, location[1]-location[0])
						xAxis.configure([(1,(location[1]-location[0])*feature_scaling_factor)])
						xAxis.valueMax=(location[1]-location[0])*feature_scaling_factor
						xAxis.valueMin=0
						if self.tick_marks:
							xAxis.tickDown=((float(self.track_height)/8)*self.track_draw_proportion)
						else:
							xAxis.visibleTicks=0
						if self.minor_tick_marks:
							xAxis.visibleSubTicks=1
						if not self.tick_mark_labels:
							xAxis.visibleLabels=0
						xAxis.tickDown=1
						xAxis.tickUp=0
						xAxis.labels.boxAnchor = 'n'
						#xAxis.labels.dy=((float(self.track_height)/6)*self.track_draw_proportion)+(inch*0.02)
						#xAxis.labels.dy=(((float(self.track_height)/8)*self.track_draw_proportion)+(inch*0.02))*-1
						xAxis.labels.dy=0
						#xAxis.labels.angle=self.tick_mark_label_angle*-1
						xAxis.labels.angle=0
						xAxis.labels.fontName=self.tick_mark_label_font
						if feature.fillcolour==colors.orange:
							xAxis.labels.fillColor=colors.brown
							xAxis.labels.strokeColor=colors.brown
						elif feature.fillcolour==colors.brown:
							xAxis.labels.fillColor=colors.orange
							xAxis.labels.strokeColor=colors.orange
		
						xAxis.labels.fontSize=((float(self.track_height)/5)*self.track_draw_proportion)
						if xAxis.labels.fontSize>self.tick_mark_label_size:
							xAxis.labels.fontSize=self.tick_mark_label_size
							
						xAxis.configure([(1,((location[1]-get_text_width(xAxis.labels.fontName, xAxis.labels.fontSize, str(feature_scaling_factor*location[1])))-location[0])*feature_scaling_factor)])
		
						
						max_text_len=get_text_width(xAxis.labels.fontName, xAxis.labels.fontSize, str(feature_scaling_factor*(location[1]-location[0])))
						
						maxticks=round_to_1sf((location[1]-location[0])/max_text_len)
						tmp=(location[1]-location[0])/max_text_len
						
						if maxticks>self.feature_maxticks:
							maxticks=self.feature_maxticks
						if maxticks>0:
							tickdistance=round_up_to_1sf((feature_scaling_factor*(location[1]-location[0]))/maxticks)
							
							valueSteps=[]
							x=tickdistance
							while x<feature_scaling_factor*(((location[1]-location[0]))-(get_text_width(xAxis.labels.fontName, xAxis.labels.fontSize, str(x))/2)):
								valueSteps.append(x)
								x+=tickdistance
						else:
							tickdistance=0
							valueSteps=[]
						
						
						if len(valueSteps)>0 and xAxis.labels.fontSize>4:
						
							#print xAxis.maximumTicks
							#xAxis.maximumTicks=int(floor((location[1]-location[0])/get_text_width(xAxis.labels.fontName, xAxis.labels.fontSize, str(feature_scaling_factor*location[1]))))
							#print xAxis.maximumTicks
							#xAxis.minimumTickSpacing=get_text_width(xAxis.labels.fontName, xAxis.labels.fontSize, str(feature_scaling_factor*location[1]))
							xAxis.valueSteps=valueSteps
							xAxis.configure([(1,((location[1]-location[0])*feature_scaling_factor))])
							#if (location[1]-location[0])>get_text_width(xAxis.labels.fontName, xAxis.labels.fontSize, str(feature_scaling_factor*location[1])):
							d.add(xAxis) 
							
						
					
			
			
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
		
		
		if len(self.key_data)==0:
			return
		
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
		
		if hasattr(self, "blast_output"):
			self.blast_output.print_blast_matches(self.track_position[0], self.track_position[1]-((float(self.track_height)/2)*self.track_draw_proportion), self.track_height*self.track_draw_proportion, self.track_length, self.beginning, self.end)
		
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
		self.line_colours=[]
		for color in default_plot_colours:
			self.line_colours.append(colourconverter[color])
#		self.line_colours=[colors.red, colors.blue, colors.green, colors.black, colors.magenta, colors.cyan, colors.yellow]#default_rgb_colours
		self.line_label_colours=[]
		for color in default_plot_colours:
			self.line_label_colours.append(colourconverter[color])
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
		self.circular=False
		self.labels=[]
		self.legend=True
		self.autolegend=True
		self.legend_font="Helvetica"
		self.legend_font_size=8
		self.reorder_data=True
		self.transparency=0.60
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
		
		if len(data)>1 and self.reorder_data and self.plot_type not in ["stackedarea", "heat"] :
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
			self.line_colours=[]
			self.line_label_colours=[]
			for x in xrange(len(colours)):
				colours[x].alpha=self.transparency
			for x, y in enumerate(data_order):
#				self.labels[x]=labels[x]
				z=x
				while z>=len(colours):
					z-=len(colours)
				self.line_label_colours.append(colours[z])
				
				z=y
				while z>=len(colours):
					z-=len(colours)
				self.line_colours.append(colours[z])	
		else:
			
			if self.plot_type=="stackedarea":
				data_order=xrange(len(data)-1,-1,-1)
			else:
				data_order=xrange(len(data))
		
		
		self.data=[]
		self.xdata=[]
		
		
		for x in data_order:
			windowdata=[]
			xdata=[]
			
			
			if self.end>len(data[x]):
				end=len(data[x])
			else:
				end=self.end
				
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
					
					bardata.append(mean(region))
					barxdata.append(y)
				
				if len(bardata)>0:
					if self.circular:
						windowdata.append(mean([bardata[0], bardata[-1]]))
					else:
						windowdata.append(bardata[0])
					xdata.append(barxdata[0])
				
				for y in xrange(len(bardata)-1):
					#print y, y+1
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
			try:
				valueMin = float(options.plot_min)
			except ValueError:
				mincoverage=float(options.plot_min[:-1])
				if self.plot_type in ["heat"] and self.heat_colour!="redandblue":
					valueMin = median(remove_zeros_from_list(self.data[0]))*mincoverage
					#print median(self.data[0]), mincoverage, valueMin
				else:
					valueMin = float("Inf")
					for x in xrange(len(self.data)):
						curmin=median(remove_zeros_from_list(self.data[x]))*mincoverage
						if curmin<valueMin:
							valueMin=curmin
#					valueMin = min(map(median,self.data))*mincoverage	

		elif self.plot_type in ["heat"] and self.heat_colour!="redandblue":
			try:
				valueMin = min(self.data[0])
			except TypeError:
				print self.data
				print "Error reading plot file"
				sys.exit()
		else:
			valueMin = min(map(min,self.data))
		if options.plot_max!=float("Inf"):
			valueMax = options.plot_max
			try:
				valueMax = float(options.plot_max)
			except ValueError:
				maxcoverage=float(options.plot_max[:-1])
				if self.plot_type in ["heat"] and self.heat_colour!="redandblue":
					valueMax = median(remove_zeros_from_list(self.data[0]))*maxcoverage
					print median(self.data[0]), maxcoverage, valueMax
				else:
					valueMax = float("-Inf")
					for x in xrange(len(self.data)):
						curmax=median(remove_zeros_from_list(self.data[x]))*maxcoverage
						if curmax>valueMax:
							valueMax=curmax
#					valueMax = max(map(median,self.data))*maxcoverage

		elif self.plot_type in ["heat"] and self.heat_colour!="redandblue":
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
				end=self.beginning
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
		
		if len(data)==0 or end==self.beginning:
			return
		
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
			if self.heat_colour=="blackwhite":
				my_legend.fillColorStart=colors.Color(0,0,0)
				my_legend.fillColorEnd=colors.Color(1,1,1)
			elif self.heat_colour=="whiteblack":
				my_legend.fillColorStart=colors.Color(1,1,1)
				my_legend.fillColorEnd=colors.Color(0,0,0)
			elif self.heat_colour=="bluered":
				my_legend.fillColorStart=colors.Color(0,0,1)
				my_legend.fillColorEnd=colors.Color(1,0,0)
			elif self.heat_colour=="redblue":
				my_legend.fillColorStart=colors.Color(1,0,0)
				my_legend.fillColorEnd=colors.Color(0,0,1)
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
		if self.heat_colour=="redandblue":
			lastredvalue=data[0][0]
			try:
				lastbluevalue=data[1][0]
			except StandardError:
				DoError("Can only use redandblue heatmap when there are at least 2 columns in the plot file")
		draw_width=feature_width
		draw_start=0
		
		for i,datum in enumerate(data[0][1:]):
		
			if self.heat_colour=="redandblue":
				try:
					bluedatum=data[1][i+1]
				except StandardError:
					DoError("Can only use redandblue heatmap when there are at least 2 columns in the plot file")
				if valueMax-valueMin>0:
					redvalue=round(float(datum-valueMin)/(valueMax-valueMin),2)
					bluevalue=round(float(bluedatum-valueMin)/(valueMax-valueMin),2)
				else:
					redvalue=0.0
					bluevalue=0.0
				if redvalue>1:
					redvalue=1.0
				if bluevalue>1:
					bluevalue=1.0
				if bluevalue!=lastbluevalue or redvalue!=lastredvalue:
					#rcolour=colors.Color(1.0,1.0-lastredvalue,1.0-lastredvalue)
					#bcolour=colors.Color(1.0-lastbluevalue,1.0-lastbluevalue,1.0)
					
					
					red=max([lastredvalue, 1.0-lastbluevalue])
					blue=max([lastbluevalue, 1.0-lastredvalue])
					
					if lastredvalue<0.5 and lastbluevalue<0.5:
						green=mean([1-lastredvalue, 1-lastbluevalue])
					else:
						green=mean([1-lastredvalue, 1-lastbluevalue])
						
					colour=colors.Color(red,green,blue)
					colour=colors.Color(lastredvalue, 0, lastbluevalue)
					
					colour=colors.PCMYKColor(lastbluevalue*100, lastredvalue*100,0,0)
					
					#print lastredvalue, green, lastbluevalue, data[0][i], data[1][i], colour, (1.0-lastbluevalue),((1.0-lastredvalue)+(1.0-lastbluevalue))/2,(1.0-lastredvalue)
					if colour!=colors.Color(1,1,1):
						d.add(Rect(x+(draw_start*feature_width), y, draw_width, height, fillColor=colour, strokeColor=None, stroke=False, strokeWidth=0))

					draw_width=0.0
					draw_start=i+1
					lastbluevalue=bluevalue
					lastredvalue=redvalue		
			
			else:
				
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
					draw_start=i+1
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
		elif self.heat_colour=="redandblue":
			#colour=colors.Color((1.0-bluevalue),((1.0-redvalue)+(1.0-bluevalue))/2,(1.0-redvalue))
			#colour=colors.Color((1.0-bluevalue),min([(1.0-redvalue),(1.0-bluevalue)]),(1.0-redvalue))
			#colour=colors.Color(lastredvalue,0,lastbluevalue,1)
			colour=colors.PCMYKColor(lastbluevalue*100, lastredvalue*100,0,0)
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
#			if self.plot_type=="stackedarea":
#				legend.colorNamePairs  = []
#				for i in xrange(len(data)):
#					legend.colorNamePairs.append((self.line_label_colours[i], self.labels[len(data)-(i+1)]))
#					#lp.lines[i].strokeColor=self.line_colours[len(lp.data)-(j+1)]
#			else:
#				legend.colorNamePairs  = [(self.line_label_colours[i], self.labels[i]) for i in xrange(len(data))]
			legend.colorNamePairs  = [(self.line_label_colours[i], self.labels[i]) for i in xrange(len(data))]
			
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
		
		if self.min_yaxis!=float("Inf"):
			lp.yValueAxis.valueMin = round_to_n(self.min_yaxis, 2)
		else:
			lp.yValueAxis.valueMin = float("Inf")
			for datum in data:
				datum_min= min(map(lambda x: x[1], datum))
				if datum_min<lp.yValueAxis.valueMin:
					lp.yValueAxis.valueMin=datum_min

		if self.max_yaxis!=float("Inf"):
			lp.yValueAxis.valueMax = round_to_n(self.max_yaxis, 2)
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
			legend.colorNamePairs  = [(self.line_label_colours[i], self.labels[i]) for i in xrange(len(data))]
			
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
#		if options.plot_min!=float("Inf"):
#			bc.valueAxis.valueMin = options.plot_min
		if self.min_yaxis!=float("Inf"):
			bc.valueAxis.valueMin = round_to_n(self.min_yaxis, 2)
		else:
			bc.valueAxis.valueMin = round_to_n(min(map(min,data)), 2)
#		if options.plot_max!=float("Inf"):
#			bc.valueAxis.valueMax = options.plot_max
		if self.max_yaxis!=float("Inf"):
			bc.valueAxis.valueMax = round_to_n(self.max_yaxis, 2)
		else:
			maxsum=0.0
			for datum in data:
				maxsum+=max(datum)
			bc.valueAxis.valueMax = round_to_n(maxsum, 2)
			
		
		bc.valueAxis.tickLeft=0
		bc.valueAxis.tickRight=5
		
		print self.min_yaxis, self.max_yaxis, bc.valueAxis.valueMin, bc.valueAxis.valueMax, max(map(max,self.data))
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



###################################################################################################################
# Class to set default control options and change them depending on command line options and control file options #
###################################################################################################################

class control_options:

		
	def __init__(self):
	
		#General drawing options
		self.drawing_pagesize="A4"
		self.drawing_orientation="landscape"
		
		self.track_beginning=0
		self.track_end=-1
		self.track_fragments=1
		self.track_pages=1
		
		self.use_default_colourlist=False
		
		# track options specific for embl files
		self.embl_track_draw_proportion=0.9#proportion of the track that should be used for drawing features
		self.embl_scale=True#show scale on embl tracks. i.e. the horizontal scale line (True or False)
		self.embl_scale_position="middle"#embl tracks scale position (top, middle or bottom)
		self.embl_tick_marks=True#show tick marks on embl tracks scale (True or False)
		self.embl_tick_mark_number=5#Number of tick marks to show on embl tracks (per fragment, per page)
		self.embl_tick_mark_labels=True#Show tick mark labels on embl tracks - i.e. label scale (True or False)
		self.embl_tick_mark_label_font="Helvetica"#Font for embl track tick mark labels (Choose from ?)
		self.embl_tick_mark_label_size=8#Size of tick mark label font
		self.embl_tick_mark_label_angle=45#Angle of tick mark labels on embl tracks (clockwise from vertical)
		self.embl_minor_tick_marks=False#show minor tick marks on embl tracks scale (True or False)
		self.embl_minor_tick_mark_number=3#Number of minor tick marks to show on embl tracks (number between major tick marks)
		self.embl_draw_feature_labels=False#Draw feature labels on embl tracks (overridden by -l command line option)
		self.embl_feature_label_size=8#Font size of embl track feature labels
		self.embl_feature_label_angle=0#Angle of embl track feature labels (clockwise from vertical)
		self.embl_feature_label_font="Helvetica"#Font for embl track labels (Choose from ?)
		self.embl_greytrack=False#Show a coloured background on embl tracks (True or False)
		self.embl_grey_track_colour=colors.Color(0.25,0.25,0.25)#Colour for embl track background (R,G,B scales 0 to 1 for each)
		self.embl_grey_track_opacity_percent=10#percent opacity of background clolour of embl tracks
		self.embl_track_show_name=False#Show name on embl tracks (True or False)
		self.embl_track_name_font="Helvetica"#Font for embl track names
		self.embl_track_name_size=12#Size of embl track name font (this will be reduced if it is too big for the track)
		
		#bcf file options
		self.bcf_minimum_feature_length=1#minimum size for variants on bcf tracks (in points)
		
		#metadata column colour options
		self.metadata_colour_start_angle=240.0#angle on HSV colour circle to start metadata colour ranges (e.g. red=0, green=120, blue=240)
		self.metadata_colour_end_angle=0.0#angle on HSV colour circle to end metadata colour ranges (e.g. red=0, green=120, blue=240)
		self.metadata_colour_start_saturation=0.8#Saturation value for metadata colour ranges
		self.metadata_colour_end_saturation=0.8#Saturation value for metadata colour ranges
		self.metadata_colour_start_value=0.9#Value parameter to use for metadata colour ranges
		self.metadata_colour_end_value=0.9#Value parameter to use for metadata colour ranges
		self.metadata_colour_direction="clockwise"#Direction of metadata colour range (e.g. clockwise = blue -> green -> red -> blue, anticlockwise = blue -> red -> green -> blue)
		self.metadata_missing_colour="black"
		
		#metadata label options
		self.metadata_column_labels=True#show metadata column labels (True or False)
		self.metadata_column_label_font="Helvetica"
		self.metadata_column_label_size=10
		self.metadata_column_label_angle=45
		
		
		#plot options
		
		self.plot_tick_marks=True
		self.plot_tick_mark_number=5
		self.plot_tick_mark_labels=True
		self.plot_tick_mark_label_size=8
		self.plot_tick_mark_label_angle=45
		self.plot_scale=False
		self.plot_transparency=0.60
		
#		self.plot_colours=[]
#scale=False, tick_marks=True, tick_mark_number=5, tick_mark_labels=True, beginning=0, end=-1
		
#		self.label=""
#		self.primary_colour=colors.red
#		self.alternate_colour=colors.blue
#		#self.heat_colour="whiteblack"
#		self.heat_colour="bluered"
#		self.line_colours=[]
#		for color in default_plot_colours:
#			self.line_colours.append(colourconverter[color])
##		self.line_colours=[colors.red, colors.blue, colors.green, colors.black, colors.magenta, colors.cyan, colors.yellow]#default_rgb_colours
#		self.line_label_colours=[]
#		for color in default_plot_colours:
#			self.line_label_colours.append(colourconverter[color])
#		self.strokeweight=0.5
#		self.raw_data=[]
#		self.data=[]
#		self.xdata=[]
#		self.max_yaxis=float("-Inf")
#		self.min_yaxis=float("Inf")
#		self.window_size=1
#		self.number_of_windows=options.windows
#		self.plot_type="line"
#		self.beginning=-1
#		self.end=-1
#		self.circular=True
#		self.labels=[]
#		self.legend=True
#		self.autolegend=True
#		self.legend_font="Helvetica"
#		self.legend_font_size=8
#		self.reorder_data=True
#		self.transparency=0.60
#		self.data_order=[]
#		self.max_feature_length=0
	
	
	
	def read_control_file(self, control_file_name):
		
		def get_int(input_value, minimum=0, maximum=1):
			
			try: newint=int(input_value)
			except ValueError:
				print "illegal int in control file"
				print input_value
				sys.exit()
			
			if newint<minimum or newint>maximum:
				print "int is outside min/max range"
				print input_value
				sys.exit()
			
			return newint
				
		def get_float(input_value, minimum=0, maximum=1):
		
			try: newfloat=float(input_value)
			except ValueError:
				print "illegal float in control file"
				print input_value
				sys.exit()
			
			if newfloat<minimum or newfloat>maximum:
				print "float is outside min/max range"
				print input_value
				sys.exit()
				
			return newfloat
		
		def get_choice(input_value, choices=[]):
			
			if not input_value in choices:
				print "illegal choice in control file"
				print input_value
				print "choose from:", ', '.join(choices)
				sys.exit()
			else:
				return input_value
		
		def get_boolean(input_value):
			
			if input_value=='True':
				return True
			elif input_value=="False":
				return False
			else:
				print "illegal boolean in control file"
				print input_value
				sys.exit()
				
		variable_types={
			
#			drawing_pagesize
#			drawing_orientation
#			
#			track_beginning=0
#			track_end=-1
#			track_fragments=1
#			track_pages=1
			
			
			
			
			# track options specific for embl files
			'embl_track_draw_proportion': 'float',
			'embl_scale': 'boolean',
			'embl_scale_position': "choice",
			'embl_tick_marks': 'boolean',
			'embl_tick_mark_number': 'int',
			'embl_tick_mark_labels': 'boolean',
			'embl_tick_mark_label_font': 'choice',
			'embl_tick_mark_label_size': 'float',
			'embl_tick_mark_label_angle': 'float',
			'embl_minor_tick_marks': 'boolean',
			'embl_minor_tick_mark_number': 'int',
			'embl_draw_feature_labels': 'boolean',
			'embl_feature_label_size': 'float',
			'embl_feature_label_angle': 'float',
			'embl_feature_label_font': "choice",
			'metadata_column_label_angle': 'float',
#			embl_greytrack=False#Show a coloured background on embl tracks (True or False)
#			embl_grey_track_colour=colors.Color(0.25,0.25,0.25)#Colour for embl track background (R,G,B scales 0 to 1 for each)
#			embl_grey_track_opacity_percent=10#percent opacity of background clolour of embl tracks
#			embl_track_show_name=False#Show name on embl tracks (True or False)
#			embl_track_name_font="Helvetica"#Font for embl track names
#			embl_track_name_size=12#Size of embl track name font (this will be reduced if it is too big for the track)
#			
			#bcf file options
			'bcf_minimum_feature_length': "int",
			
			
			#plot file options
			'plot_tick_marks': 'boolean',
			'plot_tick_mark_number': 'int',
			'plot_tick_mark_labels': 'boolean',
			'plot_tick_mark_label_size': "float",
			'plot_tick_mark_label_angle': "float",
			'plot_scale': 'boolean',
			'plot_transparency': "float",
			
			
			#metadata column colour options

			# use the default colour list if possible
			'use_default_colourlist': 'boolean',
			'metadata_colour_start_angle': "float",
			'metadata_colour_end_angle': "float",
			
			'metadata_colour_start_saturation': "float",
			'metadata_colour_end_saturation': "float",
			'metadata_colour_start_value': "float",
			'metadata_colour_end_value': "float",
			'metadata_colour_direction': "choice",
			
			'metadata_column_label': 'boolean',
			'metadata_column_label_font': 'choice',
			'metadata_column_label_size': 'float',
			'metadata_column_label_angle': 'float',
			'metadata_missing_colour': 'choice'
			
			
			
			}
		
		maximums={
			
			'embl_track_draw_proportion': 1,
			'embl_tick_mark_number': 1000,
			'embl_tick_mark_label_size': 20,
			'embl_tick_mark_label_angle': 360,
			'embl_minor_tick_mark_number': 1000,
			'embl_feature_label_size': 20,
			'embl_feature_label_angle': 360,
			'bcf_minimum_feature_length': float("Inf"),
			'metadata_colour_start_angle': 360,
			'metadata_colour_end_angle': 360,
			'metadata_colour_start_saturation': 1,
			'metadata_colour_end_saturation': 1,
			'metadata_colour_start_value': 1,
			'metadata_colour_end_value': 1,
			'metadata_column_label_size': 20,
			'metadata_column_label_angle': 360,
			'plot_tick_mark_number': 100,
			'plot_tick_mark_label_size': 20,
			'plot_tick_mark_label_angle': 360,
			'plot_transparency': 1.0
			}
		
		minimums={
			
			'embl_track_draw_proportion': 0,
			'embl_tick_mark_number': 0,
			'embl_tick_mark_label_size': 1,
			'embl_tick_mark_label_angle': 0,
			'embl_minor_tick_mark_number': 0,
			'embl_feature_label_size': 1,
			'embl_feature_label_angle': 0,
			'bcf_minimum_feature_length': 0,
			'metadata_colour_start_angle': 0,
			'metadata_colour_end_angle': 0,
			'metadata_colour_start_saturation': 0,
			'metadata_colour_end_saturation': 0,
			'metadata_colour_start_value': 0,
			'metadata_colour_end_value': 0,
			'metadata_column_label_size': 1,
			'metadata_column_label_angle': 0,
			'plot_tick_mark_number': 0,
			'plot_tick_mark_label_size': 0,
			'plot_tick_mark_label_angle': 0,
			'plot_transparency': 0
			}
		
		choices={
			
			'embl_scale_position': ['top', 'middle', 'bottom'],
			'embl_tick_mark_label_font': gfont_list,
			'embl_feature_label_font': gfont_list,
			'metadata_colour_direction':["clockwise", "anticlockwise", "counterclockwise"],
			'metadata_column_label_font': gfont_list,
			'metadata_missing_colour': ["white", "black"]
			
			}
		
		for line in open(control_file_name):
			line=line.strip().split("#")[0]
			words=line.split("=")
			variable=words[0].strip()
			value='='.join(map(string.strip,words[1:]))
			
			
			if variable=="metadata_colour_saturation":
				variables=['metadata_colour_start_saturation', 'metadata_colour_end_saturation']
			elif variable=="metadata_colour_value":
				variables=['metadata_colour_start_value', 'metadata_colour_end_value']
			elif variable=="metadata_colour_angle":
				variables=['metadata_colour_start_angle', 'metadata_colour_end_angle']
			else:
				variables=[variable]
				
			for variable in variables:
				if variable in variable_types:
					if variable_types[variable]=="int":
						vars(self)[variable]=get_int(value, minimum=minimums[variable], maximum=maximums[variable])
					if variable_types[variable]=="float":
						vars(self)[variable]=get_float(value, minimum=minimums[variable], maximum=maximums[variable])
					if variable_types[variable]=="choice":
						vars(self)[variable]=get_choice(value, choices=choices[variable])
					if variable_types[variable]=="boolean":
						vars(self)[variable]=get_boolean(value)
			
				else:
					print "Invalid variable name in control file:", variable
			
			#print vars(self)
			
			
			
			




################
# Main program #
################		

if __name__ == "__main__":

	#starttime=time.clock()

	#Get command line arguments

	(options, args) = main()
	
	if options.printcolours:
		list_of_colours=[]
		for key in colourconverter:
			list_of_colours.append(key)
		list_of_colours.sort()
		for colour in list_of_colours:
			print colour
		print "See http://html-color-codes.info/color-names/"
		sys.exit()
	
	if (options.colour_by_annotation!="" and options.figtree and options.transformation in ["acctran", "deltran"]):
		print "Can only colour tree once. You have chosen to colour by Figtree colours, annotations and parsimony!"
		sys.exit()
	elif options.colour_by_annotation!="" and options.figtree:
		print "Can only colour tree once. You have chosen to colour by Figtree colours and annotations!"
		sys.exit()
	elif options.colour_by_annotation!="" and  options.transformation in ["acctran", "deltran"]:
		print "Can only colour tree once. You have chosen to colour by annotations and parsimony!"
		sys.exit()
	elif options.figtree and options.transformation in ["acctran", "deltran"]:
		print "Can only colour tree once. You have chosen to colour by Figtree colours and parsimony!"
		sys.exit()
	
	try:
		options.plot_max=float(options.plot_max)
	except ValueError:
		plotmaxvalue=options.plot_max[:-1]
		plotmaxx=options.plot_max[-1]
		try:
			tmp=float(plotmaxvalue)
		except ValueError:
			DoError("Cannot understand your value for the plot maximum. Must be a float or float followed by x/X")
		if plotmaxx not in ["x", "x"]:
			DoError("Cannot understand your value for the plot maximum. Must be a float or float followed by x/X")
		
	
	try:
		options.plot_min=float(options.plot_min)
	except ValueError:
		plotminvalue=options.plot_min[:-1]
		plotminx=options.plot_min[-1]
		try:
			tmp=float(plotminvalue)
		except ValueError:
			DoError("Cannot understand your value for the plot minimum. Must be a float or float followed by x/X")
		if plotmaxx not in ["x", "x"]:
			DoError("Cannot understand your value for the plot minimum. Must be a float or float followed by x/X")
	
	
	
	if options.names_as_shapes!="auto":
		options.taxon_names=False
	
	
	#get a list of available fonts
	gfont_list= Canvas(options.outputfile).getAvailableFonts()
	
	#Set other options or get them from the control file
	control=control_options()
	
	
	if os.path.isfile("control_file.txt"):
		control.read_control_file("control_file.txt")
	
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
	x,y=pagesize
	if options.page=="custom":
		width, height = pagesize
	elif options.orientation=="landscape":
		if x<y:
			height, width = pagesize
		else:
			width, height = pagesize
	else:
		if x>y:
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
	colour_columns=[]
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
				if cols[0]in ["","-","0"]:
					if len(cols)==1:
						cols=[]
					else:
						colour_columns.append(0)
						cols=cols[1:]
						
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
		
		
		try: lines=open(options.metadata,"rU").readlines()
		except StandardError:
			print "Could not open metadatafile:", options.metadata
			sys.exit()
		colourslist=[]
		allcolourslist=[]
		
		if len(colour_columns)>0:
			colour_column_names=[]
			for column in colour_columns:
				if column==0:
					colour_column_names.append("")
					colourslist.append([])
					continue
				elif column<1:
					print "column numbers must be greater than zero"
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
						if column==0:
							continue
						if column>len(words):
							DoError("Metadata row does not have enough columns: "+','.join(words))
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
				if colour_column!=0:
					if len(colour_column_names[x].split(":")[0])>0:
						newtrack.key_data=[["Key ("+colour_column_names[x].split(":")[0]+"):", colors.Color(0, 0, 0)]]
					else:
						newtrack.key_data=[["Key:", colors.Color(0, 0, 0)]]
				words=colour_column_names[x].split(":")
				
				if get_text_width(control.metadata_column_label_font, control.metadata_column_label_size, words[0])>max_length_col_name:
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
					elif words[1] in ["U", "u"]:
						newtrack.datatype="colour"
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
				
				
				
				s_start=control.metadata_colour_start_saturation
				s_end=control.metadata_colour_end_saturation
				v_start=control.metadata_colour_start_value
				v_end=control.metadata_colour_end_value
				start_angle=control.metadata_colour_start_angle
				end_angle=control.metadata_colour_end_angle
				direction=control.metadata_colour_direction
				
				
				if direction=="clockwise" and start_angle>end_angle:
					rotation_degrees=start_angle-end_angle
					direction_multiplier=-1
				elif direction=="clockwise":
					rotation_degrees=start_angle+(360-end_angle)
					direction_multiplier=-1
				elif direction in ["anticlockwise", "counterclockwise"] and start_angle>end_angle:
					rotation_degrees=(360-start_angle)+end_angle
					direction_multiplier=1
				else:
					rotation_degrees=end_angle-start_angle
					direction_multiplier=1
					
				use_default_colourlist=False
				
				if newtrack.datatype=="colour":
					try:
						for name in colourslist[x]:
							colour_dict[x][name]=(float(colourconverter[name].red)*255, float(colourconverter[name].green)*255, float(colourconverter[name].blue)*255)
							newtrack.key_data.append([name, colourconverter[name]])
					except StandardError:
						print "Unknown colour in", colour_column_names[x].split(":")[0], "column:", name
						sys.exit()
						print "Failed to read column colours. Treating as discrete character"
						newtrack.datatype=="discrete"
						colour_dict[x]={}
						if "" in colourslist[x]:
							colour_dict[x][""]=(0,0,0)
						if len(colour_column_names[x].split(":")[0])>0:
							newtrack.key_data=[["Key ("+colour_column_names[x].split(":")[0]+"):", colors.Color(0, 0, 0)]]
						else:
							newtrack.key_data=[["Key:", colors.Color(0, 0, 0)]]
							
				if newtrack.datatype!="colour":
					if len(colourslist[x])==1 and newtrack.datatype=="discrete":
	#					newtrack.key_data.append([colourslist[x][0], colors.Color(1, 0, 0)])
	#					colour_dict[x][colourslist[x][0]]=(255, 0, 0)
						use_default_colourlist=True
					elif len(colourslist[x])==2 and newtrack.datatype=="discrete":
	#					newtrack.key_data.append([colourslist[x][0], colors.Color(0, 0, 1)])
	#					newtrack.key_data.append([colourslist[x][1], colors.Color(1, 0, 0)])
	#					colour_dict[x][colourslist[x][0]]=(0, 0, 255)
	#					colour_dict[x][colourslist[x][1]]=(255, 0, 0)
						use_default_colourlist=True
					if newtrack.datatype=="continuous":
						
						if control.use_default_colourlist:
							use_default_colourlist=True
						if use_default_colourlist and newtrack.datamax-newtrack.datamin<=len(default_rgb_colours):
								for y, name in enumerate(colourslist[x]):
									colourcode=name-1
									colour_dict[x][name]=(float(colourconverter[default_rgb_colours[int(colourcode)]].red)*255, float(colourconverter[default_rgb_colours[int(colourcode)]].green)*255, float(colourconverter[default_rgb_colours[int(colourcode)]].blue)*255)
									newtrack.key_data.append([name, colourconverter[default_rgb_colours[int(colourcode)]]])
						else:
							for y, name in enumerate(colourslist[x]):
								value=name
								if value<newtrack.datamin:
									value=newtrack.datamin
								elif value>newtrack.datamax:
									value=newtrack.datamax
								
								if newtrack.datamax==newtrack.datamin:
									proportion=newtrack.datamin=(newtrack.datamax-(newtrack.datamax/1000000))
									
								proportion=((float(value)-newtrack.datamin)/((newtrack.datamax-newtrack.datamin)))
								
								h=(start_angle/360)+(direction_multiplier*(((proportion/360)*rotation_degrees)))
								v=v_start+(proportion*(v_end-v_start))
								s=s_start+(proportion*(s_end-s_start))
								if h<0:
									h=1.0+h
								red, green, blue = hsv_to_rgb(h,s,v)
								colour_dict[x][name]=(float(red)*255, float(green)*255, float(blue)*255)
								
							minh=(start_angle/360)+(direction_multiplier*(((0.0/360)*rotation_degrees)))
							minv=v_start+(0.0*(v_end-v_start))
							mins=s_start+(0.0*(s_end-s_start))
							minred, mingreen, minblue = hsv_to_rgb(minh,mins,minv)
							newtrack.key_data.append([newtrack.datamin, colors.Color(float(minred), float(mingreen), float(minblue))])
							
							maxh=(start_angle/360)+(direction_multiplier*(((1.0/360)*rotation_degrees)))
							maxv=v_start+(1.0*(v_end-v_start))
							maxs=s_start+(1.0*(s_end-s_start))
							maxred, maxgreen, maxblue = hsv_to_rgb(maxh,maxs,maxv)
							
							
							numbits=7
							redbit=(maxred-minred)/(numbits+1)
							bluebit=(maxblue-minblue)/(numbits+1)
							greenbit=(maxgreen-mingreen)/(numbits+1)
							bit=(newtrack.datamax-newtrack.datamin)/(numbits+1)
#							print minred, mingreen, minblue
#							print maxred, maxgreen, maxblue
#							print redbit, greenbit, bluebit
							
							for i in xrange(0,numbits):
								#print colors.Color(minred+((i+1)*redbit),mingreen+((i+1)*greenbit),minblue+((i+1)*bluebit))
								value=newtrack.datamin+((i+1)*bit)
								proportion=((float(value)-newtrack.datamin)/((newtrack.datamax-newtrack.datamin)))
								
								h=(start_angle/360)+(direction_multiplier*(((proportion/360)*rotation_degrees)))
								v=v_start+(proportion*(v_end-v_start))
								s=s_start+(proportion*(s_end-s_start))
								if h<0:
									h=1.0+h
								red, green, blue = hsv_to_rgb(h,s,v)
								newtrack.key_data.append(["#", colors.Color(red,green,blue)])
							
#							value=newtrack.datamin+((numbits)*bit)
#							proportion=((float(value)-newtrack.datamin)/((newtrack.datamax-newtrack.datamin)))
#							
#							h=(start_angle/360)+(direction_multiplier*(((proportion/360)*rotation_degrees)))
#							v=v_start+(proportion*(v_end-v_start))
#							s=s_start+(proportion*(s_end-s_start))
#							if h<0:
#								h=1.0+h
#							red, green, blue = hsv_to_rgb(h,s,v)
#							newtrack.key_data.append([">", colors.Color(red,green,blue)])
								#newtrack.key_data.append(["=", colors.Color(minred+((i+1)*redbit),mingreen+((i+1)*greenbit),minblue+((i+1)*bluebit))])
							#print colors.Color(minred+(numbits*redbit),mingreen+(numbits*greenbit),minblue+(numbits*bluebit))
							#newtrack.key_data.append([">", colors.Color(minred+(numbits*redbit),mingreen+(numbits*greenbit),minblue+(numbits*bluebit))])
							
							#newtrack.key_data.append(["==>", colors.Color(0,0,0)])
							
							newtrack.key_data.append([newtrack.datamax, colors.Color(float(maxred), float(maxgreen), float(maxblue))])
						
					elif newtrack.datatype=="discrete":
							if control.use_default_colourlist:
								use_default_colourlist=True
							if use_default_colourlist and len(colourslist[x])<=len(default_rgb_colours):
								for y, name in enumerate(colourslist[x]):
									colour_dict[x][name]=(float(colourconverter[default_rgb_colours[y]].red)*255, float(colourconverter[default_rgb_colours[y]].green)*255, float(colourconverter[default_rgb_colours[y]].blue)*255)
									newtrack.key_data.append([name, colourconverter[default_rgb_colours[y]]])
							else:
						
								for y, name in enumerate(colourslist[x]):
									proportion=(float(y)/(len(colourslist[x])-1))
									h=(start_angle/360)+(direction_multiplier*(((proportion/360)*rotation_degrees)))
									v=v_start+(proportion*(v_end-v_start))
									s=s_start+(proportion*(s_end-s_start))
									if h<0:
										h=1.0+h
									red, green, blue = hsv_to_rgb(h,s,v)
									colour_dict[x][name]=(float(red)*255, float(green)*255, float(blue)*255)
									newtrack.key_data.append([name, colors.Color(float(red), float(green), float(blue))])
									
								
					
					
				if options.show_metadata_key and not colour_column_names[x].split(":")[0] in found_keys:
					track_count+=5	
					my_tracks["metadata_key"+str(x)]=newtrack
					found_keys.append(colour_column_names[x].split(":")[0])
	
	
	if control.metadata_column_labels:
		topmargin=topmargin+(float(max_length_col_name)*sin(control.metadata_column_label_angle))
	
	for arg in args[::-1]:
		if arg.lower() in ["tree", "list"]:
			input_order.append(arg.lower())
			continue
		if arg.split('.')[-1].lower() in ["plot", "hist", "heat", "bar", "line", "graph", "area", "stackedarea","embl", "gb", "gbk", "tab", "bam", "bcf", "fas", "fa", "fasta", "mfa", "dna", "fst", "phylip", "phy", "nexus", "nxs", "blast"]:
			
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
				newtrack.minimum_feature_length=control.bcf_minimum_feature_length
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
					name='.'.join(arg.split('/')[-1].split('.')[:-1])+"_"+str(x)
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
			elif arg.split('.')[-1].lower()=="blast":
				newtrack=add_blast_track(arg)
				newtrack.scale=False
				newtrack.track_height=5
				track_count+=5
				newtrack.track_draw_proportion=1.0
				newtrack.beginning=options.beginning
				if options.end!=-1:
					newtrack.end=options.end
				name='.'.join(arg.split('/')[-1].split('.')[:-1])
				x=1
				while name in my_tracks:
					name='.'.join(arg.split('/')[-1].split('.')[:-1])+"_"+str(x)
					x+=1
				if not newtrack.name in track_names:
					track_names[newtrack.name]=[]
				input_order.append(name)
				my_tracks[name]=newtrack
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
				
				newtrack.track_draw_proportion=control.embl_track_draw_proportion
				newtrack.scale=control.embl_scale
				newtrack.scale_position=control.embl_scale_position
				newtrack.tick_marks=control.embl_tick_marks
				newtrack.tick_mark_number=control.embl_tick_mark_number
				newtrack.tick_mark_labels=control.embl_tick_mark_labels
				newtrack.tick_mark_label_font=control.embl_tick_mark_label_font
				newtrack.tick_mark_label_size=control.embl_tick_mark_label_size
				newtrack.tick_mark_label_angle=control.embl_tick_mark_label_angle
				newtrack.minor_tick_marks=control.embl_minor_tick_marks
				newtrack.minor_tick_mark_number=control.embl_minor_tick_mark_number
				newtrack.draw_feature_labels=control.embl_draw_feature_labels
				newtrack.feature_label_size=control.embl_feature_label_size
				newtrack.feature_label_angle=control.embl_feature_label_angle
				newtrack.feature_label_font=control.embl_feature_label_font
				
								
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
			
			elif arg.split('.')[-1].lower() in ["fas", "fasta", "mfa", "dna", "fst", "phylip", "phy", "nexus", "nxs", "fa"]:
				track_count+=options.emblheight
				
				newtrack = Track()
				try:
					fastarecord=read_seq_file(arg)
				except (StandardError, SimonError):
					DoError("Cannot open sequence file "+arg+" please check the format")
				shortname='.'.join(arg.split("/")[-1].split('.')[:-1])
				name=shortname
				x=1
				while name in my_tracks:
					name=shortname+"_"+str(x)
					x+=1
				if options.end!=-1:
					newtrack.end=options.end
				newtrack.beginning=options.beginning
				newtrack=add_sequence_file_to_diagram(fastarecord, name)
				newtrack.scale=True
				if arg.split('.')[-1].lower()=="mfa":
					newtrack.feature_scale=True
				newtrack.scale_position="middle"
				newtrack.name=name
				newtrack.track_height=options.emblheight
				
				newtrack.track_draw_proportion=control.embl_track_draw_proportion
				newtrack.scale=control.embl_scale
				newtrack.scale_position=control.embl_scale_position
				newtrack.tick_marks=control.embl_tick_marks
				newtrack.tick_mark_number=control.embl_tick_mark_number
				newtrack.tick_mark_labels=control.embl_tick_mark_labels
				newtrack.tick_mark_label_font=control.embl_tick_mark_label_font
				newtrack.tick_mark_label_size=control.embl_tick_mark_label_size
				newtrack.tick_mark_label_angle=control.embl_tick_mark_label_angle
				newtrack.minor_tick_marks=control.embl_minor_tick_marks
				newtrack.minor_tick_mark_number=control.embl_minor_tick_mark_number
				newtrack.draw_feature_labels=control.embl_draw_feature_labels
				newtrack.feature_label_size=control.embl_feature_label_size
				newtrack.feature_label_angle=control.embl_feature_label_angle
				newtrack.feature_label_font=control.embl_feature_label_font
				
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
						if options.labels>0:
							track_count+=options.labels
							newtrack.feature_label_angle=options.label_angle
							newtrack.feature_label_track_height=options.labels
							newtrack.draw_feature_labels=True

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
					if options.labels>0:
						track_count+=options.labels
						newtrack.feature_label_angle=options.label_angle
						newtrack.feature_label_track_height=options.labels
						newtrack.draw_feature_labels=True
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
#			d = Drawing(width, height)
			
			tree=read_dendropy_tree(options.tree)
			
#			renderPDF.drawToFile(d, options.outputfile)
#			sys.exit()
			
#			if options.midpoint:
#				tree=midpoint_root(tree)
#			
#			
#			print dir(tree)
#			print list(tree.taxon_set)
#			order=get_leaf_nodes(tree)
			
			for terminal_node in tree.leaf_iter():
				terminal=str(terminal_node.taxon)
				terminal=terminal.strip("'")
				terminal=terminal.strip('"')
				treenames.append(terminal)
				if not terminal in track_names:
					track_count+=1
				tree_name_to_node[terminal]=terminal_node
				
				if terminal in namecolours and len(namecolours[terminal])>0:
					terminal_node.name_colour=[]
					if len(colour_columns)>0:
						for x in xrange(len(colour_columns)):
							if x in namecolours[terminal]:
								namecolour=namecolours[terminal][x]
								terminal_node.name_colour.append(colour_dict[x][namecolour])
							else:
								if control.metadata_missing_colour=="white":
									terminal_node.name_colour.append((255,255,255))
								else:
									terminal_node.name_colour.append((0,0,0))
								namecolours[terminal][x]="-"
					else:
						terminal_node.name_colour.append((0,0,0))
						namecolours[terminal][0]="-"
				else:
					namecolours[terminal]={}
					terminal_node.name_colour=[]
					if len(colour_columns)>0:
						for x in xrange(len(colour_columns)):
							if control.metadata_missing_colour=="white":
								terminal_node.name_colour.append((255,255,255))
							else:
								terminal_node.name_colour.append((0,0,0))
							namecolours[terminal][x]="-"
#					else:
#						terminal_node.name_colour.append((0,0,0))
#						namecolours[terminal][0]="-"
				
				if terminal in metadatanames:
					terminal=metadatanames[terminal]
					terminal_node.taxon=terminal
#				totalbr+=tree.sum_branchlength(root=tree.root, node=terminal_node)
#			if totalbr==0:
#				
#				def make_branches_equal(node):
#					for daughter in tree.node(node).succ:
#						make_branches_equal(daughter)
#						tree.node(daughter).data.branchlength=1
#				make_branches_equal(tree.root)
		
			if options.transformation in ["acctran", "deltran"]:
				#parsimony_reconstruction(tree, namecolours, colour_dict[0], transformation=options.transformation)
				deltran_parsimony_reconstruction(tree, transformation=options.transformation)
			
				
				
				
				
				
				
				
		
		
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
					plot.data=[[minplotheight]]
				try:
					plot.max_yaxis=max(map(max,plot.data))
					plot.min_yaxis=min(map(min,plot.data))
				except ValueError:
					plot.max_yaxis=minplotheight
					plot.min_yaxis=minplotheight
					plot.data=[[minplotheight]]
			
			if options.plot_min!=float("Inf"):
				try:
					plot.min_yaxis= float(options.plot_min)				
				except ValueError:
					mincoverage=float(options.plot_min[:-1])
					if plot.plot_type in ["heat"]:
						plot.min_yaxis = median(remove_zeros_from_list(plot.data[0]))*mincoverage
					else:
						plot.min_yaxis = float("Inf")
						for x in xrange(len(plot.data)):
							curmin=median(remove_zeros_from_list(plot.data[x]))*mincoverage
							if curmin<plot.min_yaxis:
								plot.min_yaxis=curmin
								
			if options.plot_max!=float("Inf"):
				try:
					plot.max_yaxis= float(options.plot_max)						
				except ValueError:
					maxcoverage=float(options.plot_max[:-1])
					if plot.plot_type  in ["heat"]:
						plot.max_yaxis = median(remove_zeros_from_list(plot.data[0]))*maxcoverage
					else:
						plot.max_yaxis = float("-Inf")
						for x in xrange(len(plot.data)):
							curmax=median(remove_zeros_from_list(plot.data[x]))*maxcoverage
							if curmax>plot.max_yaxis:
								plot.max_yaxis=curmax
			print plot.max_yaxis, plot.min_yaxis
				
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
			tree_name_to_node[track].vertpos=margin+((track_number)*vertical_scaling_factor)+float((my_tracks[track].track_height)/2)
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
				if len(data)>0 and data[-1]>max_feature_length:
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
				draw_dendropy_tree(tree, height-(margin*2), (width-(margin*2))*left_proportion, margin, ((((options.fragments-fragment)*track_count)*vertical_scaling_factor)), maxplot_scale_text+3)
				#drawtree(tree, height-(margin*2), (width-(margin*2))*left_proportion, margin, ((((options.fragments-fragment)*track_count)*vertical_scaling_factor)), maxplot_scale_text+3)

		
		renderPDF.draw(d, c, 0, 0)
		c.showPage()
	c.save()
	print "Done"
	
	
	
