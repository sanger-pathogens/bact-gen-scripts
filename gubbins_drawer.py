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


################################
# Get the command line options #
################################


def main():
	usage = "usage: %prog [options] args"
	parser = OptionParser(usage=usage)
	
	
	group = OptionGroup(parser, "General Output Options")
	
	group.add_option("-o", "--output", action="store", dest="outputfile", help="output file name [default= %default]", type="string", metavar="FILE", default="test.pdf")
	group.add_option("-O", "--orientation", action="store", choices=['landscape', 'portrait'], dest="orientation", help="page orientation [default= %default]", type="choice", default="landscape")
	group.add_option("-p", "--pagesize", action="store", choices=['A0', 'A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'B0', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'LEGAL', 'LETTER', 'legal', 'letter'], dest="page", help="page size [default= %default]", type="choice", default="A4")
	
	parser.add_option_group(group)
	
	
	group = OptionGroup(parser, "Tree Output Options")
	
	group.add_option("-t", "--tree", action="store", dest="tree", help="tree file to align tab files to", default="")
	group.add_option("-P", "--proportion", action="store", dest="treeproportion", help="Proportion of page to take up with the tree", default=0.3, type='float')
	group.add_option("-M", "--midpoint", action="store_true", dest="midpoint", help="Midpoint root tree", default=False)
	group.add_option("-L", "--ladderise", action="store", choices=['right', 'left'], dest="ladderise", help="page size [default= %default]", type="choice", default=None)
	
	
	parser.add_option_group(group)
	
	
	group = OptionGroup(parser, "Output Options For Tracks Without Tree")
	
	group.add_option("-n", "--names", action="store_true", dest="track_names", help="Show track names", default=False)
	
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
	
	parser.add_option_group(group)
	
	return parser.parse_args()


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
		locations.append((feature.location.start.position, feature.location.end.position))
	
	
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




###########################################
# Function to add an embl file to a track #
###########################################

def add_embl_to_diagram(record, incfeatures=["CDS", "feature", "tRNA", "rRNA"], emblfile=True, name=""):
	
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
	track=add_embl_to_diagram(record, incfeatures=["i", "d", "li", "del", "snp", "misc_feature", "core", "cds", "insertion", "deletion", "recombination", "feature", "blastn_hit", "fasta_record", "contig"], emblfile=False)
	
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
	new_tracks=add_ordered_embl_to_diagram(record, incfeatures=["i", "d", "li", "del", "snp", "misc_feature", "core", "cds", "insertion", "deletion", "recombination", "feature", "blastn_hit", "fasta_record"], emblfile=False)
	return new_tracks





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
				if treeObject.node(node).data.comment==None:
					treeObject.node(node).data.comment={}
				treeObject.node(node).data.comment["vertpos"]=float(treeObject.node(daughters[0]).data.comment["vertpos"]+treeObject.node(daughters[-1]).data.comment["vertpos"])/2
				
#				treeObject.node(node).data.comment["height"]=vertical_scaling_factor
				
		
		node=treeObject.root
		get_node_vertical_position(node)
		
	
	def drawbranch(node,horizontalpos):
		
		vertpos=treeObject.node(node).data.comment["vertpos"]+yoffset
		
		horizontalpos+=xoffset
		
		branchlength=treeObject.node(node).data.branchlength*horizontal_scaling_factor
		
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
			namewidth=get_text_width('Helvetica', fontsize, treeObject.node(node).data.taxon)+name_offset
			gubbins_length += namewidth
			colpos=1
			
			for x in xrange(colpos,len(name_colours)):
				gubbins_length += block_length
				if x!=0:
					gubbins_length += vertical_scaling_factor
			
			#Add the taxon names
			d.add(String(treewidth+xoffset+(max_name_width-gubbins_length)+(fontsize/2), vertpos-(fontsize/3), treeObject.node(node).data.taxon, textAnchor='start', fontSize=fontsize, fillColor=name_colours[0], fontName='Helvetica'))
			block_xpos=treewidth+xoffset+(max_name_width-gubbins_length)+(fontsize/2)+namewidth
							
			
			# draw dashed lines
			
			d.add(Line(horizontalpos+branchlength, vertpos, treewidth+xoffset+(max_name_width-gubbins_length), vertpos, strokeDashArray=[1, 2], strokeWidth=linewidth/2, strokeColor=name_colours[0]))
			

			
			
		
	def recurse_subtree(node, horizontalpos):
		
		daughters=treeObject.node(node).succ
		
		daughterhorizontalpos=horizontalpos+(treeObject.node(node).data.branchlength*horizontal_scaling_factor)
		drawbranch(node,horizontalpos)
		for daughter in daughters:
			recurse_subtree(daughter,daughterhorizontalpos)

		
		
	
	def get_max_name_width(name_offset, fontsize):
		max_width=0.0
		for taxon in treeObject.get_terminals():
			curwidth= get_text_width("Helvetica", fontsize, treeObject.node(taxon).data.taxon)
			if curwidth>max_width:
				max_width=curwidth
		
		return max_width
	
	
	
#	vertical_scaling_factor=float(treeheight)/(treeObject.count_terminals(node=treeObject.root)+2)
	fontsize=vertical_scaling_factor
	if fontsize>12:
		fontsize=12
	
	while get_max_name_width(name_offset, fontsize)+name_offset>treewidth/3:
		fontsize-=0.2
	max_name_width=get_max_name_width(name_offset, fontsize)+name_offset
	colblockstart=1
	
	block_length=0
	
	treewidth-=(max_name_width+(fontsize/2))
	
	max_branch_depth=get_max_branch_depth()
	horizontal_scaling_factor=float(treewidth)/max_branch_depth
	
	get_node_vertical_positions()
	
	recurse_subtree(treeObject.root, 0)
	
	treebase=treeObject.node(treeObject.get_terminals()[-1]).data.comment["vertpos"]+yoffset
	
	draw_scale()
	
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
		self.name_size=10
		self.name_length=0
		self.is_key=False
		self.key_data=[]

	
	def draw_greytrack(self):
		self.grey_track_colour.alpha=float(self.grey_track_opacity_percent)/100
		d.add(Rect(self.track_position[0], self.track_position[1]-((float(self.track_height)/2)*self.track_draw_proportion), self.track_length, float(self.track_height)*self.track_draw_proportion, fillColor=self.grey_track_colour, strokeColor=self.grey_track_colour, strokeWidth = 0))
	
	
	def add_plot(self, filename, plot_type="line", fragments=1):
		
		newplot=Plot()
		
		newplot.number_of_windows=newplot.number_of_windows*fragments
		newplot.plot_type=plot_type
		datalines=newplot.read_plot_from_file(filename)
		
		newplot.beginning=self.beginning
		if self.end==-1:
			if len(newplot.raw_data)>0:
				newplot.end=len(newplot.raw_data[0])
		else:
			newplot.end=self.end

		
		newplot.read_data(datalines)
		
		self.plots.append(newplot)	
	
	
	def draw_name(self):
		d.add(String(self.track_position[0]-(self.name_length+vertical_scaling_factor), self.track_position[1]-(self.name_size/3), self.name, textAnchor='start', fontSize=self.name_size, fontName='Helvetica'))
	
	
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
			#if the feature is white, outline it in black so we can see it
			if feature.strokecolour==colors.Color(1,1,1,1):
				feature.strokecolour=colors.Color(0,0,0,1)
			
			#outline features in black if selected in the options
			if options.outline_features:
				feature.strokecolour=colors.Color(0,0,0,1)
			
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
					
					if feature.strand==0:
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
					
					
					if feature.arrows==0:
						d.add(Rect(self.track_position[0]+location[0], y, location[1]-location[0], height, fillColor=feature.fillcolour, strokeColor=feature.strokecolour, strokeWidth=feature.strokeweight))
						
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
						
				
				
			
		for join in joins:
			d.add(join)
			
		self.scaled_features=[]
			
		
	
		
	def draw_track(self):
		
		
			
		if self.greytrack and not self.is_key and len(self.plots)==0:
			self.draw_greytrack()
		
		self.draw_features()
		#self.scale=False
		if self.scale:
			self.draw_scale()
		
		if self.show_name:
			self.draw_name()
		


			
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
	
	
	pagesize=pagesizeconverter[options.page]
	
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
	
	d = Drawing(width, height)
	
	
	margin=0.5*inch
		
	
	
	
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

	
	for arg in args[::-1]:
		if arg.lower() in ["tree", "list"]:
			input_order.append(arg.lower())
			continue
		if arg.split('.')[-1].lower() in ["plot", "hist", "heat", "bar", "line", "graph", "area","embl", "gb", "tab", "bam", "fas", "fasta", "mfa", "dna", "fst", "phylip", "phy", "nexus", "nxs"]:
			
			if arg.split('.')[-1].lower() in ["embl", "gb"]:
				track_count+=options.emblheight
				
				newtrack = Track()
				try:
					emblrecord=open_annotation(arg)
				except (StandardError, SimonError):
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
				treestring=tree_to_string(tree,plain=False,ladderize=options.ladderise)
				tree=Trees.Tree(treestring, rooted=True)
				
			tree.root
			
			treeterminals=tree.get_terminals()
			totalbr=0.0
			
			
			for terminal_node in treeterminals:
				terminal=tree.node(terminal_node).data.taxon
				treenames.append(terminal)
				if not terminal in track_names:
					track_count+=1
				tree_name_to_node[terminal]=terminal_node
				tree.node(terminal_node).data.comment={}
				tree.node(terminal_node).data.comment["name_colour"]=[(0,0,0)]
				
				totalbr+=tree.sum_branchlength(root=tree.root, node=terminal_node)
			if totalbr==0:
				
				def make_branches_equal(node):
					for daughter in tree.node(node).succ:
						make_branches_equal(daughter)
						tree.node(daughter).data.branchlength=1
				make_branches_equal(tree.root)
		
			
				

	
	#from this we can work out a constant for the height of a track which takes into account the height of the page and margin sizes
	
	vertical_scaling_factor=float(height-(margin*2))/(track_count)
	
	#to make sure names can be printed in the space of a track, we can scale the name to the same size as the vertical scaling factor, but limit it to 12pt so it doesn't get crazily big
	
	name_font_size=vertical_scaling_factor
	if name_font_size>12:
		name_font_size=12

	
	left_proportion=options.treeproportion
	
	
	#Set the orders of the tracks	
	
	treetrack=0
	output_order=treenames[::-1]
	
	
	for name in input_order[::-1]:
		if not name in treenames:
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
		track_height=my_tracks[track].track_height
		
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
			my_tracks[track].grey_track_colour=colors.Color(0,0,0)
		
		track_number+=track_height
		


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
		
	#tell each track what the max feature length is
	for track in my_tracks:
		if my_tracks[track].max_feature_length<max_feature_length:
				my_tracks[track].max_feature_length=max_feature_length
	

	
	
	#We can now start to print the tracks
	
	
	beginning=0
	end=max_feature_length
	
	for track in output_order:
		
		if not track in my_tracks or (my_tracks[track].is_key and fragment!=options.fragments):
			continue
		
		my_tracks[track].beginning=beginning
		my_tracks[track].end=end
		
		
		
		my_tracks[track].track_position[1]=margin+(((my_tracks[track].track_number)*vertical_scaling_factor)+(my_tracks[track].track_height)/2)
		
		if options.greytracks:
			my_tracks[track].greytrack=True
		my_tracks[track].sort_features_by_length()
		my_tracks[track].draw_track()
	
	if options.tree!="":
		drawtree(tree, height-(margin*2), (width-(margin*2))*left_proportion, margin, 0, 5)
			
	
	
	 
	renderPDF.drawToFile(d, options.outputfile) 

	
	
	
