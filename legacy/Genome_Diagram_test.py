#!/usr/bin/env python
import string, re, numpy
import os, sys
import math
import shlex, subprocess
from Bio import SeqIO
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram._Colors import ColorTranslator
from Bio.SeqUtils import GC
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import GenBank
from optparse import OptionParser, OptionGroup
from Bio.Nexus import Trees, Nodes
sys.path.extend(map(os.path.abspath, ['/nfs/pathogen/sh16_scripts/modules/']))
from Si_general import *
from Si_SeqIO import *
import gzip


colourconverter={'aliceblue':colors.aliceblue, 'antiquewhite':colors.antiquewhite, 'aqua':colors.aqua, 'aquamarine':colors.aquamarine, 'azure':colors.azure, 'beige':colors.beige, 'bisque':colors.bisque, 'black':colors.black, 'blanchedalmond':colors.blanchedalmond, 'blue':colors.blue, 'blueviolet':colors.blueviolet, 'brown':colors.brown, 'burlywood':colors.burlywood, 'cadetblue':colors.cadetblue, 'chartreuse':colors.chartreuse, 'chocolate':colors.chocolate, 'coral':colors.coral, 'cornflower':colors.cornflower, 'cornflowerblue':colors.cornflowerblue, 'cornsilk':colors.cornsilk, 'crimson':colors.crimson, 'cyan':colors.cyan, 'darkblue':colors.darkblue, 'darkcyan':colors.darkcyan, 'darkgoldenrod':colors.darkgoldenrod, 'darkgray':colors.darkgray, 'darkgreen':colors.darkgreen, 'darkgrey':colors.darkgrey, 'darkkhaki':colors.darkkhaki, 'darkmagenta':colors.darkmagenta, 'darkolivegreen':colors.darkolivegreen, 'darkorange':colors.darkorange, 'darkorchid':colors.darkorchid, 'darkred':colors.darkred, 'darksalmon':colors.darksalmon, 'darkseagreen':colors.darkseagreen, 'darkslateblue':colors.darkslateblue, 'darkslategray':colors.darkslategray, 'darkslategrey':colors.darkslategrey, 'darkturquoise':colors.darkturquoise, 'darkviolet':colors.darkviolet, 'deeppink':colors.deeppink, 'deepskyblue':colors.deepskyblue, 'dimgray':colors.dimgray, 'dimgrey':colors.dimgrey, 'dodgerblue':colors.dodgerblue, 'fidblue':colors.fidblue, 'fidlightblue':colors.fidlightblue, 'fidred':colors.fidred, 'firebrick':colors.floralwhite, 'floralwhite':colors.floralwhite, 'forestgreen':colors.forestgreen, 'fuchsia':colors.fuchsia, 'gainsboro':colors.gainsboro, 'ghostwhite':colors.ghostwhite, 'gold':colors.gold, 'goldenrod':colors.goldenrod, 'gray':colors.gray, 'green':colors.green, 'greenyellow':colors.greenyellow, 'grey':colors.grey, 'honeydew':colors.honeydew, 'hotpink':colors.hotpink, 'indianred':colors.indianred, 'indigo':colors.indigo, 'ivory':colors.ivory, 'khaki':colors.khaki, 'lavender':colors.lavender, 'lavenderblush':colors.lavenderblush, 'lawngreen':colors.lawngreen, 'lemonchiffon':colors.lemonchiffon, 'lightblue':colors.lightblue, 'lightcoral':colors.lightcoral, 'lightcyan':colors.lightcyan, 'lightgoldenrodyellow':colors.lightgoldenrodyellow, 'lightgreen':colors.lightgreen, 'lightgrey':colors.lightgrey, 'lightpink':colors.lightpink, 'lightsalmon':colors.lightsalmon, 'lightseagreen':colors.lightseagreen, 'lightskyblue':colors.lightskyblue, 'lightslategray':colors.lightslategray, 'lightslategrey':colors.lightslategrey, 'lightsteelblue':colors.lightsteelblue, 'lightyellow':colors.lightyellow, 'lime':colors.lime, 'limegreen':colors.limegreen, 'linen':colors.linen, 'magenta':colors.magenta, 'maroon':colors.maroon, 'math':colors.math, 'mediumaquamarine':colors.mediumaquamarine, 'mediumblue':colors.mediumblue, 'mediumorchid':colors.mediumorchid, 'mediumpurple':colors.mediumpurple, 'mediumseagreen':colors.mediumseagreen, 'mediumslateblue':colors.mediumslateblue, 'mediumspringgreen':colors.mediumspringgreen, 'mediumturquoise':colors.mediumturquoise, 'mediumvioletred':colors.mediumvioletred, 'midnightblue':colors.midnightblue, 'mintcream':colors.mintcream, 'mistyrose':colors.mistyrose, 'moccasin':colors.moccasin, 'navajowhite':colors.navajowhite, 'navy':colors.navy , 'oldlace':colors.oldlace, 'olive':colors.olive, 'olivedrab':colors.olivedrab, 'orange':colors.orange, 'orangered':colors.orangered, 'orchid':colors.orchid, 'palegoldenrod':colors.palegoldenrod, 'palegreen':colors.palegreen, 'paleturquoise':colors.paleturquoise, 'palevioletred':colors.palevioletred, 'papayawhip':colors.papayawhip, 'peachpuff':colors.peachpuff, 'peru':colors.peru, 'pink':colors.pink, 'plum':colors.plum, 'powderblue':colors.powderblue, 'purple':colors.purple, 'red':colors.red, 'rosybrown':colors.rosybrown, 'royalblue':colors.royalblue, 'saddlebrown':colors.saddlebrown, 'salmon':colors.salmon, 'sandybrown':colors.sandybrown, 'seagreen':colors.seagreen, 'seashell':colors.seashell, 'sienna':colors.sienna, 'silver':colors.silver, 'skyblue':colors.skyblue, 'slateblue':colors.slateblue, 'slategray':colors.slategray, 'slategrey':colors.slategrey, 'snow':colors.snow, 'springgreen':colors.springgreen, 'steelblue':colors.steelblue, 'tan':colors.tan, 'teal':colors.teal, 'thistle':colors.thistle, 'tomato':colors.tomato, 'turquoise':colors.turquoise, 'violet':colors.violet, 'wheat':colors.wheat, 'white':colors.white, 'whitesmoke':colors.whitesmoke, 'yellow':colors.yellow, 'yellowgreen':colors.yellowgreen}


SAMTOOLS_DIR="/nfs/users/nfs_s/sh16/samtools-0.1.17/"
BCFTOOLS_DIR="/nfs/users/nfs_s/sh16/samtools-0.1.17/bcftools/"
OLD_SAMTOOLS_DIR=""
OLD_BCFTOOLS_DIR=""




##########################################
# Function to Get command line arguments #
##########################################


def main():
	usage = "usage: %prog [options] args"
	parser = OptionParser(usage=usage)
	
	group = OptionGroup(parser, "General output options")
	group.add_option("-c", "--circular", action="store_true", dest="circular", help="write output as a circular diagram [default is linear diagram]", default=False)
	group.add_option("-o", "--output", action="store", dest="outputfile", help="output file name [default= %default]", type="string", metavar="FILE", default="Genome_diagram")
	group.add_option("-f", "--fragments", action="store", dest="fragments", help="number of fragments (lines on which to split genome) for linear diagram [default= %default]", type="int", default=1, metavar="int")
	group.add_option("-p", "--pagesize", action="store", choices=["A1","A2","A3","A4","A5","LEGAL","LETTER", "estimate"], dest="page", help="page size [default= %default]", type="choice", default="A4")
	group.add_option("-O", "--orientation", action="store", choices=["landscape", "portrait"], dest="orientation", help="page orientation [default= %default]", type="choice", default="landscape")
	group.add_option("-b", "--beginning", action="store", dest="beginning", help="Start position", default=0, metavar="int", type="int")
	group.add_option("-e", "--end", action="store", dest="end", help="End position", default=-1, metavar="int", type="int")
	group.add_option("-T", "--tracksize", action="store", dest="tracksize", help="Proportion of space available to each track that should be used in drawing", default=0.75, metavar="float", type="float")
	group.add_option("-G", "--greytracks", action="store_true", dest="greytracks", help="Make tracks behind features grey", default=False)
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "Track sorting options")
	
	group.add_option("-q", "--qualifier", action="store", dest="qualifier", help="qualifier in tab file containing strain names to sort by (must be the same as in the tree provided) [e.g. note]", default="")
	group.add_option("-t", "--tree", action="store", dest="tree", help="tree file to align tab files to", default="")
	group.add_option("-n", "--taxon_list", action="store", dest="taxon_list", help="File with ordered taxa", default="")
	group.add_option("-s", "--suffix", action="store", dest="suffix", help="suffix of tab files to remove to translate file names to those in provided tree [e.g. _rec.tab]", default="")
#	parser.add_option("-r", "--reference", action="store", dest="alignment", help="multifasta file corresponding to embl files (necesary if embl file isn't in true embl format). sequence names must correspond to the embl file prefix.", default="", metavar="FILE")
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "embl file options")
	group.add_option("-E", "--emblheight", action="store", dest="emblheight", help="Relative track height of embl tracks to other tracks [default= %default]", default=2, metavar="int", type="int")
	group.add_option("-l", "--labels", action="store_true", dest="labels", help="Show CDS labels", default=False)
	group.add_option("-A", "--label_angle", action="store", dest="labelangle", help="Angle of labels in degrees from horizontal. (linear diagram only) [default= %default]", default=90, metavar="int", type="int")
	group.add_option("-S", "--label_size", action="store", dest="labelsize", help="Size of labels (CDS names and scales) [default= %default]", default=10, metavar="int", type="int")
	group.add_option("-a", "--arrows", action="store_true", dest="arrows", help="Show genes as arrows", default=False)
	group.add_option("-N", "--name_qualifier", action="store", dest="name_qualifier", help="Embl file qualifier to use for label names (e.g. gene) [Default will try to guess]", default="")
	group.add_option("-F", "--misc_features", action="store_true", dest="misc_features", help="Show misc_features in embl files [default= %default]", default=False)
	group.add_option("-g", "--gc", action="store_true", dest="gc", help="Add GC plot", default=False)
	group.add_option("-d", "--gc_dev", action="store_true", dest="gc_dev", help="Add GC deviation plot", default=False)
	group.add_option("-k", "--kyte", action="store_true", dest="kyte", help="Add Kyte-Doolittle hydrophobicity plot", default=False)
	parser.add_option_group(group)
	
	
	group = OptionGroup(parser, "Bam options")
	
	group.add_option("-u", "--base_qual_filter", action="store", dest="base_qual_filter", help="Base quality filter for bam file mapping plots [default= %default]", default=0, type="float")
	group.add_option("-U", "--mapping_qual_filter", action="store", dest="mapping_qual_filter", help="Mapping quality filter for bam file plots [default= %default]", default=0, type="float")
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "Bcf options")
	group.add_option("-B", "--bcffilter", action="store", dest="bcffilter", help="bcf filter list. Comma separated list, which must contain 10 filter values (minimum depth, mimimum strand depth, ratio of first to second base call, minimum base quality, minimum mapping quality, minimum allele frequency, strand bias p-value cutoff, base quality p-value cutoff bias, mapping bias p-value cutoff, tail distance bias p-value cutoff) [default= %default]", default="4,2,0.8,60,30,0.95,0.001, 0.0, 0.001, 0.001")
	group.add_option("-v", "--varianttype", action="store", dest="bcfvariants", choices=["a","A","i","I","s","S"], help="bcf variants to show. Letter code for types of variant to include: a/A=include all i/I:just indels, s/S:just SNPs. Lower case means exclude non-variant sites, upper case means include them [default= %default]", default="A", type="choice")
	
	parser.add_option_group(group)
	
	
	
	group = OptionGroup(parser, "Plot options")
	group.add_option("-H", "--plotheight", action="store", dest="plotheight", help="Relative track height of plot tracks to other tracks [default= %default]", default=2, metavar="int", type="int")
	
	group.add_option("-P", "--plottype", action="store", choices=["line","bar","heat"], dest="plottype", help="Type of plot to use (line, bar or heat) [default= %default]", type="choice", default="line")
	
	group.add_option("-L", "--log", action="store_true", dest="log", help="Log user plot data", default=False)
	group.add_option("-C", "--scale", action="store_true", dest="scale", help="Show scale on plots", default=False)
	group.add_option("-M", "--plot_max", action="store", dest="plotmax", help="Cutoff plots at maximum value. (linear diagram only) [default= no maximum]", default=0.4242424242, metavar="float", type="float")
	group.add_option("-m", "--plot_min", action="store", dest="plotmin", help="Cutoff plots at minimum value. (linear diagram only) [default= no minimum]", default=0.4242424242, metavar="float", type="float")
	group.add_option("-x", "--auto_max", action="store_true", dest="automax", help="Cutoff plots at maximum value calculated as the variance plus the mean", default=False)
	group.add_option("-w", "--max_window", action="store", dest="maxwindow", help="Maximum window size for plots", default=-1, metavar="int", type="int")
	group.add_option("-W", "--window", action="store", dest="windowsize", help="Window size for plots", default=-1, metavar="int", type="int")
#	group.add_option("-T", "--thin", action="store_true", dest="thin", help="Thin features to pixel width", default=False)
	parser.add_option_group(group)	
	group.add_option("-r", "--primarycolour", action="store", dest="primarycolour", choices=['aliceblue', 'antiquewhite', 'aqua', 'aquamarine', 'azure', 'beige', 'bisque', 'black', 'blanchedalmond', 'blue', 'blueviolet', 'brown', 'burlywood', 'cadetblue', 'chartreuse', 'chocolate', 'coral', 'cornflower', 'cornflowerblue', 'cornsilk', 'crimson', 'cyan', 'darkblue', 'darkcyan', 'darkgoldenrod', 'darkgray', 'darkgreen', 'darkgrey', 'darkkhaki', 'darkmagenta', 'darkolivegreen', 'darkorange', 'darkorchid', 'darkred', 'darksalmon', 'darkseagreen', 'darkslateblue', 'darkslategray', 'darkslategrey', 'darkturquoise', 'darkviolet', 'deeppink', 'deepskyblue', 'dimgray', 'dimgrey', 'dodgerblue', 'fidblue', 'fidlightblue', 'fidred', 'firebrick', 'floralwhite', 'forestgreen', 'fuchsia', 'gainsboro', 'ghostwhite', 'gold', 'goldenrod', 'gray', 'green', 'greenyellow', 'grey', 'honeydew', 'hotpink', 'indianred', 'indigo', 'ivory', 'khaki', 'lavender', 'lavenderblush', 'lawngreen', 'lemonchiffon', 'lightblue', 'lightcoral', 'lightcyan', 'lightgoldenrodyellow', 'lightgreen', 'lightgrey', 'lightpink', 'lightsalmon', 'lightseagreen', 'lightskyblue', 'lightslategray', 'lightslategrey', 'lightsteelblue', 'lightyellow', 'lime', 'limegreen', 'linen', 'magenta', 'maroon', 'mediumaquamarine', 'mediumblue', 'mediumorchid', 'mediumpurple', 'mediumseagreen', 'mediumslateblue', 'mediumspringgreen', 'mediumturquoise', 'mediumvioletred', 'midnightblue', 'mintcream', 'mistyrose', 'moccasin', 'navajowhite', 'navy', 'oldlace', 'olive', 'olivedrab', 'orange', 'orangered', 'orchid', 'palegoldenrod', 'palegreen', 'paleturquoise', 'palevioletred', 'papayawhip', 'peachpuff', 'peru', 'pink', 'plum', 'powderblue', 'purple', 'red', 'rosybrown', 'royalblue', 'saddlebrown', 'salmon', 'sandybrown', 'seagreen', 'seashell', 'sienna', 'silver', 'skyblue', 'slateblue', 'slategray', 'slategrey', 'snow', 'springgreen', 'steelblue', 'tan', 'teal', 'thistle', 'tomato', 'transparent', 'turquoise', 'violet', 'wheat', 'white', 'whitesmoke', 'yellow', 'yellowgreen'], help="Primary colour for user plots [default= %default]", default="red", type="choice")
	group.add_option("-R", "--alternatecolour", action="store", dest="alternatecolour", choices=['aliceblue', 'antiquewhite', 'aqua', 'aquamarine', 'azure', 'beige', 'bisque', 'black', 'blanchedalmond', 'blue', 'blueviolet', 'brown', 'burlywood', 'cadetblue', 'chartreuse', 'chocolate', 'coral', 'cornflower', 'cornflowerblue', 'cornsilk', 'crimson', 'cyan', 'darkblue', 'darkcyan', 'darkgoldenrod', 'darkgray', 'darkgreen', 'darkgrey', 'darkkhaki', 'darkmagenta', 'darkolivegreen', 'darkorange', 'darkorchid', 'darkred', 'darksalmon', 'darkseagreen', 'darkslateblue', 'darkslategray', 'darkslategrey', 'darkturquoise', 'darkviolet', 'deeppink', 'deepskyblue', 'dimgray', 'dimgrey', 'dodgerblue', 'fidblue', 'fidlightblue', 'fidred', 'firebrick', 'floralwhite', 'forestgreen', 'fuchsia', 'gainsboro', 'ghostwhite', 'gold', 'goldenrod', 'gray', 'green', 'greenyellow', 'grey', 'honeydew', 'hotpink', 'indianred', 'indigo', 'ivory', 'khaki', 'lavender', 'lavenderblush', 'lawngreen', 'lemonchiffon', 'lightblue', 'lightcoral', 'lightcyan', 'lightgoldenrodyellow', 'lightgreen', 'lightgrey', 'lightpink', 'lightsalmon', 'lightseagreen', 'lightskyblue', 'lightslategray', 'lightslategrey', 'lightsteelblue', 'lightyellow', 'lime', 'limegreen', 'linen', 'magenta', 'maroon', 'mediumaquamarine', 'mediumblue', 'mediumorchid', 'mediumpurple', 'mediumseagreen', 'mediumslateblue', 'mediumspringgreen', 'mediumturquoise', 'mediumvioletred', 'midnightblue', 'mintcream', 'mistyrose', 'moccasin', 'navajowhite', 'navy', 'oldlace', 'olive', 'olivedrab', 'orange', 'orangered', 'orchid', 'palegoldenrod', 'palegreen', 'paleturquoise', 'palevioletred', 'papayawhip', 'peachpuff', 'peru', 'pink', 'plum', 'powderblue', 'purple', 'red', 'rosybrown', 'royalblue', 'saddlebrown', 'salmon', 'sandybrown', 'seagreen', 'seashell', 'sienna', 'silver', 'skyblue', 'slateblue', 'slategray', 'slategrey', 'snow', 'springgreen', 'steelblue', 'tan', 'teal', 'thistle', 'tomato', 'transparent', 'turquoise', 'violet', 'wheat', 'white', 'whitesmoke', 'yellow', 'yellowgreen'], help="Alternate colour for user plots (used for bar plots and heatmaps) [default= %default]", default="blue", type="choice")
	

	
	
	return parser.parse_args()


def sort_features_by_size(features):
	newfeatures=[]
	sizes=[]
	for x, feature in enumerate(features):
		sizes.append([feature[0].location.nofuzzy_end-feature[0].location.nofuzzy_start,x])
	
	sizes.sort()
	for size in sizes[::-1]:
		newfeatures.append(features[size[1]])
	return newfeatures



	
def plot_2_moving_window(data, windowsize, log=False):
	
	movingwindowdata=[]
	mwresults=[]
	window=float(windowsize)/2
	window=int(math.floor(window))

	
	if windowsize==1:
		window=1
		for x in range(0,len(data)):
			movingwindowdata.append((x+1,data[x]))
			mwresults.append(data[x])
	else:	
		#print windowsize
		
		#print numpy.mean(data), numpy.std(data), numpy.median(data), numpy.max(data), numpy.min(data)
		
#		if options.automax:
#			#options.plotmax=numpy.mean(data)*3
#			options.plotmax=(numpy.std(data[options.beginning:options.end])*numpy.std(data[options.beginning:options.end]))+numpy.mean(data[options.beginning:options.end])
		
		for x in range(0,len(data),windowsize):
			sum=0.0
			count=0
			
			if x<beginning:
				movingwindowdata.append((x+1,"x"))
				continue
			if end!=-1 and x>end:
				movingwindowdata.append((x+1,"x"))
				continue
			
			if x<window:
				for y in data[(len(data)-(x-window)):]:
					sum=sum+y
					count+=1
				for y in data[:(x+window)]:
					sum=sum+y
					count+=1
			elif x+window>len(data):
				for y in data[x-window:]:
					sum=sum+y
					count+=1
				for y in data[:x+window-len(data)]:
					sum=sum+y
					count+=1
			else:
				for y in data[x-window:x+window]:
					sum=sum+y
					count+=1
			
			if log:
				if sum>0:
					movingwindowdata.append((x+1,math.log((float(sum)/count),1000)))
					mwresults.append(math.log((float(sum)/count),1000))
				else:
					movingwindowdata.append((x+1,math.log((1.0/windowsize),1000)))
					mwresults.append(math.log((1.0/windowsize),1000))
			else:
				result=float(sum)/count
				movingwindowdata.append((x+1,result))
				mwresults.append(result)
	
	#print numpy.mean(mwresults[options.beginning/windowsize:options.end/windowsize]), numpy.std(mwresults[options.beginning/windowsize:options.end/windowsize]), numpy.median(mwresults[options.beginning/windowsize:options.end/windowsize]), numpy.max(mwresults[options.beginning/windowsize:options.end/windowsize]), numpy.min(mwresults[options.beginning/windowsize:options.end/windowsize])
	
	if len(mwresults)>0:
		mean=numpy.mean(mwresults)
		mwmax=numpy.max(mwresults)
	else:
		mean=0
		mwmax=0
	
	if options.automax:
			#options.plotmax=numpy.mean(data)*3
			options.plotmax=(numpy.std(mwresults)*numpy.std(mwresults))+numpy.mean(mwresults[options.beginning/windowsize:options.end/windowsize])
			if options.plotmax>mwmax:
				options.plotmax=mwmax

	for x, result in enumerate(movingwindowdata):
		if result[1]=="x":
			movingwindowdata[x]=(result[0],mean)
			result=(result[0],mean)
		if options.plotmax!=0.4242424242 and result[1]>options.plotmax:
			movingwindowdata[x]=(result[0],options.plotmax)
		elif options.plotmin!=0.4242424242 and result[1]<options.plotmin:
			movingwindowdata[x]=(result[0],options.plotmin)
	return movingwindowdata



def seq_2_moving_window_gc(sequence, windowsize):
	
	movingwindowdata=[]
	window=float(windowsize)/2
	window=int(math.floor(window))
	if window<1:
		window=1
	for x in range(0,len(sequence),windowsize):
		sum=0.0
		if x==0:
			for y in sequence[(len(sequence)-window):]:
				if y.lower() in ["g", "c"]:
					sum=sum+1
			for y in sequence[:window]:
				if y.lower() in ["g", "c"]:
					sum=sum+1
		elif x+window>len(sequence):
			for y in sequence[x-window:]:
				if y.lower() in ["g", "c"]:
					sum=sum+1
			for y in sequence[:x+window-len(sequence)]:
				if y.lower() in ["g", "c"]:
					sum=sum+1
		else:
			for y in sequence[x-window:x+window]:
				if y.lower() in ["g", "c"]:
					sum=sum+1
			
		movingwindowdata.append((x+1,sum/(windowsize+1)))
		
	return movingwindowdata

def seq_2_moving_window_gc_deviation(sequence, windowsize):
	
	movingwindowdata=[]
	window=float(windowsize)/2
	window=int(math.floor(window))
	if window<1:
		window=1
	for x in range(0,len(sequence),windowsize):
		g=0.0
		c=0.0
		if x==0:
			for y in sequence[(len(sequence)-window):]:
				if y.lower()=="g":
					g=g+1
				if y.lower()=="c":
					c=c+1
			for y in sequence[:window]:
				if y.lower()=="g":
					g=g+1
				if y.lower()=="c":
					c=c+1
		elif x+window>len(sequence):
			for y in sequence[x-window:]:
				if y.lower()=="g":
					g=g+1
				if y.lower()=="c":
					c=c+1
			for y in sequence[:x+window-len(sequence)]:
				if y.lower()=="g":
					g=g+1
				if y.lower()=="c":
					c=c+1
		else:
			for y in sequence[x-window:x+window]:
				if y.lower() =="g":
					g=g+1
				if y.lower() =="c":
					c=c+1
		if (g+c)>0:
			sum=(g-c)/(g+c)
		else:
			sum=0
		movingwindowdata.append((x+1,sum/(windowsize+1)))
	
	return movingwindowdata



def add_plot_to_diagram(filename, track_number, log=False, plotdata=[]):		
	
	if len(plotdata)==0:
		print filename
		if filename.split(".")[-1]=="gz":
			handle=gzip.open(filename,"r")
			lines=handle.readlines()
			handle.close()
			plotdata=map(strip, lines)
			lines=[]
		else:
			newplot=False
			currpos=0
			for line in open(filename,"r"):
				if line.strip()[0]=='#':
					newplot=True
					print "Plot is in new format. Will show first variable only: ", line.strip().split()[1]
					continue
				else:
					currpos+=1
				if newplot:
					while float(line.strip().split()[0])>currpos:
						plotdata.append(0.0)
						currpos+=1
					plotdata.append(float(line.strip().split()[1]))
				else:
					plotdata.append(float(line.strip().split()[0]))
		
	#windowsize=len(plotdata)/500
	if end==-1:
		windowsize=(len(plotdata)-beginning)/500
	else:
		windowsize=(end-beginning)/500
	
	if windowsize<1:
		windowsize=1
	elif options.maxwindow>0 and windowsize>options.maxwindow:
		windowsize=options.maxwindow
	if options.windowsize>0:
		windowsize=options.windowsize
	print "window =", windowsize
	graphdata=plot_2_moving_window(plotdata, windowsize, log=log)
	
	testdata=[]
	for x in graphdata[options.beginning/windowsize:options.end/windowsize]:
		testdata.append(x[1])
	if len(testdata)>0:
		print "max =", numpy.max(testdata), "min =",  numpy.min(testdata)
	else:
		testdata=[0]
	

	if options.plotmax!=0.4242424242:
		graphdata.insert(0,(0,options.plotmax))
	if options.plotmin!=0.4242424242:
		graphdata.append((graphdata[-1][0]+1,options.plotmin))
	else:
		graphdata.append((graphdata[-1][0]+1,numpy.min(testdata)))
	
	
	if graphdata[-1][0]<options.beginning or (options.end!=-1 and graphdata[0][0]>options.end):
		print filename, "Contains no plot information in range specified. Skipping."
		return track_number
	
		
#	if windowsize>1:
#		graphdata.insert(0,(0,(float(graphdata[1][1])+float(graphdata[-1][1]))/2))
	
	if options.scale:
		gdt = GenomeDiagram.Track("Plot Track", scale_fontsize=int(options.labelsize), scale_fontangle=0, scale_ticks=0, scale=1, greytrack=0, greytrack_labels=0, scale_largetick_interval=len(plotdata), height=options.plotheight)
	else:
		gdt = GenomeDiagram.Track("Plot Track", scale_fontsize=int(options.labelsize), scale_fontangle=0, scale_ticks=0, scale=0, greytrack=0, greytrack_labels=0, scale_largetick_interval=len(plotdata), height=options.plotheight)
	gdgs = GenomeDiagram.GraphSet("Plot")
	gdgs.new_graph(graphdata, name="Graph Name", style=options.plottype, color=colourconverter[options.primarycolour], altcolour=colourconverter[options.alternatecolour], linewidth=1)
	gdt.add_set(gdgs)
	gd_diagram.add_track(gdt, track_number)
	track_number=track_number+1
	return track_number


def add_bam_to_diagram(filename, track_number, log=False):		
	plotdata=[]

	print filename

	samtoolssarg = shlex.split(SAMTOOLS_DIR+"samtools depth -q "+str(options.base_qual_filter)+" -Q "+str(options.mapping_qual_filter)+" "+filename)
	returnval = subprocess.Popen(samtoolssarg, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

	stdout, stderr  = returnval.communicate()

	if len(stderr)>0:
		print "Failed to open ", filename, "samtools error:", stderr
		return
	
	
	lines=stdout.split("\n")
	
	currpos=0
	for line in lines:
		currpos+=1
		if len(line.strip().split())==0:
			continue
		while float(line.strip().split()[-2])>currpos:
			#print float(line.strip().split()[-2]), currpos
			plotdata.append(0.0)
			currpos+=1
		plotdata.append(int(line.strip().split()[-1]))
	
	
	
	track_number=add_plot_to_diagram(filename, track_number, log=log, plotdata=plotdata)
	return track_number



def add_bcf_to_diagram(filename, track_number):		

	
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

		if BASEINFO["ALT"]==".":
			SNP=False
			if options.bcfvariants not in ["A", "I", "S"]:
				continue
		
		
		#Calculate the ref/alt ratios
		BASEINFO["INFO"]["DP4ratios"]={}
		if not "DP4" in BASEINFO["INFO"]:
			BASEINFO["INFO"]["DP4"]=[0,0,0,0]
			BASEINFO["INFO"]["DP4ratios"]["fref"]=0.0
			BASEINFO["INFO"]["DP4ratios"]["rref"]=0.0
			BASEINFO["INFO"]["DP4ratios"]["falt"]=0.0
			BASEINFO["INFO"]["DP4ratios"]["ralt"]=0.0
			BASEINFO["INFO"]["AF1"]=0
			BASEINFO["INFO"]["MQ"]=0
		elif "DP4" in BASEINFO["INFO"]:
			try: BASEINFO["INFO"]["DP4ratios"]["fref"]=float(BASEINFO["INFO"]["DP4"][0])/(BASEINFO["INFO"]["DP4"][0]+BASEINFO["INFO"]["DP4"][2])
			except ZeroDivisionError:
				BASEINFO["INFO"]["DP4ratios"]["fref"]=0.0
			try: BASEINFO["INFO"]["DP4ratios"]["rref"]=float(BASEINFO["INFO"]["DP4"][1])/(BASEINFO["INFO"]["DP4"][1]+BASEINFO["INFO"]["DP4"][3])
			except ZeroDivisionError:
				BASEINFO["INFO"]["DP4ratios"]["rref"]=0.0
			try: BASEINFO["INFO"]["DP4ratios"]["falt"]=float(BASEINFO["INFO"]["DP4"][2])/(BASEINFO["INFO"]["DP4"][0]+BASEINFO["INFO"]["DP4"][2])
			except ZeroDivisionError:
				BASEINFO["INFO"]["DP4ratios"]["falt"]=0.0
			try: BASEINFO["INFO"]["DP4ratios"]["ralt"]=float(BASEINFO["INFO"]["DP4"][3])/(BASEINFO["INFO"]["DP4"][1]+BASEINFO["INFO"]["DP4"][3])
			except ZeroDivisionError:
				BASEINFO["INFO"]["DP4ratios"]["ralt"]=0.0
		
	
		
		
		
		if BASEINFO["QUAL"]<options.QUAL:
			#print options.QUAL, BASEINFO["QUAL"]
			keep=False
		elif  BASEINFO["INFO"]["MQ"]<options.MQUAL:
			#print options.MQUAL, BASEINFO["INFO"]["MQ"]
			keep=False
		elif not SNP and BASEINFO["INFO"]["DP4"][0]+BASEINFO["INFO"]["DP4"][1]<options.depth:
			keep=False
		elif not SNP and BASEINFO["INFO"]["DP4"][0]<options.stranddepth:
			keep=False
		elif not SNP and BASEINFO["INFO"]["DP4"][1]<options.stranddepth:
			keep=False
		elif not SNP and BASEINFO["INFO"]["DP4ratios"]["fref"]<options.ratio:
			keep=False
		elif not SNP and BASEINFO["INFO"]["DP4ratios"]["rref"]<options.ratio:
			keep=False
		elif SNP and BASEINFO["INFO"]["DP4"][2]+BASEINFO["INFO"]["DP4"][3]<options.depth:
			keep=False
		elif SNP and BASEINFO["INFO"]["DP4"][2]<options.stranddepth:
			keep=False
		elif SNP and BASEINFO["INFO"]["DP4"][3]<options.stranddepth:
			keep=False
		elif SNP and BASEINFO["INFO"]["DP4ratios"]["falt"]<options.ratio:
			keep=False
		elif SNP and BASEINFO["INFO"]["DP4ratios"]["ralt"]<options.ratio:
			keep=False
		elif BASEINFO["ALT"]=="." and BASEINFO["INFO"]["AF1"]>(1-options.AF1):
			keep=False
		elif BASEINFO["ALT"]!="." and BASEINFO["INFO"]["AF1"]<options.AF1:
			keep=False
		elif SNP and "PV4" in BASEINFO["INFO"]:
			if BASEINFO["INFO"]["PV4"][0]<=options.strand_bias:
				keep=False
			if BASEINFO["INFO"]["PV4"][1]<=options.baseq_bias:
				keep=False
			if BASEINFO["INFO"]["PV4"][2]<=options.mapping_bias:
				keep=False
			if BASEINFO["INFO"]["PV4"][3]<=options.tail_bias:
				keep=False
			
		
		HETERO=False
		#find hetrozygous SNP calls and INDELS
		if len(BASEINFO["ALT"].split(","))>1:
			HETERO=True
			keep=False
		elif (len(BASEINFO["ALT"].split(",")[0])>1 or len(BASEINFO["REF"].split(",")[0])>1) and "INDEL" in BASEINFO['INFO']:
			INDEL=True
		elif "INDEL" in BASEINFO['INFO']:
			keep=False
		
		if options.bcfvariants in ["s", "S"]:
			INDEL=False
		
		if keep:
			if not SNP and not INDEL:
				if inmapped==False:
					colour=(230,230,230)
					inmapped=True
					start=BASEINFO["POS"]
			elif SNP:
				if inmapped:
					end=BASEINFO["POS"]
					feature = SeqFeature(FeatureLocation(start, end), strand=None)
					inmapped=False
					try: features.append((feature,colour))
					except NameError: pass #print "here"
				if INDEL and len(BASEINFO["ALT"])>len(BASEINFO["REF"]):
					colour='6'
					start=BASEINFO["POS"]
					end=BASEINFO["POS"]
					feature = SeqFeature(FeatureLocation(start, end), strand=None)
					try: features.append((feature,colour))
					except NameError: pass #print "here"
				elif INDEL and len(BASEINFO["ALT"])<len(BASEINFO["REF"]):
					colour='1'
					start=BASEINFO["POS"]
					end=BASEINFO["POS"]+(len(BASEINFO["REF"])-len(BASEINFO["ALT"]))
					feature = SeqFeature(FeatureLocation(start, end), strand=None)
					try: features.append((feature,colour))
					except NameError: pass #print "here"
				elif SNP and options.bcfvariants not in ["i", "I"] and BASEINFO["ALT"]=="A":
					colour='3'
					start=BASEINFO["POS"]
					end=BASEINFO["POS"]
					feature = SeqFeature(FeatureLocation(start, end), strand=None)
					try: features.append((feature,colour))
					except NameError: pass #print "here"
				elif SNP and options.bcfvariants not in ["i", "I"] and BASEINFO["ALT"]=="C":
					colour='2'
					start=BASEINFO["POS"]
					end=BASEINFO["POS"]
					feature = SeqFeature(FeatureLocation(start, end), strand=None)
					try: features.append((feature,colour))
					except NameError: pass #print "here"
				elif SNP and options.bcfvariants not in ["i", "I"] and BASEINFO["ALT"]=="G":
					colour='4'
					start=BASEINFO["POS"]
					end=BASEINFO["POS"]
					feature = SeqFeature(FeatureLocation(start, end), strand=None)
					try: features.append((feature,colour))
					except NameError: pass #print "here"
				elif SNP and options.bcfvariants not in ["i", "I"] and BASEINFO["ALT"]=="T":
					colour='14'
					start=BASEINFO["POS"]
					end=BASEINFO["POS"]
					feature = SeqFeature(FeatureLocation(start, end), strand=None)
					try: features.append((feature,colour))
					except NameError: pass #print "here"
		elif inmapped:
			end=BASEINFO["POS"]
			feature = SeqFeature(FeatureLocation(start, end), strand=None)
			inmapped=False
			try: features.append((feature,colour))
			except NameError: pass #print "here"
	
			
	print len(features), "features found in", filename
	
	gd_track_for_features = gd_diagram.new_track(track_number, scale=0, height=1)
	gd_feature_set = gd_track_for_features.new_set()
	
#	for x in range(len(features)-1,-1,-1):
#		print x, len(features)
#		feature=features[x]

	features=sort_features_by_size(features)

	for x, feature in enumerate(features):
		#print feature, isinstance(feature[1], str)
		if isinstance(feature[1], str):
			color = translator.artemis_color(feature[1])
		else:
			color = translator.int255_color(feature[1])
		#color = translator.artemis_color(feature[1])
		gd_feature_set.add_feature(feature[0], color=color, label=options.labels)
	
	track_number=track_number+1
	return track_number	
	

def iterate_subfeatures(feature, locations):
	if len(feature.sub_features)>0:
		for subfeature in feature.sub_features:
			locations=iterate_subfeatures(subfeature, locations)
	else:
		locations.append((feature.location.start.position, feature.location.end.position))
	
	
	return locations	


def add_embl_to_diagram(record, track_number, filetype, incfeatures=["CDS", "feature"], emblfile=True):

	incfeatures=map(string.lower,incfeatures)
	
		
	print len(record.features), "features found for", record.name
	
	if options.end==-1:
		windowsize=(len(record.seq)-options.beginning)/500
	else:
		windowsize=(options.end-options.beginning)/500
	
	if windowsize<10:
		windowsize=10
	elif options.maxwindow>0 and windowsize>options.maxwindow:
		windowsize=options.maxwindow
		
	if emblfile and options.kyte==True:
		
		
		graphdata=plot_2_moving_window(kyte_doolittle(record.seq), windowsize)
		
		graphdata[0]=(1,(float(graphdata[1][1])+float(graphdata[-1][1]))/2)
		if options.scale:
			gdt = GenomeDiagram.Track("Kyte Doolittle Track", scale_fontsize=int(options.labelsize), scale_fontangle=0, scale_ticks=0, scale=1, greytrack=0, greytrack_labels=0, scale_largetick_interval=len(record.seq), height=options.plotheight)
		else:
			gdt = GenomeDiagram.Track("Kyte Doolittle Track", scale_fontsize=int(options.labelsize), scale_fontangle=0, scale_ticks=0, scale=0, greytrack=0, greytrack_labels=0, scale_largetick_interval=len(record.seq), height=options.plotheight)
		gdgs = GenomeDiagram.GraphSet("Plot")
		gdgs.new_graph(graphdata, name="Kyte Doolittle", style=options.plottype, color=colors.darkblue, altcolour=colors.lightblue)
		gdt.add_set(gdgs)
		gd_diagram.add_track(gdt, track_number)
		track_number=track_number+1
	
	if emblfile and options.gc_dev==True:
			
		graphdata=seq_2_moving_window_gc_deviation(record.seq, windowsize)
		graphdata[0]=(1,(float(graphdata[1][1])+float(graphdata[-1][1]))/2)
		if options.scale:
			gdt = GenomeDiagram.Track("GC deviation Track", scale_fontsize=int(options.labelsize), scale_fontangle=0, scale_ticks=0, scale=1, greytrack=0, greytrack_labels=0, scale_largetick_interval=len(record.seq), height=options.plotheight)
		else:
			gdt = GenomeDiagram.Track("GC deviation Track", scale_fontsize=int(options.labelsize), scale_fontangle=0, scale_ticks=0, scale=0, greytrack=0, greytrack_labels=0, scale_largetick_interval=len(record.seq), height=options.plotheight)
		gdgs = GenomeDiagram.GraphSet("Plot")
		gdgs.new_graph(graphdata, name="GC deviation", style=options.plottype, color=colourconverter[options.primarycolour], altcolour=colourconverter[options.alternatecolour])#color=colors.violet, altcolour=colors.purple)
		gdt.add_set(gdgs)
		gd_diagram.add_track(gdt, track_number)
		track_number=track_number+1
	
	if emblfile and options.gc==True:
		
		
		graphdata=seq_2_moving_window_gc(record.seq, windowsize)
		graphdata[0]=(1,(float(graphdata[1][1])+float(graphdata[-1][1]))/2)
		
		if options.scale:
			gdt = GenomeDiagram.Track("GC content Track", scale_fontsize=int(options.labelsize), scale_fontangle=0, scale_ticks=0, scale=1, greytrack=0, greytrack_labels=0, scale_largetick_interval=len(graphdata), height=options.plotheight)
		else:
			gdt = GenomeDiagram.Track("GC content Track", scale_fontsize=int(options.labelsize), scale_fontangle=0, scale_ticks=0, scale=0, greytrack=0, greytrack_labels=0, scale_largetick_interval=len(graphdata), height=options.plotheight)
			
		gdgs = GenomeDiagram.GraphSet("Plot")
		gdgs.new_graph(graphdata, name="GC content", style=options.plottype, color=colors.lightgreen, altcolour=colors.darkgreen)
		gdt.add_set(gdgs)
		gd_diagram.add_track(gdt, track_number)
		track_number=track_number+1
	
	
	if len(record.seq)>500000:
		scale_largetick_interval=int(round((len(record.seq)/10),-5))
		scale_smalltick_interval=int(round((len(record.seq)/10),-5)/5)
	else:
		scale_largetick_interval=len(record.seq)
		scale_smalltick_interval=len(record.seq)/5
	
	if emblfile:
		gd_track_for_features = gd_diagram.new_track(track_number, scale_fontsize=int(options.labelsize), scale_largetick_interval=scale_largetick_interval, scale_smalltick_interval=scale_smalltick_interval, greytrack=greytracks, greytrack_labels=0, name="Annotated Features", height=options.emblheight)
	else:
		if track_number %2:
			gd_track_for_features = gd_diagram.new_track(track_number, scale=0, height=1, greytrack=greytracks, greytrack_labels=0)
		else:
			gd_track_for_features = gd_diagram.new_track(track_number, scale=0, height=1)
	gd_feature_set = gd_track_for_features.new_set()
	
	
	if options.misc_features:
		incfeatures.append("misc_feature")
	
	for x, feature in enumerate(record.features):
		#print feature.type, feature.location.nofuzzy_start, feature.location.nofuzzy_end, options.beginning, options.end
		if feature.type.lower() not in  incfeatures or feature.location.nofuzzy_end<options.beginning or (feature.location.nofuzzy_start>options.end and options.end!=-1):#,"tRNA","repeat_unit"] :
			#Exclude this feature
			#print "here"
			continue
		
		if feature.qualifiers.has_key("colour"):
			colourline=feature.qualifiers["colour"][0]
		elif feature.qualifiers.has_key("color"):
			colourline=feature.qualifiers["color"][0]
		else:
			colourline = translator.artemis_color("5")
			
		if len(colourline.split())==1:
			colour=translator.artemis_color(colourline)
		elif len(colourline.split())==3:
			colour=translator.int255_color((int(colourline.split()[0]),int(colourline.split()[1]),int(colourline.split()[2])))
		else:
			print "Can't understand colour code!"
			sys.exit()
		
			
		locations=[]
		#get gene locations (including subfeatures)
		locations=iterate_subfeatures(feature, locations)
		if feature.type=="CDS":	
			gd_feature_set.add_feature(feature, color=colour, label=0, sigil=sigiltype, arrowhead_length=0.25, locations=locations)
		else:
			gd_feature_set.add_feature(feature, color=colour, label=0, strand=0, locations=locations)
	
	track_number=track_number+1
	if emblfile and options.labels==True:# and options.circular==True:
		gd_track_for_features = gd_diagram.new_track(track_number, scale_ticks=0, scale=0, greytrack=0, greytrack_labels=0, name="Labels")
		gd_feature_set = gd_track_for_features.new_set()
		
		for x, feature in enumerate(record.features):
			
			if feature.type.lower() not in incfeatures or feature.location.nofuzzy_end<options.beginning or (feature.location.nofuzzy_start>options.end and options.end!=-1):#,"tRNA","repeat_unit"] :
				#Exclude this feature
				continue
			
			if options.beginning!=0:
			
				if feature.location.nofuzzy_start<options.beginning and feature.location.nofuzzy_end>options.beginning:
					feature.location=FeatureLocation(options.beginning, feature.location.nofuzzy_end)
			
			if options.end!=-1:
				if feature.location.nofuzzy_start<options.end and feature.location.nofuzzy_end>options.end:
					feature.location=FeatureLocation(feature.location.nofuzzy_start, options.end)
					

			#If the user has asked for a specific qualifier to be used as the name, find this qualifier and set it to name
			if options.name_qualifier!="":
				name=''
				for key in feature.qualifiers:
					if key==options.name_qualifier:
						name=feature.qualifiers[key][0]
			
				gd_feature_set.add_feature(feature, color=translator.int255_color((254,255,255)), label=options.labels, label_size=int(options.labelsize), label_angle=int(options.labelangle)+180, sigil=sigiltype, arrowhead_length=0.25, label_position=label_position, strand=-1, name=name)
			else:
				#If the user hasn't specified a specific qualifier to use as a name, use the default list
				gd_feature_set.add_feature(feature, color=translator.int255_color((254,255,255)), label=options.labels, label_size=int(options.labelsize), label_angle=int(options.labelangle)+180, sigil=sigiltype, arrowhead_length=0.25, label_position=label_position, strand=-1)
		track_number=track_number+1
	return track_number
	









def add_ordered_embl_to_diagram(record, track_number, filetype, nameorder, incfeatures=["CDS", "feature"], emblfile=True):

	incfeatures=map(string.lower,incfeatures)
		
	print len(record.features), "features found for", record.name
	
	if options.end==-1:
		windowsize=(len(record.seq)-options.beginning)/500
	else:
		windowsize=(options.end-options.beginning)/500
	
	if windowsize<10:
		windowsize=10
	elif options.maxwindow>0 and windowsize>options.maxwindow:
		windowsize=options.maxwindow
	
	
	if len(record.seq)>500000:
		scale_largetick_interval=int(round((len(record.seq)/10),-5))
		scale_smalltick_interval=int(round((len(record.seq)/10),-5)/5)
	else:
		scale_largetick_interval=len(record.seq)
		scale_smalltick_interval=len(record.seq)/5
	
	
	for featurename in nameorder[::-1]:
	
		if emblfile:
			gd_track_for_features = gd_diagram.new_track(track_number, scale_fontsize=int(options.labelsize), scale_largetick_interval=scale_largetick_interval, scale_smalltick_interval=scale_smalltick_interval, greytrack=greytracks, greytrack_labels=0, name="Annotated Features")
		else:
			if track_number %2:
				gd_track_for_features = gd_diagram.new_track(track_number, scale=0, height=1, greytrack=greytracks, greytrack_labels=0)
			else:
				gd_track_for_features = gd_diagram.new_track(track_number, scale=0, height=1)
		gd_feature_set = gd_track_for_features.new_set()
		
		
		if options.misc_features:
			incfeatures.append("misc_feature")
		
		for x, feature in enumerate(record.features):
			if feature.type.lower() not in incfeatures or feature.location.nofuzzy_end<options.beginning or (feature.location.nofuzzy_start>options.end and options.end!=-1) or options.qualifier not in feature.qualifiers or not featurename in feature.qualifiers[options.qualifier]:#,"tRNA","repeat_unit"] :
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
				colourline = translator.artemis_color("5")
			if len(colourline.split())==1:
				colour=translator.artemis_color(colourline)
			elif len(colourline.split())==3:
				colour=translator.int255_color((int(colourline.split()[0]),int(colourline.split()[1]),int(colourline.split()[2])))
			else:
				print "Can't understand colour code!"
				sys.exit()
				
			locations=[]
			#get gene locations (including subfeatures)
			locations=iterate_subfeatures(feature, locations)
			if feature.type=="CDS":	
				gd_feature_set.add_feature(feature, color=colour, label=0, sigil=sigiltype, arrowhead_length=0.25, locations=locations)
			else:
				gd_feature_set.add_feature(feature, color=colour, label=0, strand=0, locations=locations)
		
		track_number=track_number+1
		if emblfile and options.labels==True:# and options.circular==True:
			gd_track_for_features = gd_diagram.new_track(track_number, scale_ticks=0, scale=0, greytrack=0, greytrack_labels=0, name="Labels")
			gd_feature_set = gd_track_for_features.new_set()
			
			for x, feature in enumerate(record.features):
				
				if feature.type.lower() not in incfeatures or feature.location.nofuzzy_end<options.beginning or (feature.location.nofuzzy_start>options.end and options.end!=-1):#,"tRNA","repeat_unit"] :
					#Exclude this feature
					continue
	
				#If the user has asked for a specific qualifier to be used as the name, find this qualifier and set it to name
				if options.name_qualifier!="":
					name=''
					for key in feature.qualifiers:
						if key==options.name_qualifier:
							name=feature.qualifiers[key][0]
				
					gd_feature_set.add_feature(feature, color=translator.int255_color((254,255,255)), label=options.labels, label_size=int(options.labelsize), label_angle=int(options.labelangle)+180, sigil=sigiltype, arrowhead_length=0.25, label_position=label_position, strand=-1, name=name)
				else:
					#If the user hasn't specified a specific qualifier to use as a name, use the default list
					gd_feature_set.add_feature(feature, color=translator.int255_color((254,255,255)), label=options.labels, label_size=int(options.labelsize), label_angle=int(options.labelangle)+180, sigil=sigiltype, arrowhead_length=0.25, label_position=label_position, strand=-1)
			track_number=track_number+1
	return track_number











		
def add_tab_to_diagram(filename, track_number):
	
	features=[]
	
	if filename.split(".")[-1]=="gz":
		print "need to add gzip functionality to tab reader"
		return
	
	
	record=tab_parser(open(filename,"r"))
	record.name=filename
	track_number=add_embl_to_diagram(record, track_number, filetype, incfeatures=["i", "d", "li", "del", "snp", "misc_feature", "core", "cds", "insertion", "deletion", "recombination", "feature", "blastn_hit", "fasta_record"], emblfile=False)
	return track_number
	
	
	
	
def add_ordered_tab_to_diagram(filename, track_number, nameorder):
	
	features={"":[]}
	
	featurename=""
	names_to_add_feature_to=[]
	
	if filename.split(".")[-1]=="gz":
		print "need to add gzip functionality to tab reader. Until then, please unzip your input files."
		return
	
	record=tab_parser(open(filename,"r"))
	record.name=filename
	track_number=add_ordered_embl_to_diagram(record, track_number, filetype, nameorder, incfeatures=["i", "d", "li", "del", "snp", "misc_feature", "core", "cds", "insertion", "deletion", "recombination", "feature", "blastn_hit", "fasta_record"], emblfile=False)
	return track_number




			
################
# Main program #
################		

if __name__ == "__main__":


	#starttime=time.clock()

	#Get command line arguments

	(options, args) = main()


	if options.arrows==True:
		sigiltype="ARROW"
	else:
		sigiltype="BOX"
		
	if options.greytracks==True:
		greytracks=1
	else:
		greytracks=0
		
	
#	if options.alignment!="":
#		try:
#			alignment=read_alignment(options.alignment)
#		except StandardError:
#			DoError("Cannot open alignment file")
#			
#		
#		sequencenames={}
#		idtoindex={}
#		for x, record in enumerate(alignment):
#			sequencenames[x]=record.id
#			idtoindex[record.id]=x
	
	options.labelangle=int(options.labelangle)
	
	while options.labelangle>360:
		options.labelangle=options.labelangle-360
	while options.labelangle<0:
		options.labelangle=options.labelangle+360
		
	if (options.labelangle>=45 and options.labelangle<=135) or (options.labelangle>=225 and options.labelangle<=325):
		label_position="middle"
	elif options.labelangle<45 or options.labelangle>325:
		label_position="end"
	elif options.labelangle>135 and options.labelangle<225:
		label_position="start"
	
	#label_position="middle"
	
	#if options.labels==True:
	#	labels=1
	#else:
	#	labels=0
	
	if len(args)==0:
		print "No input files selected! Use -h or --help for help"
		sys.exit()
	
	end=int(options.end)
	beginning=int(options.beginning)	
	
	if end!=-1 and beginning>end:
		print "Beginning cannot be greater than end"
		sys.exit()
	
	
	#if a tree file is specified, open it and extract the names:
	
	if options.tree!="":
		if not os.path.isfile(options.tree):
			print "Cannot find file:", options.tree
			sys.exit()
		
		treestring=open(options.tree,"rU").read().strip()
		tree=Trees.Tree(treestring, rooted=True)
		taxonorder = tree.get_taxa(0)
	
	elif options.taxon_list!="":
		if not os.path.isfile(options.taxon_list):
			print "Cannot find file:", options.taxon_list
			sys.exit()
		
		taxonorder=[]
		for name in open(options.taxon_list,"rU").read().split("\n"):
			name=name.strip()
			if len(name)>0:
				taxonorder.append(name)
	
	bcffilter=options.bcffilter.split(',')
	if len(bcffilter)!=10:
		print "bcf filter list must contain 10 filters (minimum depth, mimimum strand depth, ratio of first to second base call, minimum base quality, minimum mapping quality, minimum allele frequency, strand bias p-value cutoff, base quality p-value cutoff bias, mapping bias p-value cutoff, tail distance bias p-value cutoff) followed by a letter code for filtering types of variant: a=include all i:just indels, s:just SNPs", options.taxon_list
		sys.exit()
	else:
		try: options.depth=int(bcffilter[0])
		except StandardError:
			DoError("Minimum mapping depth (bcffilter first value) must be an integer >= 0")
		if options.depth<0:
			print "Minimum mapping depth (bcffilter second value) must be >=0. Resetting to 0"
			options.depth=0
			
		try: options.stranddepth=int(bcffilter[1])
		except StandardError:
			DoError("Minimum strand mapping depth (bcffilter second value) must be an integer >= 0")
		if options.stranddepth<0:
			print "Minimum strand mapping depth (bcffilter second value) must be >=0. Resetting to 0"
			options.stranddepth=0
		if options.depth<(options.stranddepth*2):
			print "Minimum mapping depth must be at least double that for each strand. Resetting to", options.stranddepth*2
			options.stranddepth=options.stranddepth*2
			
		try: options.ratio=float(bcffilter[2])
		except StandardError:
			DoError("Ratio of first to second base call (bcffilter third value) must be between 0.5 and 1")
		if options.ratio<0.5 or options.ratio>1:
			DoError("Ratio of first to second base (bcffilter third value) must be greater than 0.5 and less than or equal to 1")
			
		try: options.QUAL=float(bcffilter[3])
		except StandardError:
			DoError("Minimum base quality (bcffilter fourth value) must be between 0 and 99")
		if options.QUAL<0 or options.QUAL>99:
			DoError("Base quality (bcffilter fourth value) must be between 0 and 99")
			
		try: options.MQUAL=float(bcffilter[4])
		except StandardError:
			DoError("Minimum mapping quality (bcffilter fifth value) must be between 0 and 99")
		if options.MQUAL<0 or options.MQUAL>99:
			DoError("Mapping quality (bcffilter fifth value) must be between 0 and 99")
			
		try: options.AF1=float(bcffilter[5])
		except StandardError:
			DoError("Minimum allele frequency for SNPs (bcffilter sixth value) must be between 0 and 1")
		if options.AF1<0 or options.AF1>1:
			DoError("Minimum allele frequency for SNPs (bcffilter sixth value) must be between 0 and 1")
		
		try: options.strand_bias=float(bcffilter[6])
		except StandardError:
			DoError("Strand bias p-value (bcffilter seventh value) must be between 0 and 1")
		if options.strand_bias<0 or options.strand_bias>1:
			DoError("p-value cutoff for strand bias (bcffilter seventh value) must be between 0 and 1")
		
		try: options.baseq_bias=float(bcffilter[7])
		except StandardError:
			DoError("Base quality bias p-value (bcffilter eigth value) must be between 0 and 1")
		if options.baseq_bias<0 or options.baseq_bias>1:
			DoError("p-value cutoff for base quality bias (bcffilter eigth value) must be between 0 and 1")
		
		try: options.mapping_bias=float(bcffilter[8])
		except StandardError:
			DoError("Mapping bias p-value (bcffilter ninth value) must be between 0 and 1")
		if options.mapping_bias<0 or options.mapping_bias>1:
			DoError("p-value cutoff for mapping bias (bcffilter ninth value) must be between 0 and 1")
		
		try: options.tail_bias=float(bcffilter[9])
		except StandardError:
			DoError("Tail bias p-value (bcffilter tenth value) must be between 0 and 1")
		if options.tail_bias<0 or options.tail_bias>1:
			DoError("p-value cutoff for tail distance bias (bcffilter tenth value) must be between 0 and 1")

	if not options.bcfvariants in ["a", "i", "s", "A", "I", "S"]:
		DoError("bcf filter code (bcffilter 11th value) must be a, i or s")
	
	
	#create translator object for translating artemis colours to GenomeDiagram colours
	
	translator = ColorTranslator()
	
	
	#get embl filename from command line and check the file exists
	
	
	gd_diagram = GenomeDiagram.Diagram()
	track_number=1
	
	filenames=[]
	suffix=options.suffix
	qualifier=options.qualifier
	
	if (options.tree=="" and options.taxon_list=="") or options.suffix=="":
		filenames=args
	else:
		for filename in args:
			
			suffixlen=0-len(suffix)
			
			if filename[suffixlen:]!=suffix or filename[:suffixlen] not in taxonorder:
				filenames.append(filename)
				print filename+" does not correspond to any taxon in the tree. Will add to the top of the diagram."
				
		for taxon in taxonorder:
			
			if taxon=="":
				continue
			elif taxon+suffix in args:
				filenames.append(taxon+suffix)
			else:
				filenames.append("")
				print "Warning! "+taxon+" is in tree, but you did not specify the corresponding "+taxon+suffix+" tab file. Will add blank track for this taxon."
		
		
	
	
	#print filenames
	for filename in filenames[::-1]:
	
		if filename=="":
			if track_number %2:
				gd_track_for_features = gd_diagram.new_track(track_number, scale=0, height=1, greytrack=greytracks, greytrack_labels=0)
			else:
				gd_track_for_features = gd_diagram.new_track(track_number, scale=0, height=1)
			
			#gd_track_for_features = gd_diagram.new_track(track_number, scale=0, height=1, greytrack=greytracks, greytrack_labels=0)
			gd_feature_set = gd_track_for_features.new_set()
			#gd_track_for_features = gd_diagram.new_track(track_number, scale=0, height=1, greytrack=greytracks, greytrack_labels=0)
			print "Adding blank track"
			#gd_track_for_features = gd_diagram.new_track(track_number, scale=0, height=1, greytrack=greytracks, greytrack_labels=0)
			#gd_feature_set = gd_track_for_features.new_set()
			track_number+=1
			continue
		
		if not os.path.isfile(filename):
			print "Cannot find file:", filename
			continue
		
		if filename.split(".")[-1]=="gz":
			filetype=filename.split(".")[-2].lower()
		else:
			filetype=filename.split(".")[-1].lower()
		
		
		if filetype=="plot":
			track_number=add_plot_to_diagram(filename,track_number, options.log, plotdata=[])
		elif filetype=="bam":
			track_number=add_bam_to_diagram(filename,track_number, options.log)
		elif filetype=="bcf":
			track_number=add_bcf_to_diagram(filename,track_number)
		elif filetype in ["embl", "gb"]:
			
#			if options.alignment!="" and filename.split(".")[0] in idtoindex.keys():
#				sequence=remove_gaps_from_sequence(alignment[idtoindex[filename.split(".")[0]]].seq)
#			
#				try:
#					emblrecord=open_annotation(filename, sequence)
#				except (StandardError, SimonError):
#					DoError("Cannot open annotation file "+filename+" please check the format")
#			else:
			try:
				emblrecord=open_annotation(filename)
			except (StandardError, SimonError):
				DoError("Cannot open annotation file "+filename+" please check the format")
				
			track_number=add_embl_to_diagram(emblrecord, track_number, filetype)
				
#			if (options.tree=="" and options.taxon_list=="") or options.qualifier=="":
#				track_number=add_embl_to_diagram(emblrecord, track_number, filetype)
#			else:
#				track_number=add_ordered_embl_to_diagram(emblrecord, track_number, filetype, taxonorder)
			
		elif filetype in ["tab", "art"]:
			if (options.tree=="" and options.taxon_list=="") or options.qualifier=="":
				track_number=add_tab_to_diagram(filename, track_number)
			else:
				track_number=add_ordered_tab_to_diagram(filename, track_number, taxonorder)
		
		elif filetype in ["fasta", "fna", "mfa", "dna", "faa", "fas", "phy", "phylip", "phylip2", "aln", "nxs", "nex", "nexus", "fastq", "ace"]:
			try:
				sequences=read_seq_file(filename)
			except StandardError:
				print "Cannot open file...skipping"
				continue
			
			colours=["10","11"]
			colourtouse=0
			start=1
			finish=0
	
			if track_number %2:
				gd_track_for_features = gd_diagram.new_track(track_number, scale=0, height=1, greytrack=greytracks, greytrack_labels=0)
			else:
				gd_track_for_features = gd_diagram.new_track(track_number, scale=0, height=1)
			gd_feature_set = gd_track_for_features.new_set()
			for x, record in enumerate(sequences):
				color = translator.artemis_color(colours[colourtouse])
				if colourtouse==0:
					colourtouse=1
				else:
					colourtouse=0
				finish=start+len(record)-1
				
				gd_feature_set.add_feature(SeqFeature(FeatureLocation(start, finish), strand=None), color=color, label=options.labels)
				start=finish+1
			print x, "features found in", filename
			track_number=track_number+1
			
	
	if track_number==1:
		DoError("No tracks found")
	
	if end==-1:
		end=gd_diagram.range()[1]
	pagewidth=29.7
	
	
	
#	
#	if options.thin:
#	
#		thin_max=int(float(end-beginning)*(float(pagewidth/200000)))
#		
#		
#		for track in gd_diagram.get_tracks():
#			
#			for sets in  track.get_sets():
#				features=sets.features
#				for feature_key in features.keys():
#					#print feature_key, features[feature_key].location.start, features[feature_key].location.end
#					if int(features[feature_key].location.nofuzzy_end)-int(features[feature_key].location.nofuzzy_start)<thin_max:
#						#print ""
#						sets.del_feature(feature_key)
#					else:
#						#features[feature_key].location.nofuzzy_end=int(features[feature_key].location.nofuzzy_end)-(thin_max/2)
#						#features[feature_key].location.nofuzzy_start=int(features[feature_key].location.nofuzzy_start)+(thin_max/2)
#						features[feature_key].location=FeatureLocation(start = features[feature_key].location._start._shift((thin_max/2)),end = features[feature_key].location._end._shift(((thin_max*-1)/2)))
#						#print features[feature_key].location
#					#print dir(features[feature_key].location)
#	#		
#	#
#	#print beginning, end
	
	if options.circular:
		#gd_diagram.move_track(1,3) # move track to make an empty space in the middle
		#print gd_diagram.track_offsets()
		gd_diagram.renumber_tracks(3)
		#print gd_diagram.track_offsets()
		gd_diagram.draw(format="circular", circular=True, pagesize=(20*cm,20*cm), fragments=options.fragments, start=beginning, end=end, track_size=options.tracksize)#, start=0, end=len(record))
		gd_diagram.write(options.outputfile+"_circular.pdf", "PDF")
	elif options.page=="estimate":
	
		estimated_width=(float(end-beginning)/10000)/options.fragments
		
		#print end, beginning, end-beginning, estimated_width
		
		if estimated_width<29.7:
			estimated_width=29.7
	
		x=2.0/estimated_width
		y=2.0/21
		if x>0.05:
			x=0.05
		if y>0.05:
			y=0.05
	
		gd_diagram.draw(format="linear", pagesize=(estimated_width*cm,21*cm), fragments=options.fragments, start=beginning, end=end, x=x, y=y, track_size=options.tracksize)#, start=0, end=len(record))
		gd_diagram.write(options.outputfile+"_linear.pdf", "PDF")
	else:
		gd_diagram.draw(format="linear", orientation=options.orientation, pagesize=options.page, fragments=options.fragments, start=beginning, end=end, track_size=options.tracksize)#, start=0, end=len(record))
		
		
		#gd_diagram.draw(format="linear", pagesize=(pagewidth*cm,10*cm), fragments=options.fragments, start=beginning, end=end, x=0, y=0)
		gd_diagram.write(options.outputfile+"_linear.pdf", "PDF")
	



