#!/usr/bin/env python
import string, re, numpy
import os, sys
from Bio.Nexus import Trees, Nodes
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_nexus import *
from Bio.Graphics.GenomeDiagram._Colors import ColorTranslator

node_order=["UW-3", "5", "SOTOND6_270710", "4", "J6276", "3", "G11222", "9", "G11074", "11", "G9301", "13", "G9768", "2", "JALI", "16", "CTA", "18", "CTB", "20", "CTR", "1", "F70", "23", "Bour_250310", "26", "E11023", "25", "SOTOND1_280710", "31", "Ds2923", "30", "SW4_210709", "34", "SW5_210709", "36", "SOTONF3_280710", "38", "SWFP-_140410", "29", "EC599_280710", "41", "E150", "43", "SW3_160410", "45", "SOTONE8_280710", "47", "SW2_Carma", "0", "RS2_160410", "56", "D4_140410", "58", "RS5_200410", "61", "L2b_SW_190410", "63", "Fr4_080410", "60", "RS1_160410", "55", "D10_140410", "54", "Fr1_2707101", "69", "Fr2_270710", "68", "RS6_190410", "73", "D7_140410", "72", "D6_140410", "76", "D8_140410", "53", "SA16_270710", "52", "RS4_160410", "80", "RS3_200410", "51", "L2", "83", "L2P-_080410", "50", "RS7_190410", "86", "RS8_160410"]

#Just LGVs
#node_order=["RS2_160410", "56", "D4_140410", "58", "RS5_200410", "61", "L2b_SW_190410", "63", "Fr4_080410", "60", "RS1_160410", "55", "D10_140410", "54", "Fr1_2707101", "69", "Fr2_270710", "68", "RS6_190410", "73", "D7_140410", "72", "D6_140410", "76", "D8_140410", "53", "SA16_270710", "52", "RS4_160410", "80", "RS3_200410", "51", "L2", "83", "L2P-_080410", "50", "RS7_190410", "86", "RS8_160410"]



donor_db={}


filename=sys.argv[1]

features=[]

if filename.split(".")[-1]=="gz":
	print "need to add gzip functionality to tab reader"
	sys.exit()
lines=open(filename,"r")

maxvalue=0

for line in lines:
	words=line.split()
	
#	if len(words)==3 and words[1].lower() =="recombination":
#		try: features.append((feature,colour))
#		except NameError: pass #print "here"
		
		
		
#perhaps add something back in here to get lengths/number of snps etc in recombinations			
#		position=words[2].replace("complement(","").replace(")","")
#		if len(position.split(".."))>1:
#			start=int(position.split("..")[0])
#			end=int(position.split("..")[1])
#			if start<options.beginning and end>options.beginning:
#				start=options.beginning
#			if end>options.end and start<options.end:
#				end=options.end
#			feature = SeqFeature(FeatureLocation(start, end), strand=None)
#		else:
#			feature = SeqFeature(FeatureLocation(int(position), int(position)), strand=None)
	
	
	if len(words)==2 and words[1].split("=")[0][1:]=="node":
		node=words[1].split("=")[1].replace('"','').split(">")[1]
		for x in range(0,2):
			line=lines.next()
		words=line.split()
		if len(words)>=2 and words[1].split("=")[0][1:]=="homoplasy":
			donor_all_list=" ".join(words[1:]).split("=")[1].replace('"','').split(',')
			for x in range(0,len(donor_all_list)):
				donor_all_list[x]=donor_all_list[x].split()[-1]
		
		
			if not node in donor_db:
				donor_db[node]={}
			for donor in donor_all_list:
				if not donor in donor_db[node]:
					donor_db[node][donor]=1
				else:
					donor_db[node][donor]+=1
				if donor_db[node][donor]>maxvalue and node in node_order and donor in node_order:
					maxvalue=donor_db[node][donor]
					maxpair=[node,donor]
					
print maxvalue, maxpair

	
Rinfile=open("tmp_r_in.txt", "w")

nodenames={}
nodenamelist=[]

for node in node_order:
	try:
		int(node)
		nodename="Node_"+node
	except ValueError:
		nodename=node
	nodenames[node]=nodename
	nodenamelist.append(nodename)


print >> Rinfile, "Node\t"+"\t".join(nodenamelist)
	
for node in node_order:
	linelist=[nodenames[node]]
	for nodeb in node_order:
		if node in donor_db:
			if nodeb in donor_db[node]:
				linelist.append(str(donor_db[node][nodeb]))
			else:
				linelist.append("0")
		else:
				linelist.append("0")
	
	print >> Rinfile, "\t".join(linelist)

Rinfile.close()

Rfile=open("tmp_r.R", "w")
print >> Rfile, 'library(RColorBrewer)'
print >> Rfile, 'library(gplots)'
print >> Rfile, 'd <- read.table("tmp_r_in.txt", header=TRUE)'
print >> Rfile, 'row.names(d) <- d$Node'
print >> Rfile, 'd$Node <- NULL'
print >> Rfile, 'd_matrix <- data.matrix(d)'
print >> Rfile, 'pdf(file="homoplasy_heatmap_all.pdf")'
#BLUE to RED:
#print >> Rfile, 'mycolors <- c("#FFFFFF","#0000FF", "#1100EE", "#2200DD", "#3300CC", "#4400BB", "#5500AA", "#660099", "#770088", "#880077", "#990066", "#AA0055", "#BB0044", "#CC0033", "#DD0022", "#EE0011", "#FF0000")'
#YELLOW to RED
#print >> Rfile, 'mycolors <- c("#FFFFFF", "#FFFF00", "#FFEE00", "#FFDD00", "#FFCC00", "#FFBB00", "#FFAA00", "#FF9900", "#FF8800", "#FF7700", "#FF6600", "#FF5500", "#FF4400", "#FF3300", "#FF2200", "#FF1100", "#FF0000")'


colors=[]


#add white zero
colors.append("#FFFFFF")

x=0

#BLUE magenta RED:
#bit=510.0/maxvalue
#while x<510:
#	green=0
#	if x<=255:
#		blue=255
#		red=x
#	else:
#		red=255
#		blue=510-x
#	
#	colors.append(RGBToHTMLColor((red, green, blue)))
#	x+=bit
	
	
#BLUE cyan green yellow RED:
#bit=1020.0/maxvalue
#while x<1020:
#	
#	if x<=255:
#		red=0
#		green=x
#		blue=255
#	elif x<=510:
#		red=0
#		green=255
#		blue=510-x	
#	elif x<=765:
#		red=x-510
#		green=255
#		blue=0
#	else:
#		red=255
#		green=1020-x
#		blue=0
#		
#	
#	
#	colors.append(RGBToHTMLColor((red, green, blue)))
#	x+=bit



#YELLOW to RED:
bit=255.0/maxvalue
while x<255:
	
	blue=0
	red=255
	green=255-x
	
	colors.append(RGBToHTMLColor((red, green, blue)))
	x+=bit





print >> Rfile, 'mycolors <- c("'+'",\n"'.join(colors)+'")'

print >> Rfile, 'd_heatmap <- heatmap.2(d_matrix, Rowv=NA, Colv=NA, dendrogram= c("none"), key=TRUE, keysize=1, trace="none", density.info=c("none"), margins=c(5, 5), col=mycolors)'
print >> Rfile, 'dev.off()'

Rfile.close()
os.system("R CMD BATCH tmp_r.R")
#
##maxvalue=0
##
##for node in node_order:
##	for nodeb in node_order:
##		if donor_all_db[node][nodeb]>maxvalue:
##			maxvalue
#Rinfile=open("tmp_r_in.txt", "w")
#
#nodenames={}
#nodenamelist=[]
#
#for node in node_order:
#	try:
#		int(node)
#		nodename="Node_"+node
#	except ValueError:
#		nodename=node
#	nodenames[node]=nodename
#	nodenamelist.append(nodename)
#
#
#print >> Rinfile, "Node\t"+"\t".join(nodenamelist)
#for node in node_order:
#	linelist=[nodenames[node]]
#	for nodeb in node_order:
#		if node in donor_all_db:
#			if nodeb in donor_all_db[node]:
#				linelist.append(str(donor_all_db[node][nodeb]))
#			else:
#				linelist.append("0")
#		else:
#				linelist.append("0")
#	
#	print >> Rinfile, "\t".join(linelist)
#
#Rinfile.close()
#
#Rfile=open("tmp_r.R", "w")
#print >> Rfile, 'library(RColorBrewer)'
#print >> Rfile, 'library(gplots)'
#print >> Rfile, 'd <- read.table("tmp_r_in.txt", header=TRUE)'
#print >> Rfile, 'row.names(d) <- d$Node'
#print >> Rfile, 'd$Node <- NULL'
#print >> Rfile, 'd_matrix <- data.matrix(d)'
#print >> Rfile, 'pdf(file="tmp_all.pdf")'
##print >> Rfile, 'd_heatmap <- heatmap.2(d_matrix, Rowv=NA, Colv=NA, col = heat.colors(10), scale="column", margins=c(5,5))'
##, xlab = "X data", ylab = "Y data"
##BLUE to RED:
##print >> Rfile, 'mycolors <- c("#FFFFFF","#0000FF", "#1100EE", "#2200DD", "#3300CC", "#4400BB", "#5500AA", "#660099", "#770088", "#880077", "#990066", "#AA0055", "#BB0044", "#CC0033", "#DD0022", "#EE0011", "#FF0000")'
##YELLOW to RED
#print >> Rfile, 'mycolors <- c("#FFFFFF", "#FFFF00", "#FFEE00", "#FFDD00", "#FFCC00", "#FFBB00", "#FFAA00", "#FF9900", "#FF8800", "#FF7700", "#FF6600", "#FF5500", "#FF4400", "#FF3300", "#FF2200", "#FF1100", "#FF0000")'
#print >> Rfile, 'd_heatmap <- heatmap.2(d_matrix, Rowv=NA, Colv=NA, dendrogram= c("none"), key=TRUE, keysize=1, trace="none", density.info=c("none"), margins=c(5, 5), col=mycolors)'
#print >> Rfile, 'dev.off()'
#
#Rfile.close()
#os.system("R CMD BATCH tmp_r.R")





Rinfile=open("tmp_r_in.txt", "w")

try:
	tree_string = open(sys.argv[2]).read()
except IOError:
	DoError("Cannot open tree file "+options.tree)
tree = Trees.Tree(tree_string, rooted=True)

pairwise_distances={}

root=tree.root

allnodes=get_downstream_nodes(tree, root)
nametonum={}
maxvalue=0
for node in allnodes:
	if tree.is_terminal(node):
		nodename=tree.node(node).get_data().taxon
		nametonum[nodename]=node
	else:
		nodename=str(node)
		nametonum[nodename]=node
	
	pairwise_distances[nodename]={}
	
	for nodeb in allnodes:
		distance=tree.distance(node,nodeb)
		
		
		if tree.is_terminal(nodeb):
			nodebname=tree.node(nodeb).get_data().taxon
		else:
			nodebname=str(nodeb)

		pairwise_distances[nodename][nodebname]=distance
		if distance>maxvalue and nodename in node_order and nodebname in node_order:
			maxvalue=distance
print maxvalue

print >> Rinfile, "Node\t"+"\t".join(nodenamelist)
	
for node in node_order:
	linelist=[nodenames[node]]
	for nodeb in node_order:
		if node in pairwise_distances:
			if nodeb in pairwise_distances[node]:
				linelist.append(str(pairwise_distances[node][nodeb]))
			else:
				linelist.append("0")
		else:
				linelist.append("0")
	
	print >> Rinfile, "\t".join(linelist)

Rinfile.close()

Rfile=open("tmp_r.R", "w")
print >> Rfile, 'library(RColorBrewer)'
print >> Rfile, 'library(gplots)'
print >> Rfile, 'd <- read.table("tmp_r_in.txt", header=TRUE)'
print >> Rfile, 'row.names(d) <- d$Node'
print >> Rfile, 'd$Node <- NULL'
print >> Rfile, 'd_matrix <- data.matrix(d)'
print >> Rfile, 'pdf(file="homoplasy_heatmap_all_distances.pdf")'
#BLUE to RED:
#print >> Rfile, 'mycolors <- c("#FFFFFF","#0000FF", "#1100EE", "#2200DD", "#3300CC", "#4400BB", "#5500AA", "#660099", "#770088", "#880077", "#990066", "#AA0055", "#BB0044", "#CC0033", "#DD0022", "#EE0011", "#FF0000")'
#YELLOW to RED
#print >> Rfile, 'mycolors <- c("#FFFFFF", "#FFFF00", "#FFEE00", "#FFDD00", "#FFCC00", "#FFBB00", "#FFAA00", "#FF9900", "#FF8800", "#FF7700", "#FF6600", "#FF5500", "#FF4400", "#FF3300", "#FF2200", "#FF1100", "#FF0000")'


colors=[]


#add white zero
colors.append("#FFFFFF")

x=0

#BLUE magenta RED:
bit=510.0/maxvalue
while x<510:
	green=0
	if x<=255:
		blue=255
		red=x
	else:
		red=255
		blue=510-x
	
	colors.append(RGBToHTMLColor((red, green, blue)))
	x+=bit
	
	
#BLUE cyan green yellow RED:
#bit=1020.0/maxvalue
#while x<1020:
#	
#	if x<=255:
#		red=0
#		green=x
#		blue=255
#	elif x<=510:
#		red=0
#		green=255
#		blue=510-x	
#	elif x<=765:
#		red=x-510
#		green=255
#		blue=0
#	else:
#		red=255
#		green=1020-x
#		blue=0
#		
#	
#	
#	colors.append(RGBToHTMLColor((red, green, blue)))
#	x+=bit



#YELLOW to RED:
#bit=255.0/maxvalue
#while x<255:
#	
#	blue=0
#	red=255
#	green=255-x
#	
#	colors.append(RGBToHTMLColor((red, green, blue)))
#	x+=bit





print >> Rfile, 'mycolors <- c("'+'",\n"'.join(colors)+'")'

print >> Rfile, 'd_heatmap <- heatmap.2(d_matrix, Rowv=NA, Colv=NA, dendrogram= c("none"), key=TRUE, keysize=1, trace="none", density.info=c("none"), margins=c(5, 5), col=mycolors)'
print >> Rfile, 'dev.off()'

Rfile.close()
os.system("R CMD BATCH tmp_r.R")

