#!/usr/bin/env python
import string, re, numpy
import os, sys
from Bio.Nexus import Trees, Nodes
sys.path.extend(map(os.path.abspath, ['/nfs/pathogen/sh16_scripts/modules/']))
from Si_nexus import *

#node_order=["UW-3", "5", "SOTOND6_270710", "4", "J6276", "3", "G11222", "9", "G11074", "11", "G9301", "13", "G9768", "2", "JALI", "16", "CTA", "18", "CTB", "20", "CTR", "1", "F70", "23", "Bour_250310", "26", "E11023", "25", "SOTOND1_280710", "31", "Ds2923", "30", "SW4_210709", "34", "SW5_210709", "36", "SOTONF3_280710", "38", "SWFP-_140410", "29", "EC599_280710", "41", "E150", "43", "SW3_160410", "45", "SOTONE8_280710", "47", "SW2_Carma", "0", "RS2_160410", "56", "D4_140410", "58", "RS5_200410", "61", "L2b_SW_190410", "63", "Fr4_080410", "60", "RS1_160410", "55", "D10_140410", "54", "Fr1_2707101", "69", "Fr2_270710", "68", "RS6_190410", "73", "D7_140410", "72", "D6_140410", "76", "D8_140410", "53", "SA16_270710", "52", "RS4_160410", "80", "RS3_200410", "51", "L2", "83", "L2P-_080410", "50", "RS7_190410", "86", "RS8_160410"]
node_order=["L2_434_BU", "L2_P", "4", "L2b_UCH-2", "L2b_s300", "L2b_UCH-1", "L2b_s121", "L2b_s906", "18", "L2b_s11", "17", "15", "13", "11", "L2b_CV204", "L2_LST", "23", "L2b_795", "22", "10", "L2b_s750", "L2b_C2", "28", "L2b_C1", "27", "9", "L2b_8200_07", "8", "L1_SA16", "7", "3", "L1_224", "L1_115", "34", "2", "L3_404_LN", "L1_440_LN", "37", "1", "D_UW-3_CX", "D_SOTOND5", "47", "G_SOTONG1", "46", "K_SOTONK1", "D_SOTOND6", "51", "45", "G_11222", "44", "G_9301", "G_9768", "56", "G_11074", "55", "43", "Ia_SOTONIa3", "Ia_SOTONIa1", "61", "J_6276", "60", "42", "A_2497", "A_5291", "69", "A_363", "A_7249", "72", "68", "B_TZ1A828_OT", "67", "B_Jali20", "66", "A_HAR-13", "65", "41", "E_SW3", "E_150", "83", "E_SW2", "E_SOTONE8", "86", "82", "E_SOTONE4", "81", "E_Bour", "E_11023", "90", "80", "F_SOTONF3", "F_SW5", "95", "Ds_2923", "D_SOTOND1", "98", "94", "F_SW4", "93", "79", "F_70", "78", "40"]


donor_ca_db={}
donor_all_db={}

filename=sys.argv[1]

features=[]

if filename.split(".")[-1]=="gz":
	print "need to add gzip functionality to tab reader"
	sys.exit()
lines=open(filename,"r")
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
		node=words[1].split("=")[1]
		line=lines.next()
		words=line.split()
		if len(words)>=2 and words[1].split("=")[0][1:]=="possible_donors":
			donor_all_list="".join(words[1:]).split("=")[1].replace('"','').split(',')
		else:
			print "expected possible_donors after node", line
		line=lines.next()
		words=line.split()
		if len(words)==2 and words[1].split("=")[0][1:]=="donor_common_ancestor":
			donor_ca=words[1].split("=")[1]
		else:
			print "expected donor_common_ancestor after possible_donors"
		
		if not node in donor_ca_db:
			donor_ca_db[node]={}
		if not donor_ca in donor_ca_db[node]:
			donor_ca_db[node][donor_ca]=1
		else:
			donor_ca_db[node][donor_ca]+=1
		
		
		if not node in donor_all_db:
			donor_all_db[node]={}
		for donor in donor_all_list:
			if not donor in donor_all_db[node]:
				donor_all_db[node][donor]=1.0/len(donor_all_list)
			else:
				donor_all_db[node][donor]+=1.0/len(donor_all_list)
	
	
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
		if node in donor_ca_db:
			if nodeb in donor_ca_db[node]:
				linelist.append(str(donor_ca_db[node][nodeb]))
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
print >> Rfile, 'pdf(file="tmp_ca.pdf")'
#BLUE to RED:
#print >> Rfile, 'mycolors <- c("#FFFFFF","#0000FF", "#1100EE", "#2200DD", "#3300CC", "#4400BB", "#5500AA", "#660099", "#770088", "#880077", "#990066", "#AA0055", "#BB0044", "#CC0033", "#DD0022", "#EE0011", "#FF0000")'
#YELLOW to RED
print >> Rfile, 'mycolors <- c("#FFFFFF", "#FFFF00", "#FFEE00", "#FFDD00", "#FFCC00", "#FFBB00", "#FFAA00", "#FF9900", "#FF8800", "#FF7700", "#FF6600", "#FF5500", "#FF4400", "#FF3300", "#FF2200", "#FF1100", "#FF0000")'
print >> Rfile, 'd_heatmap <- heatmap.2(d_matrix, Rowv=NA, Colv=NA, dendrogram= c("none"), key=TRUE, keysize=1, trace="none", density.info=c("none"), margins=c(5, 5), col=mycolors)'
print >> Rfile, 'dev.off()'

Rfile.close()
os.system("R CMD BATCH tmp_r.R")

#maxvalue=0
#
#for node in node_order:
#	for nodeb in node_order:
#		if donor_all_db[node][nodeb]>maxvalue:
#			maxvalue
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
		if node in donor_all_db:
			if nodeb in donor_all_db[node]:
				linelist.append(str(donor_all_db[node][nodeb]))
			else:
				linelist.append("0")
		else:
				linelist.append("0")
	
	print >> Rinfile, "\t".join(linelist)

Rinfile.close()
#
Rfile=open("tmp_r.R", "w")
print >> Rfile, 'library(RColorBrewer)'
print >> Rfile, 'library(gplots)'
print >> Rfile, 'd <- read.table("tmp_r_in.txt", header=TRUE)'
print >> Rfile, 'row.names(d) <- d$Node'
print >> Rfile, 'd$Node <- NULL'
print >> Rfile, 'd_matrix <- data.matrix(d)'
print >> Rfile, 'pdf(file="tmp_all.pdf")'
#print >> Rfile, 'd_heatmap <- heatmap.2(d_matrix, Rowv=NA, Colv=NA, col = heat.colors(10), scale="column", margins=c(5,5))'
#, xlab = "X data", ylab = "Y data"
#BLUE to RED:
#print >> Rfile, 'mycolors <- c("#FFFFFF","#0000FF", "#1100EE", "#2200DD", "#3300CC", "#4400BB", "#5500AA", "#660099", "#770088", "#880077", "#990066", "#AA0055", "#BB0044", "#CC0033", "#DD0022", "#EE0011", "#FF0000")'
#YELLOW to RED
print >> Rfile, 'mycolors <- c("#FFFFFF", "#FFFF00", "#FFEE00", "#FFDD00", "#FFCC00", "#FFBB00", "#FFAA00", "#FF9900", "#FF8800", "#FF7700", "#FF6600", "#FF5500", "#FF4400", "#FF3300", "#FF2200", "#FF1100", "#FF0000")'
print >> Rfile, 'd_heatmap <- heatmap.2(d_matrix, Rowv=NA, Colv=NA, dendrogram= c("none"), key=TRUE, keysize=1, trace="none", density.info=c("none"), margins=c(5, 5), col=mycolors)'
print >> Rfile, 'dev.off()'

Rfile.close()
os.system("R CMD BATCH tmp_r.R")







try:
	tree_string = open(sys.argv[2]).read()
except IOError:
	DoError("Cannot open tree file "+options.tree)
tree = Trees.Tree(tree_string, rooted=True)

pairwise_distances={}

root=tree.root

allnodes=get_downstream_nodes(tree, root)
nametonum={}

for node in allnodes:
	if tree.is_terminal(node):
		nodename=tree.node(node).get_data().taxon
		nametonum[nodename]=node
	else:
		nodename=str(node)
		nametonum[nodename]=node
	
	pairwise_distances[nodename]={}
	
	comparisonnodes=get_nodes_not_downstream(tree, node)
	
	for nodeb in comparisonnodes:
		distance=tree.distance(node,nodeb)
		
		
		if tree.is_terminal(nodeb):
			nodebname=tree.node(nodeb).get_data().taxon
		else:
			nodebname=str(nodeb)

		pairwise_distances[nodename][nodebname]=distance


#distanceout=open("tmp_distance.csv","w")
#
#print >> distanceout, "to", "from", "distance", "ca_recombinations", "all_recombinations"
#
#for node in node_order:
#	if not node in pairwise_distances:
#		continue
#	for nodeb in node_order:
#		if not nodeb in pairwise_distances[node]:
#			continue
#		if node in donor_ca_db and nodeb in donor_ca_db[node]:
#			cavalue=donor_ca_db[node][nodeb]
#		else:
#			cavalue=0
#		if node in donor_all_db and nodeb in donor_all_db[node]:
#			allvalue=donor_all_db[node][nodeb]
#		else:
#			allvalue=0
#		print >> distanceout, node, nodeb, pairwise_distances[node][nodeb], cavalue, allvalue
#
#distanceout.close()


#distanceout=open("tmp_distance_for_hist.csv","w")
#
#print >> distanceout, "to", "from", "distance", "ca_recombinations", "all_recombinations"
#
#for node in node_order:
#	if not node in pairwise_distances:
#		continue
#	for nodeb in node_order:
#		if not nodeb in pairwise_distances[node]:
#			continue
#		if node in donor_ca_db and nodeb in donor_ca_db[node]:
#			cavalue=donor_ca_db[node][nodeb]
#		else:
#			cavalue=0
#		if node in donor_all_db and nodeb in donor_all_db[node]:
#			allvalue=donor_all_db[node][nodeb]
#		else:
#			allvalue=0
#		for x in range(0, cavalue):
#			print >> distanceout, pairwise_distances[node][nodeb]
#	#		print >> distanceout, node, nodeb, pairwise_distances[node][nodeb], cavalue, allvalue
#
#distanceout.close()



distanceout=open("tmp_pairwise_distances.csv","w")

print >> distanceout, "to", "from", "distance", "ca_recombinations", "all_recombinations"

for node in node_order:
	if not node in pairwise_distances:
		continue
	for nodeb in node_order:
		if not nodeb in pairwise_distances[node]:
			continue
		
		print >> distanceout, pairwise_distances[node][nodeb]
	#		print >> distanceout, node, nodeb, pairwise_distances[node][nodeb], cavalue, allvalue

distanceout.close()


distanceout=open("tmp_brlen_recombs.csv","w")
blout=open("tmp_brlen.csv","w")

print >> distanceout, "to", "from", "brlen", "ca_recombinations", "all_recombinations"

for node in node_order:
	if not node in pairwise_distances:
		continue
	tobrlen=tree.node(nametonum[node]).get_data().branchlength
	for nodeb in node_order:
		if not nodeb in pairwise_distances[node]:
			continue
		if node in donor_ca_db and nodeb in donor_ca_db[node]:
			cavalue=donor_ca_db[node][nodeb]
		else:
			cavalue=0
		if node in donor_all_db and nodeb in donor_all_db[node]:
			allvalue=donor_all_db[node][nodeb]
		else:
			allvalue=0
		
		frombrlen=tree.node(nametonum[nodeb]).get_data().branchlength
		
		for x in range(0, cavalue):
			print >> distanceout, tobrlen, frombrlen, pairwise_distances[node][nodeb]
	#		print >> distanceout, node, nodeb, pairwise_distances[node][nodeb], cavalue, allvalue
		print >> blout, tobrlen, frombrlen, pairwise_distances[node][nodeb]

distanceout.close()


