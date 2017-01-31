#!/usr/bin/env python3

#################################
# Import some necessary modules #
#################################

import os, sys
import dendropy
from optparse import OptionParser, OptionGroup
from random import randrange
from scipy.stats import chisquare
import numpy as np


clade_colours=[["#9e0142", "#5e4fa2", "#ffffbf", "#f46d43", "#66c2a5", "#fee08b", "#abdda4", "#d53e4f", "#3288bd", "#fdae61", "#e6f598"],["#40004b", "#00441b", "#f7f7f7", "#9970ab", "#5aae61", "#c2a5cf", "#a6dba0", "#762a83", "#1b7837", "#e7d4e8", "#d9f0d3"],["#543005", "#003c30", "#f5f5f5", "#bf812d", "#35978f", "#f6e8c3", "#80cdc1", "#8c510a", "#01665e", "#dfc27d", "#c7eae5"]]


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
        
	group = OptionGroup(parser, "Required options")
	group.add_option("-t", "--tree", action="store", dest="tree", help="Input tree file in newick or nexus format", default="", metavar="FILE")
	group.add_option("-p", "--phenotype", action="store", dest="phenotype", help="Phenotype file (e.g. csv or tsv)", default="", metavar="FILE")
	group.add_option("-o", "--output_prefix", action="store", dest="prefix", help="Output prefix", default="treeBreaker")
	group.add_option("-c", "--column", action="store", dest="column", help="Number of the column in the phenotype file that contains the phenotype of interest [default= %default]", default=2, type="int")
	group.add_option("-s", "--separator", action="store", dest="separator", help="Separator for columns in input file [default= %default]", default=",")
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "treeBreaker options")
	group.add_option("-x", "--postburnin", action="store", dest="postburnin", help="Number of post-burnin iterations [default= %default]", default=500000, type="int")
	group.add_option("-y", "--burnin", action="store", dest="burnin", help="Number of burnin iterations [default= %default]", default=500000, type="int")
	group.add_option("-z", "--sample_iters", action="store", dest="sample_iters", help="Number of iterations between samples [default= %default]", default=1000, type="int")
	group.add_option("-S", "--seed", action="store", dest="seed", help="Starting seed for random number generator [default= %default]", default=-1, type="int")
	group.add_option("-m", "--missing", action="store_true", dest="include_missing", help="Include missing data as phenotype [default= %default]", default=False)
	
	parser.add_option_group(group)
	group = OptionGroup(parser, "General options")
	group.add_option("-P", "--posterior_cutoff", action="store", dest="posterior_cutoff", help="Posterior probability cutoff for defining clades based on breakpoints[default= %default]", default=0.5, type="float")
	group.add_option("-C", "--C", action="store", dest="colourset", help="Clade colour palette to use [default= %default]", default="1", type="choice", choices=["1","2","3"])
	group.add_option("-v", "--verbose", action="store_true", dest="verbose", help="Verbose [default= %default]", default=False)
	parser.add_option_group(group)
	
	return parser.parse_args()
	

def check_command_line_options(o):
	
	if  not os.path.isfile(o.tree):
		print("Cannot find file", o.tree)
		sys.exit()
		
	if  not os.path.isfile(o.phenotype):
		print("Cannot find file", o.phenotype)
		sys.exit()
	
	if o.column<2:
		print("Column number must be greater than 1")
		sys.exit()
	
	if o.postburnin<=0:
		print("treeBreaker post-brunin iterations should be greater than 0")
		sys.exit()
	if o.burnin<0:
		print("treeBreaker brunin iterations should be 0 or greater")
		sys.exit()
	if o.sample_iters<0:
		print("treeBreaker sample iterations should be 0 or greater")
		sys.exit()
	if o.posterior_cutoff<=0:
		print("Posterior probability cutoff should be greater than 0")
		sys.exit()
	if o.posterior_cutoff>1:
		print("Posterior probability cutoff should be less than or equal to 1")
		sys.exit()
	
	return


def read_tree(tfile):
	try:
		t=dendropy.Tree.get(path=tfile, schema="nexus", preserve_underscores=True)
	except dendropy.dataio.nexusreader.NexusReader.NotNexusFileError:
		t=dendropy.Tree.get(path=tfile, schema="newick", preserve_underscores=True)
	except:
		print("Cannot read tree file. It shoul dbe in newick or nexus format")
	return t

def rgbint2rgbtuple(RGBint):
	blue =  RGBint & 255
	green = (RGBint >> 8) & 255
	red =   (RGBint >> 16) & 255
	return (red, green, blue)	
	
################
# Main program #
################                

if __name__ == "__main__":

	(options, args) = main()
	check_command_line_options(options)
	
	treefile=options.tree
		
	phenotypefile=options.phenotype
	
	output_prefix=options.prefix
	
	column=options.column-1
	
	separator=options.separator
	
	postburnin=options.postburnin
	
	burnin=options.burnin
	
	posterior_cutoff=options.posterior_cutoff
	
	sample_iters=options.sample_iters
	
	colourset=int(options.colourset)-1
	
	if options.seed<0:
		seed=randrange(1,99999)
	else:
		seed=options.seed
	
	verbose=options.verbose
	
	treeBreakerOptions=["-x", postburnin, "-y", burnin, "-z", sample_iters, "-S", seed]
	if verbose:
		treeBreakerOptions.append("-v")
	treeBreakerOptionString=' '.join(map(str,treeBreakerOptions))
	
	tree=read_tree(treefile)
	
	tree_taxa=list(tree.taxon_namespace.labels())
	
	phenotypes=set([])
	taxon_phenotypes={}
	missing_phenotypes=[]
	not_in_tree=[]
	for line in open(phenotypefile, "rU", encoding='ISO-8859-1'):
		words=line.strip().split(separator)
		if len(words)<column:
			if options.verbose:
				print("Invalid line in phenotype file (Too few columns):", line.strip())
		if not words[0] in tree_taxa:
			if options.verbose:
				print(words[0], "not in tree. Skipping")
			not_in_tree.append(words[0])
			continue
		if words[column]=="":
			if options.verbose:
				print(words[0], "has missing phenotype")
			missing_phenotypes.append(words[0])
			if options.include_missing:
				taxon_phenotypes[words[0]]="missing"
				phenotypes.add("missing")
			continue
		taxon_phenotypes[words[0]]=words[column]
		phenotypes.add(words[column])
	
	phenotype={}
	for x, p in enumerate(phenotypes):
		phenotype[p]=x
	
	print(len(taxon_phenotypes), "taxa have phenotypes")
	if options.include_missing:
		pruned_tree=tree
	else:
		pruned_tree=tree.extract_tree_without_taxa_labels(labels=missing_phenotypes)
	pruned_tree.write(path=output_prefix+".nwk", schema="newick", unquoted_underscores=True)
	
	output=open(output_prefix+".tab", "w")
	for t in taxon_phenotypes:
		print(t+"\t"+str(phenotype[taxon_phenotypes[t]]), file=output)
	output.close()
	
	print("Phenotype key:")
	for p in phenotype:
		print(str(phenotype[p])+": "+p)
	
	print("Running treeBreaker")
	os.system("treeBreaker "+treeBreakerOptionString+" "+output_prefix+".nwk "+output_prefix+".tab "+output_prefix+".out")
	
	lines = open(output_prefix+".out","rU").readlines()
	treeBreaker_tree = lines[-1]
	
	treeBreaker=dendropy.Tree.get_from_string(treeBreaker_tree.replace("{", "[&").replace("}", "]").replace("|", ","), schema="newick", preserve_underscores=True)
	
	
	x=0
	clusternum=1
	leafset=set([])
	for l in treeBreaker.leaf_node_iter():
		leafset.add(l.taxon.label)
	for node in treeBreaker.postorder_node_iter():
		if node==treeBreaker.seed_node:
			continue
		for a in node.annotations:
			if a.name=="posterior":
				leaves=[]
				for l in node.leaf_iter():
					leaves.append(l.taxon.label)
				if float(a.value)>=posterior_cutoff:
					node.annotations.add_new(name="!hilight", value="{1,1,"+clade_colours[colourset][x]+"}")
					x+=1
					if x>=len(clade_colours[colourset]):
						x=0
					for l in node.leaf_iter():
						if hasattr(l, 'treeBreaker_cluster'):
							if not hasattr(l, 'treeBreaker_parent_cluster'):
								l.treeBreaker_parent_cluster=clusternum
						if not hasattr(l, 'treeBreaker_cluster'):
							l.treeBreaker_cluster=clusternum
							l.treeBreaker_cluster_posterior=float(a.value)
						
					clusternum+=1
						
				R=round(255*float(a.value))
				G=0
				B=round(255*(1.0-float(a.value)))
				
				Rint=R*65536
				Gint=0
				Bint=B
				RGBint=Rint+Gint+Bint
				
				node.annotations.add_new(name="!color", value="#"+str(RGBint))
				if G!=0:
					print(R, G, B, rgbint2rgbtuple(RGBint))
		
	
	clusters={}			
	for l in treeBreaker.leaf_node_iter():
		if not hasattr(l, 'treeBreaker_cluster'):
			cluster=0
			parent_cluster=-1
			posterior=1.0
		else:
			cluster=l.treeBreaker_cluster
			if not hasattr(l, 'treeBreaker_parent_cluster'):	
				parent_cluster=0
			else:
				parent_cluster=l.treeBreaker_parent_cluster
			posterior=l.treeBreaker_cluster_posterior
		if not cluster in clusters:
			clusters[cluster]={"posterior":posterior, "parent":parent_cluster, "members":[], "phenotype_counts":{}}
			for p in phenotypes:
				clusters[cluster]["phenotype_counts"][p]=0
		label=l.taxon.label
		clusters[cluster]["members"].append([label, taxon_phenotypes[label]])
		clusters[cluster]["phenotype_counts"][taxon_phenotypes[label]]+=1
	cluster_output=open(output_prefix+".clusters", "w")
	for cluster in clusters:
		print("cluster_"+str(cluster)+"="+','.join([x[0] for x in clusters[cluster]["members"]]), sep=",", file=cluster_output)
		if clusters[cluster]["parent"] in clusters:
			obs=[]
			exp=[]
			for p in phenotypes:
				obs.append(clusters[cluster]["phenotype_counts"][p])
				exp.append(clusters[clusters[cluster]["parent"]]["phenotype_counts"][p])
			cs=chisquare(np.array(obs), f_exp=np.array(exp))
			print(str(clusters[cluster]["parent"])+"->"+str(cluster), clusters[clusters[cluster]["parent"]]["phenotype_counts"], clusters[cluster]["phenotype_counts"], clusters[cluster]["posterior"], cs)

	
	treeBreaker.write(path=output_prefix+".nexus", schema="nexus", unquoted_underscores=True)
	
	#os.system("~sh16/scripts/iCANDY.py -t "+output_prefix+".nexus -s posterior -J posterior -O portrait -p A0 -z circle -m "+phenotypefile+" -C "+str(column+1)+" -a 2 -o "+output_prefix+".pdf")
	os.system("~sh16/scripts/iCANDY.py -t "+output_prefix+".nexus -s posterior -j -O portrait -p A0 -m "+phenotypefile+" -C ,"+str(column+1)+" -a 2 -o "+output_prefix+".pdf")
	
	
	
	sys.exit()
