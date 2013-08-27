#!/usr/bin/env python
import string, re, copy
import os, sys
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Nexus import Trees, Nodes
from Bio.Align import AlignInfo
from Bio.Align.Generic import Alignment
from Bio.Alphabet import IUPAC, Gapped
from Bio.SeqFeature import SeqFeature, FeatureLocation
from optparse import OptionParser
from random import *
from Bio.Alphabet import IUPAC
#sys.path.extend(map(os.path.abspath, ['/usr/lib/python2.4/site-packages/']))
#sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/lib/python2.5/site-packages/']))
#from scipy.stats import chi2
#from ghmm import *

sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_nexus import *
from Si_general import *
from Si_SeqIO import *
import Si_SNPs_temp

import time

#import pylab
#import numpy

#Requires Biopython and Pysam


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-a", "--alignment", action="store", dest="alignment", help="alignment file name", default="", metavar="FILE")
	parser.add_option("-t", "--tree", action="store", dest="tree", help="tree file", default="", metavar="FILE")
	parser.add_option("-o", "--outgroup", action="store", dest="outgroup", help="outgroup", default="")
	parser.add_option("-p", "--prefix", action="store", dest="prefix", help="prefix for output files", default="")
	parser.add_option("-T", "--transformation", action="store", dest="transformation", help="transformation type (acctran, deltran or ML). [Default= %default]", default="acctran", type="choice", choices=["acctran","deltran", "ML"])
	parser.add_option("-R", "--RAxML", action="store_true", dest="runtree", help="run phylogeny with RAxML [default=%default]", default=False)
	parser.add_option("-m", "--model", action="store", dest="model", help="Model of evolution to use. [Default= %default]", default="GTRGAMMA", type="choice", choices=["GTRGAMMA","GTRGAMMAI", "GTRCAT", "GTRMIX", "GTRMIXI"])
	parser.add_option("-b", "--bootstrap", action="store", dest="bootstrap", help="Number of bootstrap replicates (0 = do not run bootstrap). [Default= %default]", default=100, type="int", metavar="int")
	parser.add_option("-e", "--embl", action="store", dest="embl", help="Embl/genbank annotation file for reference strain (for dN/dS etc.)", default="", metavar="FILE")
	parser.add_option("-r", "--reference", action="store", dest="reference", help="Name of reference sequence relating to the embl file", default="", metavar="string")
	parser.add_option("-d", "--dNdS", action="store_true", dest="dnds", help="Calculate dN/dS for each branch (requires embl file) [default=%default]", default=False)
	parser.add_option("-g", "--gaps", action="store_true", dest="gaps", help="Gaps (-) are real", default=False)
	parser.add_option("-P", "--pdf", action="store_true", dest="pdf", help="Print pdfs (can be slow)", default=False)
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


#	if options.embl!="" and options.reference=='':
#		DoError('No reference selected! If you give an embl file you must specify which sequence it is linked to (i.e. the reference)')
	if options.alignment=='':
		DoError('No alignment file selected!')
	elif not os.path.isfile(options.alignment):
		DoError('Cannot find file '+options.alignment+'!')
	elif options.embl!='' and not os.path.isfile(options.embl):
		DoError('Cannot find file '+options.embl+'!')
	elif options.bootstrap>10000 or options.bootstrap<0:
		DoError('Number of bootstrap replicates (-b) must be between 0 (do not run bootstrap) and 10,000!')
	elif options.tree=="" and not options.runtree:
		options.runtree=True
	elif options.tree!="" and options.runtree:
		print "!!!Warning: Treefile provided and option to create tree with RAxML selected. Using user tree. RAxML will not be run!!!"
		options.runtree=False	
	if options.prefix=='':
		options.prefix=options.alignment.split("/")[-1].split(".")[0]



#	while os.path.isfile(options.outfile+".aln") and options.overwrite==False:
#		outopt=""
#		outopt=raw_input('\nOutput files with chosen prefix already exist.\n\nWould you like to overwrite (o), choose a new output file prefix (n) or quit (Q): ')
#		if outopt=='Q':
#			sys.exit()
#		elif outopt=="o":
#			break
#		elif outopt=="n":
#			options.outfile=raw_input('Enter a new output file prefix: ')
		
		
	return


	
	
#####################################################
# Function to identify snps and gaps from alignment #
#####################################################

def Find_SNP_and_gap_locations(alignment):
	
	gaplocations={}
	
	for record in alignment:
		gaplocations[record.id]=[]
	
	summary_align = AlignInfo.SummaryInfo(alignment)
	
	print "Identifying SNP and gap locations"
	sys.stdout.flush()
	
	count=0
	total=0.0
	
	for x in range(0,alignment.get_alignment_length()):

		count=count+1
		if count==10000:
			total=total+count
			count=0
			print "%.2f%% complete\r" % (100*(total/alignment.get_alignment_length())),
			sys.stdout.flush()

		
		foundbases=[]
		for record in alignment:
			base=record.seq[x].upper()
			
			if base=="-":
				gaplocations[record.id].append(x)
			if base!='-' and base not in foundbases:# 
				foundbases.append(base)
			if len(foundbases)>1:
				SNPlocations.append(x)
				break
	

	print "100.00% complete\n"#Found %d SNP locations" % len(SNPlocations),
	sys.stdout.flush()
	return SNPlocations, gaplocations




############################################
# Function to identify snps from alignment #
############################################

def Find_SNP_locations(alignment, startinglocations):
	
	SNPlocations=[]
	
	summary_align = AlignInfo.SummaryInfo(alignment)
	
	print "Identifying SNP locations"
	sys.stdout.flush()
	
	count=0
	total=0.0
	
	for x in startinglocations:

		count=count+1
		if count==10000:
			total=total+count
			count=0
			print "%.2f%% complete\r" % (100*(total/alignment.get_alignment_length())),
			sys.stdout.flush()

		
		foundbases=[]
		for record in alignment:
			base=record.seq[x].upper()
			
			if base!='-' and base not in foundbases:# 
				foundbases.append(base)
			if len(foundbases)>1:
				SNPlocations.append(x)
				break
	

	print "100.00% complete\n"#Found %d SNP locations" % len(SNPlocations),
	sys.stdout.flush()
	return SNPlocations





			
################
# Main program #
################		

if __name__ == "__main__":


	starttime=time.clock()

	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)
	
	prefix=options.prefix
	
	
	#Read the alignment file
	
	try:
		alignment=read_alignment(options.alignment)
	except StandardError:
		DoError("Cannot open alignment file")



	sequencenames={}
	convertnameback={}
	seqnametoindex={}
	count=1
	
	for record in alignment:

		name="seq"+str(count)
		
		sequencenames[name]=record.id
		convertnameback[record.id]=name
		seqnametoindex[name]=count-1
		count=count+1

		
#	tree=run_phyML(alignment, bootstrap=options.bootstrap, datatype="DNA", model="GTR", gamma=True, pinvar=False, cleanup=True)
#	tree.display()
#	sys.exit()

	if options.runtree:
		#If the user has chosen to run the tree, run RAxML to get the tree
		tree=run_RAxML(alignment, model=options.model, bootstrap=options.bootstrap)#, cleanup=False)
		
			
	else:
		#Else, read the tree file
		print "Reading tree file"
		sys.stdout.flush()
		
		try:
			tree_string = open(options.tree).read()
		except IOError:
			DoError("Cannot open tree file "+options.tree)
		tree = Trees.Tree(tree_string, rooted=True)
	


	#Check if the tree and alignment have the same set of taxa
	treetaxa=tree.get_taxa(tree.root)
	alignmenttaxa=[]
	for sequence in alignment:
		alignmenttaxa.append(sequence.name)
	
	treetaxa.sort()
	alignmenttaxa.sort()
	if treetaxa!=alignmenttaxa:
		print "Error! Your tree and alignment have different sets of taxa:"
		
		
		if len(alignmenttaxa)>=len(treetaxa):
			maxtaxlen=len(alignmenttaxa)
		else:
			maxtaxlen=len(treetaxa)
		
		print "Alignment            Tree"
		
		for x in range(maxtaxlen):
			if x<len(alignmenttaxa):
				print alignmenttaxa[x][:20]+" "*(20-len(alignmenttaxa[x])),
			else:
				print " "*20,
			if x<len(treetaxa):
				print treetaxa[x]
		
		sys.exit()
			

	#print alignment.get_column(52545)	
	#print alignment.get_column(52546)
	
	#print alignment.get_column(52547)
	
	#Read the annotation file

	if options.embl!='':
		
		refnum=-1
		if options.reference!="":
			refnum=-1
			for x, taxon in enumerate(alignment):
				if taxon.id==options.reference:
					refnum=x
					break
			
			if refnum==-1:
				DoError("Reference "+options.reference+" is not in your alignment")
		else:
			refnum=0
			options.reference=alignment[refnum].id
			print "No reference specified. Using first sequence in alignment as reference"
		
		
		
	
		try:
			emblrecord=open_annotation(options.embl, remove_gaps_from_sequence(alignment[refnum].seq))
		except StandardError:
			DoError("Cannot open annotation file "+options.embl+" please check the format")
		
		if options.gaps:
			alnreflen=len(remove_gaps_from_sequence(alignment[refnum].seq))
		else:	
			alnreflen=len(alignment[refnum].seq)
			
		
		if len(emblrecord.seq)!=len(remove_gaps_from_sequence(alignment[refnum].seq)):
			print len(emblrecord.seq), len(remove_gaps_from_sequence(alignment[refnum].seq))
			DoError("The length of your reference in the alignment and embl files are different")
	
	
		try:
			ref_to_alignment, alignment_to_ref=get_ref_to_alignment_translations(options.reference, alignment)
		
		except StandardError:
			DoError("Reference "+options.reference+" is not in your alignment")
	
		
		
		#SNP_types_using_reference(SNPlocations, "JJA", SeqRecordObject, alignmentObject)
	
		#add dNdS to embl features
		
	
		emblCDSrecord=SeqRecord(emblrecord.seq, id=emblrecord.id, name=emblrecord.name, description=emblrecord.description)
		
		feature_number=0
		
		
		
		for x, feature in enumerate(emblrecord.features):
			#print emblrecord.features[x]
			#print feature
			feature=change_reference_location_to_alignment_location(feature, ref_to_alignment)
			emblrecord.features[x]=feature

			if feature.type=="CDS":
				#feature.qualifiers["dNdS"]={}
				emblCDSrecord.features.append(feature)
				feature_number+=1
				#print emblCDSrecord.features[-1]
			
			#print emblrecord.features[x]
		#sys.exit()
		if len(alignment[refnum].seq)!=len(emblrecord):
			
			print "Printing extended annotation file in genbank format"
			
			emblrecord.seq=Seq(str(alignment[refnum].seq).upper(), IUPAC.IUPACAmbiguousDNA())
			emblCDSrecord.seq=Seq(str(alignment[refnum].seq).upper(), IUPAC.IUPACAmbiguousDNA())

			SeqIO.write([emblrecord], open(options.prefix+"_extended_annotation.gb","w"), "genbank")
		emblrecord=[]
		
	
	
	if options.outgroup!="" and options.outgroup!="None":
		print "Rooting tree on", options.outgroup
		sys.stdout.flush()
		tree.root_with_outgroup(outgroup=convertnameback[options.outgroup])
	elif options.outgroup!="None":
		print "Midpoint rooting tree"
		sys.stdout.flush()
		midpoint_root(tree)



	#If the tree has just been created with RAxML, print it to file

	if options.runtree:
		#print the tree to file
		tree.name="RAxML_tree"
		treestring= tree_to_string(tree, False, False, False, False)
		
		if treestring[-1]!=";":
			treestring=treestring=treestring+";"
		treestring="("+"(".join(treestring.split("(")[1:])
		handle = open(prefix+"_RAxML"+".tre", "w")
		print >> handle, treestring
		handle.close()
		

	
	#if the gaps in the alignment are not real gaps, we need to change them to unknowns

	if not options.gaps:
		alignment=gap_to_unknown(alignment)
		
	SNPlocations, consensus_sequence=Si_SNPs_temp.snp_locations_from_alignment(alignment)
	

	#blank_sequence=Seq("N"*alignment.get_alignment_length())
	
	#Add a consensus sequence to each node on the tree
	
	tree=add_object_to_all_nodes(tree, Seq(consensus_sequence), tree.root)#, Objecttype="annotation")
	
	
	print "Reconstructing sequences on tree"
	sys.stdout.flush()
	tree=parsimonious_sequence_reconstruction(tree, alignment, transformation=options.transformation, locations =SNPlocations)#,locations=range(0,50000))# locations=range(89868,89870))#, sequence_Objecttype="annotation")

#	tree=branchlengths_to_SNP_count(tree)#, lengthtype="insertion_locations")
#	tree=support_to_node_names(tree)
#	length=get_total_tree_length(tree)
#	print "Total tree length =", length

	#sys.exit()
	draw_ascii_tree(tree, show_nodes=True)
	
	print "Finding homoplasies"
	sys.stdout.flush()
	tree=identify_homoplasies(tree)
	
	
	#print tree.node(19).get_data().comment["sequence"]

	if options.embl!="":
		print "Applying annotation to each branch of the tree"
		sys.stdout.flush()
		
		
		tree=apply_annotation_to_branches(tree, emblCDSrecord)
		
		
		print "Annotating SNPS"
		sys.stdout.flush()
		
		tree=annotate_SNPs(tree)

		
		handle=open(prefix+"_changed_genes.tab","w")
		summary_handle=open(prefix+"_changed_genes.csv","w")
		print_changed_genes(tree, handle=handle, summary_handle=summary_handle)
		summary_handle.close()
		handle.close()

	
		if options.dnds:
#			print "Calculating dNdS for each branch"
#			sys.stdout.flush()
#			branch_dnds(tree)
			print "Calculating dNdS for each CDS on each branch"
			sys.stdout.flush()
			
			tree=dNdS_per_branch(tree)
		
			print "Calculating branch dNdS"
			sys.stdout.flush()
			handle=open(prefix+"_dNdS.txt","w")
			calculate_branch_dNdS(tree, handle)
			handle.close()
	

	print "Printing tab files"
	sys.stdout.flush()
	handle=open(prefix+"_snps_on_tree.tab","w")
	write_tab_output(tree, handle, node=-1, colour_snps_by="synonymous")
	handle.close()
	#sys.exit()
	handle=open(prefix+"_homoplasies_on_tree.tab","w")
	write_tab_output(tree, handle, node=-1, colour_snps_by="homoplasy")
	handle.close()
	
	
	
	tree=branchlengths_to_SNP_count(tree)#, lengthtype="insertion_locations")
	tree=support_to_node_names(tree)
	treestring= tree_to_string(tree, False, False, False, False)
	length=get_total_tree_length(tree)
	print "Total tree length =", length
	print "Printing tree with SNPs reconstructed"
	print treestring
	handle = open(prefix+"_parsimony_steps"+".tre", "w")
	print >> handle, treestring
	handle.close()
	
	#tree=branchlengths_to_SNP_count(tree, lengthtype="insertion_locations")
	
	
	if options.embl!="" and options.pdf:
		genome_diagram_for_tree(tree, prefix+"_homoplasies_on_tree.pdf", ladderize=None, colourby="homoplasy", referenceObject=emblCDSrecord)
		genome_diagram_for_tree(tree, prefix+"_snps_on_tree.pdf", ladderize=None, referenceObject=emblCDSrecord)
		genome_diagram_for_tree(tree, prefix+"_changed_genes_on_tree.pdf", ladderize=None, referenceObject=emblCDSrecord, printtype="changed_genes")
	elif options.pdf:
		genome_diagram_for_tree(tree, prefix+"_homoplasies_on_tree.pdf", ladderize=None, colourby="homoplasy")
		genome_diagram_for_tree(tree, prefix+"_snps_on_tree.pdf", ladderize=None)
	
	print time.clock()-starttime
	sys.exit()
	
	
	
