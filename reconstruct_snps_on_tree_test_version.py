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
	parser.add_option("-E", "--exclude", action="store", dest="exclude", help="Remove sequences less that float% missing data [default=%default]", default=50.0, type="float")
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
	parser.add_option("-c", "--recombination", action="store_true", dest="recombination", help="Run recombination analysis (slow)", default=False)
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
	elif options.exclude>=100 or options.exclude<=0:
		DoError('exclude below percentage (-E) must be between 0 and 100% (exclusive)')
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
			
			if base in ["N", "?"]:
				continue
			elif base=="-":
				gaplocations[record.id].append(x)
			elif base!='-' and base not in foundbases:# 
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
	toremove=[]

	
	
	for record in alignment:
		
		if len(str(record.seq).upper().replace("-", "").replace("N", "").replace("?",""))<(len(record.seq)*(options.exclude/100)):
			print "sequence record "+str(record.id)+" is > "+str(options.exclude)+"% unknown. Removing"
			toremove.append(record)
			continue
		name="seq"+str(count)
		
		sequencenames[name]=record.id
		convertnameback[record.id]=name
		seqnametoindex[name]=count-1
		count=count+1

	for record in toremove:
		alignment._records.remove(record)
	
		
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
				print treetaxa[x],
			print
		
		sys.exit()
			

	#print alignment.get_column(52545)	
	#print alignment.get_column(52546)
	
	#print alignment.get_column(52547)
	

	
	
	
	if options.outgroup!="" and options.outgroup!="None":
		print "Rooting tree on", options.outgroup
		sys.stdout.flush()
		tree.root_with_outgroup(outgroup=convertnameback[options.outgroup])
	elif options.outgroup!="None":
		print "Midpoint rooting tree"
		sys.stdout.flush()
		tree=midpoint_root(tree)
	

	colour_nodes_by_splitting(tree)
		
	#colour_nodes_by_tree_distances(tree, ladderize=None)
	#print a file with the node colours in RGB
	handle=open(prefix+"_node_colours.csv", "w")
	print_node_colours(tree, handle)
	handle.close()
	handle=open(prefix+"_recombination_key.tre","w")
	print >> handle, tree_to_figtree_string(tree, False, False, False, False, ladderize=None)
		
	handle.close()
	sys.exit()

	if options.outgroup!="None":
	
		chars = string.ascii_letters + string.digits
		tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))
	
		treestring= tree_to_string(tree, False, False, False, False)
		
		if treestring[-1]!=";":
			treestring=treestring=treestring+";"
		treestring="("+"(".join(treestring.split("(")[1:])
		handle = open(tmpname+".tre", "w")
		print >> handle, treestring
		handle.close()
		try:
			tree_string = open(tmpname+".tre").read()
		except IOError:
			DoError("Cannot open tree file "+tmpname+".tre")
		tree = Trees.Tree(tree_string, rooted=True)
		os.system("rm "+tmpname+".tre")


	

	#draw_ascii_tree(tree)
	

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
		
		
		
	
		try:
			emblrecord=open_annotation(options.embl, remove_gaps_from_sequence(alignment[refnum].seq))
		except StandardError:
			DoError("Cannot open annotation file "+options.embl+" please check the format")
		
		
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

		if len(alignment[refnum].seq)!=len(emblrecord):
			
			print "Printing extended annotation file in genbank format"
			sys.stdout.flush()
			
			emblrecord.seq=Seq(str(alignment[refnum].seq).upper(), IUPAC.IUPACAmbiguousDNA())
			emblCDSrecord.seq=Seq(str(alignment[refnum].seq).upper(), IUPAC.IUPACAmbiguousDNA())

			SeqIO.write([emblrecord], open(options.prefix+"_extended_annotation.gb","w"), "genbank")
			#sys.exit()
		#emblrecord=[]
		
	#draw_ascii_tree(tree)


	
	
	
	#if the gaps in the alignment are not real gaps, we need to change them to unknowns

	if not options.gaps:
		alignment=gap_to_unknown(alignment)
	
	
	print "Finding SNP positions in alignment"
	sys.stdout.flush()
	
	SNPlocations, consensus_sequence=Si_SNPs_temp.snp_locations_from_alignment(alignment)
	#SNPlocations, consensus_sequence=Si_SNPs_temp.snp_locations_from_alignment_tmp(alignment)
	

	#consensus_sequence="N"*alignment.get_alignment_length()
	

	
	#Add a consensus sequence to each node on the tree
	
	tree=add_object_to_all_nodes(tree, Seq(consensus_sequence), tree.root)#, Objecttype="annotation")

	
	
	print "Reconstructing sequences on tree"
	sys.stdout.flush()
	
	tree=parsimonious_sequence_reconstruction(tree, alignment, transformation=options.transformation, locations =SNPlocations)#locations=range(0,50000))#,locations=range(0,50000))# locations=range(89868,89870))#, sequence_Objecttype="annotation")

	
	#tree.display()
		
	if options.outgroup!="" and options.outgroup!="None":
		print "Pruning ", options.outgroup

		sys.stdout.flush()
		tree.prune(options.outgroup)
		recordtoremove=-1
		for x, record in enumerate(alignment):
			if record.id==options.outgroup:
				recordtoremove=x
			
		if recordtoremove>=0:
			newalignment=alignment[:recordtoremove]
			for record in alignment[recordtoremove+1:]:
				newalignment.add_sequence(record.id, str(record.seq))
			alignment=newalignment
		
	
	#tree.display()
	
	
	print "Finding homoplasies"
	sys.stdout.flush()
	
	tree=identify_homoplasies(tree)
	
	
	
	
	#Colour the nodes of the tree
	
	if options.recombination:
		#tree=colour_nodes_by_tree_position(tree, ladderize=None)
		tree=colour_nodes_by_splitting(tree)
		
		#colour_nodes_by_tree_distances(tree, ladderize=None)
		#print a file with the node colours in RGB
		handle=open(prefix+"_node_colours.csv", "w")
		print_node_colours(tree, handle)
		handle.close()
		handle=open(prefix+"_recombination_key.tre","w")
	
	
		print >> handle, tree_to_figtree_string(tree, False, False, False, False, ladderize=None)
		
		handle.close()

	if options.embl!="":
		
		
		print "Applying annotation to each branch of the tree"
		sys.stdout.flush()
		
		tree=apply_annotation_to_branches_new(tree, emblCDSrecord)
		
		#tree=apply_annotation_to_branches(tree, emblCDSrecord)
		
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
	handle=open(prefix+"_snps.txt", "w")
	print_summary(tree, handle, node=-1)
	handle.close()
	#sys.exit()
	handle=open(prefix+"_homoplasies_on_tree.tab","w")
	write_tab_output(tree, handle, node=-1, colour_snps_by="homoplasy")
	handle.close()
	handle=open(prefix+"_homoplasy_bases_on_tree.tab","w")
	write_tab_output(tree, handle, node=-1, colour_snps_by="homoplasy_bases")
	handle.close()
	handle=open(prefix+"_base_changes.tab","w")
	write_tab_output(tree, handle, node=-1, colour_snps_by="base")
	handle.close()
	
	
	#moving_window_recombination_detection(tree, prefix)
	#sys.exit()
	
	tree=branchlengths_to_SNP_count(tree)#, lengthtype="insertion_locations")
	
	tree=support_to_node_names(tree)
	treestring= tree_to_string(tree, False, False, False, False, ladderize=None)
	length=get_total_tree_length(tree)
	print "Total tree length =", length
	print "Printing tree with SNPs reconstructed"
	sys.stdout.flush()
	#print treestring
	handle = open(prefix+"_parsimony_steps"+".tre", "w")
	print >> handle, treestring
	handle.close()
	#sys.exit()
	draw_ascii_tree(tree)#, show_nodes=True)
	
	nodes=get_downstream_nodes(tree, tree.root)
	for node in nodes:
		if node==tree.root:
			continue
		if tree.is_terminal(node):
			nodename=tree.node(node).get_data().taxon
		else:
			nodename=str(node)
		nodefile=open(prefix+"_"+nodename+".csv","w")
		mydata=tree.node(node).get_data().comment
		print >> nodefile, "\t".join(['Gene Name', 'Location', 'Strand', 'Stop Codon Position Offset', 'SNPs', 'Synonymous', 'Non-synonymous', 'STOP->Non-STOP', 'Non-STOP->STOP', 'Others', 'Homoplasies', 'Insertions', 'Total Insertion Length', 'Deletions', 'Total Deletion Length', 'dNdS'])
		if 'annotation' in mydata:
			annotation=mydata['annotation']
			for x, gene in enumerate(annotation):
				name= gene['name']
				start=gene['location'][0]
				end=gene['location'][1]
				strand=gene['strand']
			
				astart=tree.node(tree.node(node).get_prev()).get_data().comment['annotation'][x]['location'][0]
				aend=tree.node(tree.node(node).get_prev()).get_data().comment['annotation'][x]['location'][1]
				diff=(aend-astart)-(end-start)
			
				SNPs=0
				hom=0
				Syn=0
				Non=0
				Stip=0
				Snop=0
				insertions=0
				ins_length=0
				deletions=0
				del_length=0
				others=0

				if 'SNP_locations' in mydata:
					for y in mydata['SNP_locations']:
						SNP=mydata['SNP_locations'][y]
						if y<=end and y>=start:
							SNPs+=1
							if SNP.homoplasy:
								hom+=1
						
							if SNP.codon_type=="S" or SNP.codon_type=="*":
								Syn+=1
							elif SNP.codon_type=="N":
								Non+=1
							elif SNP.codon_type=="M":
								Stip+=1
							elif SNP.codon_type=="O":
								Snop+=1
							else:
								#print "Invalid codon type", node, start, end, y, SNP.codon_type
								others+=1

				if 'insertion_locations' in mydata:
					for insertion in mydata['insertion_locations']:
						if (insertion[0]>=start and insertion[0]<=end):
							insertions+=1
							if insertion[1]<=end:
								ins_length+=(insertion[1]-insertion[0])+1
							else:
								ins_length+=(end-insertion[0])+1
						elif insertion[1]>=start and insertion[1]<=end:
							insertions+=1
							ins_length+=(insertion[1]-start)+1
			

				if 'deletion_locations' in mydata:
					for deletion in mydata['deletion_locations']:
						if (deletion[0]>=start and deletion[0]<=end):
							deletions+=1
							if deletion[1]<=end:
								del_length+=(deletion[1]-deletion[0])+1
							else:
								del_length+=(end-deletion[0])+1
						elif deletion[1]>=start and deletion[1]<=end:
							deletions+=1
							del_length+=(deletion[1]-start)+1


				if SNPs>0 and 'dNdS' in gene:
					try:
						dNdS=gene['dNdS']['dN']/gene['dNdS']['dS']
					except ZeroDivisionError:
						dNdS="Inf"
				else:
					dNdS="-"


			
				if SNPs>0 or insertions>0 or deletions>0 or diff!=0:
					print >> nodefile, "\t".join(map(str,[name, str(start+1)+".."+str(end+1), strand, diff, SNPs, Syn, Non, Stip, Snop, others, hom, insertions, ins_length, deletions, del_length, dNdS]))
		nodefile.close()




	nodeout=open(prefix+"_node_summary.csv","w")
	print >>nodeout, "\t".join(['Node', 'Taxon', 'SNPs', 'Synonymous', 'Non-Synonymous', 'STOP->Non-STOP', 'Non-STOP->STOP', 'Intergenic', 'Others', 'Homoplasies', 'Homoplasies with STI group 1', 'Homoplasies with STI group 2', 'Homoplasies with ocular group', 'Homoplasies with LGV group', 'Insertions', 'Total Insertion Length', 'Deletions', 'Total Deletion Length',  'dNdS'])
	for node in nodes:
		if node==tree.root:
			continue
		if tree.is_terminal(node):
			nodename=tree.node(node).get_data().taxon
		else:
			nodename="-"
		mydata=tree.node(node).get_data().comment

		SNPs=0
		hom=0
		Syn=0
		Non=0
		Stip=0
		Snop=0
		intergenic=0
		other=0
		insertions=0
		ins_length=0
		deletions=0
		del_length=0
		ocular_hom=0
		STI1_hom=0
		STI2_hom=0
		LGV_hom=0
		

		if 'dNdS' in mydata:
			if 'dNdS' in mydata['dNdS']:
				dNdS=mydata['dNdS']['dNdS']
			else:
				dNdS="?"
		else:
			dNdS="-"
		
		if 'SNP_locations' in mydata:
			for y in mydata['SNP_locations']:
				SNP=mydata['SNP_locations'][y]
				SNPs+=1
				if SNP.homoplasy:
					hom+=1
					foundLGV=False
					foundocular=False
					foundsti1=False
					foundsti2=False
					for h in SNP.homoplasies:
						if h[1]>=51 and not foundLGV:
							LGV_hom+=1
							foundLGV=True
						elif h[1]>=17 and h[1]<=22 and not foundocular:
							ocular_hom+=1
							foundocular=True
						elif h[1]>=4 and h[1]<=15 and not foundsti1:
							STI1_hom+=1
							foundsti1=True
						elif h[1]>=24 and h[1]<=49 and not foundsti2:
							STI2_hom+=1
							foundsti2=True
						
				if SNP.codon_type=="S" or SNP.codon_type=="*":
					Syn+=1
				elif SNP.codon_type=="N":
					Non+=1
				elif SNP.codon_type=="M":
					Stip+=1
				elif SNP.codon_type=="O":
					Snop+=1
				elif SNP.codon_type=="I":
					intergenic+=1
				else:
					other+=1
					#print SNP.codon_type


		if 'insertion_locations' in mydata:
			for insertion in mydata['insertion_locations']:
				insertions+=1
				ins_length+=(insertion[1]-insertion[0])+1

		if 'deletion_locations' in mydata:
			for deletion in mydata['deletion_locations']:
				deletions+=1
				del_length+=(deletion[1]-deletion[0])+1
		
				


		print >>nodeout, "\t".join(map(str,[node, nodename, SNPs, Syn, Non, Stip, Snop, intergenic, other, hom, STI1_hom, STI2_hom, ocular_hom, LGV_hom, insertions, ins_length, deletions, del_length,  dNdS]))
	nodeout.close()

	homoplasic_sites=get_homoplasic_sites(tree)
	homoplasic_sites.sort()

	if options.recombination:
		print "Finding recombination blocks"
		sys.stdout.flush()
		#moving_window_recombination_detection(tree, prefix)
		#find_clusters_on_branches(tree,prefix=prefix, homoplasy=True, locations =SNPlocations)
		tree=find_clusters_on_all_branches(tree,prefix=prefix, homoplasy=True, locations=SNPlocations)
		#sys.exit()
		#draw_ascii_tree(tree)
		#sys.exit()
	
	
	# reset homoplasic sites
	tree=reset_homoplasies(tree)
	
	# recalculate homoplasic sites
	tree=identify_homoplasies(tree, original=False)
	
	#print homoplasic_sites
	handle=open(prefix+"_homoplasies_on_tree_after_recombination_detection.tab","w")
	write_tab_output(tree, handle, node=-1, colour_snps_by="homoplasy")
	handle.close()
	#tree=branchlengths_to_SNP_count(tree, lengthtype="insertion_locations")
	
	tree=branchlengths_to_SNP_count(tree)#, lengthtype="insertion_locations")
	treestring= tree_to_string(tree, False, False, False, False, ladderize=None)
	length=get_total_tree_length(tree)
	print "Total tree length =", length
	print "Printing tree with SNPs reconstructed and recombinations removed"
	sys.stdout.flush()
	#print treestring
	handle=open(prefix+"_snps_post_recombination.txt", "w")
	print_summary(tree, handle, node=-1)
	handle.close()
	handle = open(prefix+"_parsimony_steps_minus_recombination"+".tre", "w")
	print >> handle, treestring
	handle.close()
	#sys.exit()
	draw_ascii_tree(tree)#, show_nodes=True)
	
	if options.embl!="" and options.pdf:
		
		genome_diagram_for_tree(tree, prefix+"_changed_genes_on_tree.pdf", ladderize=None, referenceObject=emblCDSrecord, printtype="changed_genes")
		genome_diagram_for_tree(tree, prefix+"_homoplasies_on_tree.pdf", ladderize=None, colourby="homoplasy", referenceObject=emblCDSrecord)
		genome_diagram_for_tree(tree, prefix+"_snps_on_tree.pdf", ladderize=None, referenceObject=emblCDSrecord)
		genome_diagram_for_tree(tree, prefix+"_base_changes.pdf", ladderize=None, referenceObject=emblCDSrecord, colourby="base", fragments=1, locations=homoplasic_sites)
		genome_diagram_for_tree(tree, prefix+"_taxa_patterns.pdf", ladderize=None, referenceObject=emblCDSrecord, colourby="taxa", fragments=1, locations=homoplasic_sites)
	elif options.pdf:
		genome_diagram_for_tree(tree, prefix+"_homoplasies_on_tree.pdf", ladderize=None, colourby="homoplasy")
		genome_diagram_for_tree(tree, prefix+"_snps_on_tree.pdf", ladderize=None)
		genome_diagram_for_tree(tree, prefix+"_base_changes.pdf", ladderize=None, colourby="base", fragments=1, locations=homoplasic_sites)
		genome_diagram_for_tree(tree, prefix+"_taxa_patterns.pdf", ladderize=None, colourby="taxa", fragments=1, locations=homoplasic_sites)
	
	
	
		#	handle=open(prefix+"_recombinations.tab","w")
		#	get_homoplasy_blocks(tree, handle)
		#	handle.close()
	
	
	

	#tree.display()
	
	print time.clock()-starttime
	sys.exit()
	
	
	
