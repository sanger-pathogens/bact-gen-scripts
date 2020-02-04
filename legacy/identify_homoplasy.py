#!/usr/bin/env python
import string, re, gzip
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
#sys.path.extend(map(os.path.abspath, ['/usr/lib/python2.4/site-packages/']))
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/lib/python2.5/site-packages/']))
from scipy.stats import chi2
from ghmm import *

import pylab
import numpy

#Requires Biopython and Pysam
#################################
# Simple Error Printing Funtion #
#################################

def DoError(ErrorString):
	print "!!!Error:", ErrorString,"!!!"
	sys.exit()


###################################################################
# Function to create a SNP alignment #
###################################################################

def Create_SNP_alignment(alignment, SNPlocations):
	
	alphabet = Gapped(IUPAC.unambiguous_dna)

	SNPalignment = Alignment(alphabet)

	
	for record in alignment:
		SNPseq=""
		
		for base in SNPlocations:
			SNPseq=SNPseq+record.seq[base].replace("-","?")
		
		SNPalignment.add_sequence(record.id, SNPseq)
		
	
	return SNPalignment
	
	
#####################################################
# Function to identify snps and gaps from alignment #
#####################################################

def Find_SNP_and_gap_locations(alignment):
	
	SNPlocations=[]
	
	summary_align = AlignInfo.SummaryInfo(alignment)
	
	print "Identifying SNP locations"
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

			if base!='-' and base not in foundbases:# 
				foundbases.append(base)
			if len(foundbases)>1:
				SNPlocations.append(x+1)
				break
	

	print "100.00% complete\n"#Found %d SNP locations" % len(SNPlocations),
	sys.stdout.flush()
	return SNPlocations


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-a", "--alignment", action="store", dest="alignment", help="alignment file name", default="")
	parser.add_option("-t", "--tree", action="store", dest="tree", help="starting tree (optional)", default="")
	parser.add_option("-o", "--output", action="store", dest="outputfile", help="output file name", default="Homoplasies.csv")
	
	return parser.parse_args()


###########################################
# Function to create paml baseml.ctl file #
###########################################

def create_baseml_control_file(datafile, treefile):

	output=open("baseml.ctl","w")

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

    fix_alpha = 1   * 0: estimate alpha; 1: fix alpha at value below
        alpha =0* initial or fixed alpha, 0:infinity (constant rate)
        Malpha = 0   * 1: different alpha's for genes, 0: one alpha
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



################
# Main program #
################		

if __name__ == "__main__":


	#Get command line arguments

	(options, args) = main()


	#Do some checking of the input files
	
	if options.alignment=="":
		DoError("No alignment file specified")
#
	sequencenames={}
	convertnameback={}
	seqnametoindex={}	
	#Read the alignment file
		
	print "Reading alignment file"
	sys.stdout.flush()
	try:
		alignment = AlignIO.read(open(options.alignment), "fasta")
	except ValueError:
		DoError("Cannot open alignment file "+options.alignment+". Is it in the correct format?")
	
	prefix=options.alignment.split('.')[0]

	SNPlocations=Find_SNP_and_gap_locations(alignment)
		
	#create alignment of just the SNP sites
	
	SNPalignment=Create_SNP_alignment(alignment, SNPlocations)

	
		


	handle = open("SNPS_"+prefix+".phy", "w")
	print >> handle, len(SNPalignment), SNPalignment.get_alignment_length()
	
	count=1
	for record in SNPalignment:
	
		name="seq"+str(count)
#			
		sequencenames[name]=record.id
		convertnameback[record.id]=name
		seqnametoindex[name]=count-1
		
		print >> handle, name+"  "+record.seq
		count=count+1
	handle.close()
	
	
	
	if options.tree=="":
			
		if os.path.isfile("RAxML_info.SNPS_"+prefix):
			print "Removing old RAxML files"
			os.system("rm RAxML_*.SNPS_"+prefix)
			sys.stdout.flush()
			
		print "Running tree with RAxML"
		sys.stdout.flush()
		os.system("/software/pathogen/external/applications/RAxML/RAxML-7.0.4/raxmlHPC -f d -s SNPS_"+prefix+".phy -m GTRGAMMA -n SNPS_"+prefix+" > "+prefix+"temp.tmp")
		options.tree="RAxML_result."+"SNPS_"+prefix
	
	
	print "Reading tree file"
	sys.stdout.flush()
	
	try:
		tree_string = open(options.tree).read()
	except IOError:
		DoError("Cannot open tree file "+options.tree)
	tree = Trees.Tree(tree_string, rooted=True)
	
	
#	if options.outgroup!="":
#		print "Rooting tree on", options.outgroup
#		sys.stdout.flush()
#		tree.root_with_outgroup(outgroup=convertnameback[options.outgroup])

	
	treestring=tree.to_string(False, True, True, True)
	for name in sequencenames.keys():
		treestring=treestring.replace(sequencenames[name]+":", name+":")
	
	handle = open(prefix+".tre", "w")
	print >> handle, treestring+";"
	handle.close()
	
	print "Running PAML to reconstruct ancestral states"
	sys.stdout.flush()

	create_baseml_control_file("SNPS_"+prefix+".phy", prefix+".tre")
	
	#run paml
	
	os.system("/nfs/phd/mh10/software/paml41/bin/baseml > temp.tmp")
	
	output=open(options.outputfile,"w")
	print >> output, "Position,Number of changes,Changes seen"
	
	lines=open("rst","rU").read().split('Counts of changes at sites')[1].split("\n")
	
	for line in lines:
		words=line.split()
		
		numchanges=0
		changes=[]
		for x in words[1:-1]:
			if not "?" in x:
				numchanges=numchanges+1
				changes.append(x)
		
		if len(changes)>1:
			print >> output, str(SNPlocations[int(words[0])-1])+","+str(len(changes))+","+" ".join(changes)
		
#		if len(words)>2 and int(words[-1].replace("(","").replace(")",""))>1:
#			print >> output, str(SNPlocations[int(words[0])-1])+","+words[-1].replace("(","").replace(")","")+","+" ".join(words[1:-1])
	
	output.close()
	
	os.system("rm rates lnf mlb rst 2base.t temp.tmp rub rst1 baseml.ctl")
