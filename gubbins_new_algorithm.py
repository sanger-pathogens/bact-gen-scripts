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
import math
#sys.path.extend(map(os.path.abspath, ['/usr/lib/python2.4/site-packages/']))
#sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/lib/python2.5/site-packages/']))
#from scipy.stats import chi2
#from scipy import factorial
#from ghmm import *
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_SeqIO import *
from Si_nexus import draw_ascii_tree, tree_to_string, midpoint_root
import random
import numpy
import time

#Requires Biopython and Pysam


gaps_and_missing_set=set(["N", "?", "-"])


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-a", "--alignment", action="store", dest="alignment", help="alignment file name", default="")
	parser.add_option("-i", "--iterations", action="store", dest="maxiterations", type="int", help="maximum number of iterations to run [default=%default]", default=5)
	parser.add_option("-m", "--minsnps", action="store", dest="minsnps", type="int", help="minimum number of snps required on branch to run recombination detection algorithm [default=%default]", default=3)
	parser.add_option("-o", "--outgroup", action="store", dest="outgroup", help="outgroup", default="")
	parser.add_option("-p", "--prefix", action="store", dest="prefix", help="prefix for output files", default="")
	parser.add_option("-R", "--RAxML", action="store_true", dest="runtree", help="run phylogeny with RAxML [default=%default]", default=False)
	parser.add_option("-w", "--wholedata", action="store_true", dest = "wholedata", help="Rerun final phylogeny on whole data (rather than just SNP sites - better but slooooower", default=False)
	parser.add_option("-b", "--bootstrap", action="store_true", dest = "bootstrap", help="bootstrap final phylogeny", default=False)
	parser.add_option("-t", "--tree", action="store", dest="tree", help="starting tree (optional)", default="")
	parser.add_option("-A", "--alpha", action="store", dest="alpha", type="float", help="alpha parameter for gamma distribution (starting tree is specified) [default=%default]", default=1)
	parser.add_option("-r", "--randomisations", action="store", dest="permutations", help="number of randomisations to use to calculate p-value [default=%default]", default=100, type="int")
	parser.add_option("-c", "--cutoff", action="store", dest="cutoff", help="P-value cutoff [default=%default]", default=0.05, type="float")
	#parser.add_option("-x", "--oldPAML", action="store_true", dest="usepreviouspamlrun", help="start from previous paml data", default=False)
	parser.add_option("-g", "--geodesic", action="store_true", dest="geodesic", help="calculate geodesic distance between trees", default=False)
	#parser.add_option("-r", "--reference", action="store", dest="reference", help="reference embl file to add to final diagram (optional)", default="")
	
	return parser.parse_args()


#################################
# Simple Error Printing Funtion #
#################################

def DoError(ErrorString):
	print "!!!Error:", ErrorString,"!!!"
	sys.exit()


######################################
# Function to create a SNP alignment #
######################################

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
		if count==alignment.get_alignment_length()/100:
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




#####################################################
# Function to identify snps and gaps from alignment #
#####################################################

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
			if base in ["N", "?"]:
				continue
			elif base!='-' and base not in foundbases:# 
				foundbases.append(base)
			if len(foundbases)>1:
				SNPlocations.append(x)
				break
	

	print "100.00% complete\n"#Found %d SNP locations" % len(SNPlocations),
	sys.stdout.flush()
	return SNPlocations




###########################################
# Function to create paml baseml.ctl file #
###########################################

def create_baseml_control_file(datafile, treefile, alpha):

	output=open("baseml.ctl","w")
	
	if alpha>10:
		alpha=0

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

    fix_alpha = 1   * 0: estimate alpha; 1: fix alpha at value below"""
	print >> output, "        alpha =", alpha, "* initial or fixed alpha, 0:infinity (constant rate)"
	print >> output, """        Malpha = 0   * 1: different alpha's for genes, 0: one alpha
        ncatG = 4   * # of categories in the dG, AdG, or nparK models of rates
        nparK = 0   * rate-class models. 1:rK, 2:rK&fK, 3:rK&MK(1/K), 4:rK&MK 

        nhomo = 0   * 0 & 1: homogeneous, 2: kappa for branches, 3: N1, 4: N2
        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 1   * (0,1,2): rates (alpha>0) or ancestral states

   Small_Diff = 7e-6
    cleandata = 0 * remove sites with ambiguity data (1:yes, 0:no)?
*        icode = 0  * (with RateAncestor=1. try "GC" in data,model=4,Mgene=4)
  fix_blength = 0  * 0: ignore, -1: random, 1: initial, 2: fixed
       method = 1  * Optimization method 0: simultaneous; 1: one branch a time
"""
	output.close()



		



def calculate_best_block_likelihoods(A, Genome_length, L):

	def get_block_likelihood(A,i,j):
		
		
		#using bionomial probability
		
		
		if i==j:
			return float("inf")
	
		n=float(A[j]+1-A[i])
		c=float(j+1-i)
		N=float(Genome_length)
		C=float(Number_of_SNPs)
		
		if (c/n)<=(C/N):
			return float("inf")
		
		try:
			part1=math.log((c/n),10)*c
		except (ValueError, ZeroDivisionError, OverflowError):
			print "part1err", n, c, N, C
		try:
			part2=math.log((((n-c)/n)),10)*(n-c)
		except (ValueError, ZeroDivisionError, OverflowError):
			#print "part2err", n, c, N, C
			part2=0
		try:
			part3=math.log((((C-c)/(N-n))),10)*(C-c)
		except (ValueError, ZeroDivisionError, OverflowError):
			#print "part3err", n, c, N, C
			part3=0
		try:
			part4=math.log(((((N-n)-(C-c))/(N-n))),10)*((N-n)-(C-c))
		except (ValueError, ZeroDivisionError, OverflowError):
			print "part4err", n, c, N, C
			part4=0
		
		likelihood=(part1+part2+part3+part4)*-1
		
		if likelihood<0:
			print A
			print n, c, N, C, part1, part2, part3, part4, i, j, A[i], A[j]
			sys.exit()
	
		
		return likelihood
	
	
		
	def calculate_block_likelihoods(A,x,y):
		p=[0]*((y+1)-x)
		ll=[0.0]*((y+1)-x)
		
		for i in xrange(y, (y+1)-x):
			p[i]=i
	
		
		for i in xrange(y,x-1,-1):
			p[i]=i
			ll[i]=[get_block_likelihood(A,i, p[i]+(L-1)), i, p[i]+(L-1)]
			
			while p[i]+(L-1)<y+1 and get_block_likelihood(A,i, p[i]+(L-1))>=get_block_likelihood(A,i, p[p[i]+1]+(L-1)):
				p[i]=p[p[i]+1]
				ll[i]=[get_block_likelihood(A, i, p[i]+(L-1)), i, p[i]+(L-1)]
	
		
		return ll
	
	
	
	Number_of_SNPs=len(A)
	#L=3
	if len(A)<L:
		print "Cannot calculate values for fewer than ", L, "SNP sites"
	
	starttime=time.clock()
	
	
	#Note that L in real SNP data refers to the number of SNPs
	
	l=b=r=y=Number_of_SNPs-L
	x=0
	
	#search for the most likely block using a method based on that of Goldwasser, Kao and Lu, Sept 2002. Fast Algorithms for Finding Maximum-Density Segments of a Sequence with Applications to Bioinformatics. Algorithms in Bioinformatics, Second international workshop, WABI 2002. Rome, Italy.
	
	ll=calculate_block_likelihoods(A,x,y)
	ll.sort()
	
	
	Minll=ll[0][:]
	
	#print Maxll, L
	
	numpermut=options.permutations
	pvaluecutoff=options.cutoff
	
	lls=[]
	
	pvalue=0.0
	#print A
	if Minll[0]!=float("inf"):
		for x in xrange(numpermut):
			try:
				B=random.sample(xrange(1,Genome_length+1), Number_of_SNPs)
			except StandardError:
				print Genome_length, Number_of_SNPs
			B.sort()
			#print B
			l=b=r=y=len(B)-L
			x=0
			ll=calculate_block_likelihoods(B,x,y)
			ll.sort()
			
			llinfo=[ll[0][0], (ll[0][2]+1)-ll[0][1], (B[ll[0][2]]+1)-B[ll[0][1]]]
			
			lls.append(llinfo)
			
			
			if ll[0][0]!=float("inf") and ll[0][0]<=Minll[0]:
				pvalue+=1
			if pvalue>=(numpermut*pvaluecutoff):
				#pvalue=float(numpermut)
				break
		pvalue=pvalue/numpermut
	else:
		pvalue=1.0
	
	Minll.append(pvalue)
	endtime=time.clock()
	
	#print Number_of_SNPs, endtime-starttime, pvalue
	
	if Minll[0]<0:
		
		print Minll[0], Minll[1], Minll[2], A[Minll[1]], A[Minll[2]], Minll[3], bettercount#, endtime
	if Minll[-1]<options.cutoff and (((A[Minll[2]]+1)-A[Minll[1]])/((Minll[2]+1)-Minll[1]))>Genome_length/len(A):
		
		lls.sort()
		print Genome_length, len(A), (Minll[2]+1)-Minll[1], (A[Minll[2]]+1)-A[Minll[1]], Minll, lls, get_block_likelihood(A,Minll[1],Minll[2])
	
	#ll structure [ll, startpos in A, endpos in A]
	
	return Minll









def find_recombinations_with_new_method(treeobject, snpposns, ungapped_to_gapped, nongaplen, node, daughter, nodenames, daughternames, downstreamtaxa):

	pvaluethreshold=options.cutoff
	
#	print snpposns
#	
#	print ungapped_to_gapped
	
	lennogaps=nongaplen
	totalsnps=float(len(snpposns))
	
	
	added=True
	
	blocks=[]
	
	lastcutoff=-1
	
	while added and totalsnps>=options.minsnps:
		#print node, daughter, snpposns[-10:]

		Best_block=calculate_best_block_likelihoods(snpposns, lennogaps, options.minsnps)
		
		
		
		snpcount=Best_block[2]+1-Best_block[1]
		added=False
		if snpcount>=options.minsnps:		
			
			blockll=Best_block[0]
			blockstart=snpposns[Best_block[1]]
			blockend=snpposns[Best_block[2]]
			genomestart=ungapped_to_gapped[Best_block[1]]
			genomeend=ungapped_to_gapped[Best_block[2]]
			pvalue=Best_block[3]
			
			if pvalue<pvaluethreshold:
			
				oldtotalsnps=totalsnps
				oldlennogaps=lennogaps
				
				lennogaps=lennogaps-((blockend+1)-blockstart)
				
				if lennogaps<0:
					print snpposns
					print oldlennogaps, lennogaps, blockll, Best_block[1], Best_block[2], blockstart, blockend, genomestart, genomeend	
					sys.exit()
				
				snpcount=(Best_block[2]+1)-Best_block[1]
				removedlength=(blockend+1)-blockstart
				tmpoutput=open("temp.output", "a")
				if removedlength>100000:
					print >> tmpoutput, oldlennogaps, lennogaps, blockll, Best_block[1], Best_block[2], blockstart, blockend, genomestart, genomeend
					print >> tmpoutput, snpposns
				tmpoutput.close()
				snpposns=snpposns[:Best_block[1]+1]+map(lambda x: x-removedlength, snpposns[Best_block[2]+1:])
				ungapped_to_gapped=ungapped_to_gapped[:Best_block[1]+1]+ungapped_to_gapped[Best_block[2]+1:]

				
				totalsnps=float(len(snpposns))

				blocks.append([genomestart, genomeend, blockll, snpcount, pvalue])
				added=True
	
	
	for block in blocks:
		downstreamnamelist=[]

		print >> tabout, "FT   misc_feature    "+str(block[0]+1)+".."+str(block[1]+1)
		if node==0 and options.outgroup!="" and options.outgroup!="None":
			
			editablealignment[options.outgroup]=editablealignment[options.outgroup][:block[0]]+"?"*(block[1]-block[0])+editablealignment[options.outgroup][block[1]:]
			
			
			print >> tabout, "FT                   /colour=3"
			print >> tabout, 'FT                   /taxa="'+options.outgroup+'"'
			print >> tabout, 'FT                   /node="'+nodenames[1]+'->'+options.outgroup+'"'	
			
		else:
			for taxon in downstreamtaxa:
				downstreamnamelist.append(sequencenames[taxon.split("_")[1]])
				editablealignment[sequencenames[taxon.split("_")[1]]]=editablealignment[sequencenames[taxon.split("_")[1]]][:block[0]]+"?"*(block[1]-block[0])+editablealignment[sequencenames[taxon.split("_")[1]]][block[1]:]

			#This bit adds each block to the rec.tab output file
	
			if treeobject.is_internal(daughter):
				print >> tabout, "FT                   /colour=2"
				print >> tabout, 'FT                   /taxa="'+', '.join(downstreamnamelist)+'"'
				print >> tabout, 'FT                   /node="'+nodenames[1]+'->'+daughternames[daughter][1]+'"'
						
			else:
				print >> tabout, "FT                   /colour=4"
				print >> tabout, 'FT                   /taxa="'+sequencenames[daughternames[daughter][1]]+'"'
				print >> tabout, 'FT                   /node="'+nodenames[1]+'->'+sequencenames[daughternames[daughter][1]]+'"'
			
		print >> tabout, 'FT                   /neg_log_likelihood='+str(block[2])
		print >> tabout, 'FT                   /SNP_count='+str(block[3])
		print >> tabout, 'FT                   /pvalue='+str(block[4])
	
	#print len(blocks)
	
	return blocks



################################
# Function to count tree nodes #
################################

def count_nodes(treeObject):


	def count_next_node(treeObject, node, nodecount):
		
		nodecount+=1
		daughters=treeObject.node(node).succ
		for daughter in daughters:
			nodecount=count_next_node(treeObject, daughter, nodecount)

		return nodecount
	
	nodecount=count_next_node(treeObject, treeObject.root, 0)
	return nodecount
	


################################################################
# Function to recurse across tree and run recombination search #
################################################################


def tree_recurse(node,treeobject,nodecount=0.0):
	
	
	if treeobject.is_internal(node):
		nodenames=(str(int(treeobject.node(node).get_data().support)),str(int(treeobject.node(node).get_data().support)))
	else:
		nodenames=(sequencenames[treeobject.node(node).get_data().taxon.split("_")[1]],treeobject.node(node).get_data().taxon.split("_")[1])
	
	daughters=treeobject.node(node).get_succ()
	
	
	#convert daughter node_ids to their names
	daughternames={}
	for daughter in daughters:
		if treeobject.is_internal(daughter):
			daughternames[daughter]=(str(int(treeobject.node(daughter).get_data().support)),str(int(treeobject.node(daughter).get_data().support)))
		else:
			daughternames[daughter]=(sequencenames[treeobject.node(daughter).get_data().taxon.split("_")[1]],treeobject.node(daughter).get_data().taxon.split("_")[1])
	

	
	for daughter in daughters:
		if not gaplocations.has_key(daughternames[daughter][0]):
			nodecount=tree_recurse(daughter,pamltree, nodecount)
	
	
	#identify missing sites in non-terminal nodes
	
	missingset=set(gaplocations[daughternames[daughters[0]][0]])
	for daughter in daughters[1:]:
		daughterset=set(gaplocations[daughternames[daughter][0]])
		comparisonset=missingset
		missingset=daughterset.intersection(comparisonset)
		
	
	#add to gaplocations dictionary
	
	gaplocations[nodenames[0]]=list(missingset)

	#create binary snp file for each parent-daughter comparison
	
	for daughter in daughters:
		nodecount+=1
	
		print "%.2f%% complete\r" % (100*(nodecount/total_number_of_nodes)),
		sys.stdout.flush()
		
		daughtergaps=set(gaplocations[daughternames[daughter][0]])
		nodegaps=set(gaplocations[nodenames[0]])
		
		gapposns=list(daughtergaps.union(comparisonset))
		
		gapposns.sort()
		
		gapposncount=0
	
		ungapped_to_gapped=[]
		SNPlocations_ungapped=[]
		numgaps=0
		numsnps=0
		nongaplen=0
		
		for x in xrange(pamlalignment.get_alignment_length()):
		
			if pamlsequences[nodenames[1]][x]!=pamlsequences[daughternames[daughter][1]][x] and pamlsequences[nodenames[1]][x] not in gaps_and_missing_set and pamlsequences[daughternames[daughter][1]][x] not in gaps_and_missing_set:
				while gapposncount<len(gapposns) and gapposns[gapposncount]<AllSNPlocations[x]:
					gapposncount+=1
				ungapped_to_gapped.append(AllSNPlocations[x])
				SNPlocations_ungapped.append(AllSNPlocations[x]-gapposncount)
				numsnps+=1
				#print SNP locations to snplocout file
				if treeobject.is_internal(daughter):
					print >> snplocout, str(AllSNPlocations[x]+1)+",node_"+str(nodenames[1])+"->node_"+str(daughternames[daughter][1])+","+pamlsequences[nodenames[1]][x]+","+pamlsequences[daughternames[daughter][1]][x]
					print >>snptabout, "FT   SNP             "+str(AllSNPlocations[x]+1)
					print >>snptabout, 'FT                   /node="'+str(nodenames[1])+'->'+str(daughternames[daughter][1])+'"'
					downstreamtaxa=treeobject.get_taxa(daughter)
					downstreamnames=[]
					for downtax in downstreamtaxa:
						downstreamnames.append(sequencenames[downtax.split("_")[1]])
					print >>snptabout, 'FT                   /taxa="'+' '.join(downstreamnames)+'"'
					print >>snptabout, 'FT                   /SNP="'+pamlsequences[nodenames[1]][x]+'->'+pamlsequences[daughternames[daughter][1]][x]+'"'
					print >>snptabout, 'FT                   /colour=1'
				else:
					print >> snplocout, str(AllSNPlocations[x]+1)+",node_"+str(nodenames[1])+"->"+str(sequencenames[daughternames[daughter][1]])+","+pamlsequences[nodenames[1]][x]+","+pamlsequences[daughternames[daughter][1]][x]
					print >>snptabout, "FT   SNP             "+str(AllSNPlocations[x]+1)
					print >>snptabout, 'FT                   /node="'+str(nodenames[1])+'->'+str(sequencenames[daughternames[daughter][1]])+'"'
					downstreamtaxa=treeobject.get_taxa(daughter)
					downstreamnames=[]
					for downtax in downstreamtaxa:
						downstreamnames.append(sequencenames[downtax.split("_")[1]])
					print >>snptabout, 'FT                   /taxa="'+' '.join(downstreamnames)+'"'
					print >>snptabout, 'FT                   /SNP="'+pamlsequences[nodenames[1]][x]+'->'+pamlsequences[daughternames[daughter][1]][x]+'"'
					print >>snptabout, 'FT                   /colour=1'
			
		
		nongaplen=alignment.get_alignment_length()-len(gapposns)
		
		if numsnps>options.minsnps:
			
			downstreamtaxa=treeobject.get_taxa(daughter)
			
			#print alignment.get_alignment_length(), nongaplen
			
			#use new approach to detect recombinations
			
			blocks=find_recombinations_with_new_method(treeobject, SNPlocations_ungapped, ungapped_to_gapped, nongaplen, node, daughter, nodenames, daughternames, downstreamtaxa)

					

	return nodecount


	
#####################################################################
# Function to create a tab file from a list of recombination blocks #
#####################################################################
	
def create_tabfile_from_blocks(blocks, treeobject, branch):	
	
	node=int(branch.split("_")[0])
	daughter=int(branch.split("_")[1])
	
	
	if treeobject.is_internal(node):
		nodenames=(str(int(treeobject.node(node).get_data().support)),str(int(treeobject.node(node).get_data().support)))
	else:
		nodenames=(sequencenames[treeobject.node(node).get_data().taxon.split("_")[1]],treeobject.node(node).get_data().taxon.split("_")[1])
	
	daughters=treeobject.node(node).get_succ()
	
	#convert daughter node_ids to their names
	if treeobject.is_internal(daughter):
		daughternames=(str(int(treeobject.node(daughter).get_data().support)),str(int(treeobject.node(daughter).get_data().support)))
	else:
		daughternames=(sequencenames[treeobject.node(daughter).get_data().taxon.split("_")[1]],treeobject.node(daughter).get_data().taxon.split("_")[1])

	downstreamtaxa=treeobject.get_taxa(daughter)
	
	for block in blocks:
		downstreamnamelist=[]
		
		if node==0:
			continue
			editablealignment[options.outgroup]=editablealignment[options.outgroup][:block[0]]+"-"*(block[1]-block[0])+editablealignment[options.outgroup][block[1]:]
			
			print >> tabout, "FT   misc_feature    "+str(block[0])+".."+str(block[1])
			print >> tabout, "FT                   /colour=3"
			print >> tabout, 'FT                   /note="'+options.outgroup+'"'
			print >> tabout, 'FT                   /node="node_'+nodenames[1]+'->node_'+options.outgroup+'"'					
			
		else:
			for taxon in downstreamtaxa:
				downstreamnamelist.append(sequencenames[taxon.split("_")[1]])
				#print editablealignment[sequencenames[taxon.split("_")[1]]][block[0]-10:block[1]+10]
				editablealignment[sequencenames[taxon.split("_")[1]]]=editablealignment[sequencenames[taxon.split("_")[1]]][:block[0]]+"-"*(block[1]-block[0])+editablealignment[sequencenames[taxon.split("_")[1]]][block[1]:]
				#print editablealignment[sequencenames[taxon.split("_")[1]]][block[0]-10:block[1]+10]

	
			#This bit adds each block to the rec.tab output file
	
#			if treeobject.is_internal(daughter):
#				print >> tabout, "FT   misc_feature    "+str(block[0])+".."+str(block[1])
#				print >> tabout, "FT                   /colour=2"
#				print >> tabout, 'FT                   /note="'+', '.join(downstreamnamelist)+'"'
#				print >> tabout, 'FT                   /node="node_'+nodenames[1]+'->node_'+daughternames[1]+'"'
#						
#			
#			else:
#				print >> tabout, "FT   misc_feature    "+str(block[0])+".."+str(block[1])
#				print >> tabout, "FT                   /colour=4"
#				print >> tabout, 'FT                   /note="'+sequencenames[daughternames[1]]+'"'
#				print >> tabout, 'FT                   /node="node_'+nodenames[1]+'->'+sequencenames[daughternames[1]]+'"'	
			
			
			
#			This bit prints different coloured features for different emissions in 3-state hmm
			
			if block[2]==1:
				print >> tabout, "FT   misc_feature    "+str(block[0])+".."+str(block[1])
				print >> tabout, "FT                   /colour=2"
				print >> tabout, 'FT                   /note="'+', '.join(downstreamnamelist)+'"'
				print >> tabout, 'FT                   /node="node_'+nodenames[1]+'->node_'+daughternames[1]+'"'
						
			
			else:
				print >> tabout, "FT   misc_feature    "+str(block[0])+".."+str(block[1])
				print >> tabout, "FT                   /colour=4"
				print >> tabout, 'FT                   /note="'+', '.join(downstreamnamelist)+'"'
				print >> tabout, 'FT                   /node="node_'+nodenames[1]+'->node_'+daughternames[1]+'"'

			
			
################
# Main program #
################		

if __name__ == "__main__":

	#Get command line arguments
	
	(options, args) = main()

	#Do some checking of the input files
	
	if options.alignment=="":
		DoError("No alignment file specified")
	elif options.tree=="" and not options.runtree:
		options.runtree=True
	elif options.tree!="" and options.runtree:
		print "!!!Warning: Treefile provided and option to create tree with RAxML selected. Using user tree. RAxML will not be run!!!"
	
	if options.prefix=="":
		prefix=options.alignment.split('.')[0]
	else:
		prefix=options.prefix
	
	if options.maxiterations<2 or options.maxiterations>1000:
		DoError("Maximum number of iterations (-i) must be between 2 and 1,000")
	
	if options.minsnps<2 or options.minsnps>1000:
		DoError("Minimum number of SNPs in a recombination block (-m) must be between 2 and 1,000")
	
	if options.permutations<0 or options.permutations>100000:
		DoError("Number of random repeats (-r) must be between 0 and 100,000")
	
	if options.cutoff<=0 or options.cutoff>=1:
		DoError("P-value cutoff (-c) must be between 0 and 1 (exclusive)")
	
	if options.cutoff<(1.0/options.permutations):
		DoError("P-value cutoff (-c) must be greater than 1/number of repeats (-r)")
		
	
	
	#Read the alignment file
		
#	print "Reading alignment file"
	
#	AlignIO.read(open(options.alignment, "rU"), "fasta")
#	
#	sys.exit()
	
	sys.stdout.flush()
	try:
		alignment=read_alignment(options.alignment)
		#alignment = AlignIO.read(open(options.alignment), "fasta")
	except ValueError:
		DoError("Cannot open alignment file "+options.alignment+". Is it in the correct format?")
		
	#Create a copy of the alignment that can be changed
	
	editablealignment={}
	for record in alignment:
		editablealignment[record.id]=Seq(str(record.seq).replace("?","N"))
	alignmentlength=len(editablealignment[editablealignment.keys()[0]])
	
	#Here we can start our iterations
	
	iteration=0
	newtree=Trees.Tree()
	

	while True: 
		iteration=iteration+1
		
		print "\nIteration", iteration
		sys.stdout.flush()
		oldtree=newtree
		
		
		
		
		
		#create newalignment from edited alignment for making the tree
		
		newalignment=Alignment(Gapped(IUPAC.unambiguous_dna))
		for record in alignment:
			newalignment.add_sequence(str(record.id),  str(editablealignment[record.id]))
			editablealignment[record.id]=record.seq
	
		#locate sites with a SNP
		
		if iteration==1:
			SNPlocations, gaplocations=Find_SNP_and_gap_locations(newalignment)
		else:
			SNPlocations=Find_SNP_locations(newalignment, AllSNPlocations)
			gaplocations=Allgaplocations.copy()

		
		
		#create alignment of just the SNP sites
		if options.wholedata:
			SNPalignment=newalignment
		else:
			SNPalignment=Create_SNP_alignment(newalignment, SNPlocations)
		
		#print a phylip file of the SNP alignment
		
		sequencenames={}
		convertnameback={}
		seqnametoindex={}
		
		handle = open("SNPS_"+prefix+".phy", "w")
		handleb = open(prefix+"_iteration"+str(iteration)+".aln", "w")	
		print >> handle, len(SNPalignment), SNPalignment.get_alignment_length()
		count=1
		for record in SNPalignment:
	
			name="seq"+str(count)
			
			sequencenames[name]=record.id
			convertnameback[record.id]=name
			seqnametoindex[name]=count-1
			
			print >> handle, name+"  "+record.seq
			print >> handleb, ">"+record.id
			print >> handleb, record.seq
	# add this back in to split the alignment into blocks of 60 bases
	#		for f in range(0,len(record.seq),60):
	#			if f+60< len(record.seq):
	#				print >> handle, record.seq[f:f+60]
	#			else:
	#				print >> handle, record.seq[f:]
				
			count=count+1
		handle.close()
		handleb.close()
		
		
		#the first time, make a copy of the SNP alignment file to be used in all future iterations of the PAML analyses
		
		if iteration==1:
			os.system("cp SNPS_"+prefix+".phy AllSNPS_"+prefix+".phy")
			AllSNPlocations=list(SNPlocations)
			Allgaplocations=gaplocations.copy()
			if options.maxiterations>1:
				treedistances=open(prefix+"_treedistances.txt", "w")
				if options.geodesic:
					print >> treedistances, "Tree1,Tree2,Recomb SNPs in 1,Recomb SNPs in 2,Recomb SNPs in both,RF,GD,Tree1 RAxML -ll, Tree2 RAxML -ll"
				else:
					print >> treedistances, "Tree1,Tree2,Recomb SNPs in 1,Recomb SNPs in 2,Recomb SNPs in both,RF,Tree1 RAxML -ll, Tree2 RAxML -ll"
		
				
		
		
		#run tree
		
		if (options.runtree or options.tree=="") or iteration>1:
			
			if os.path.isfile("RAxML_info.SNPS_"+prefix):
				print "Removing old RAxML files"
				os.system("rm RAxML_*.SNPS_"+prefix)
				sys.stdout.flush()
				
			print "Running tree with RAxML"
			sys.stdout.flush()
			os.system("/software/pathogen/external/applications/RAxML/RAxML-7.0.4/raxmlHPC -f d -s SNPS_"+prefix+".phy -m GTRGAMMA -n SNPS_"+prefix+" > "+prefix+"temp.tmp")
			options.tree="RAxML_result."+"SNPS_"+prefix
			
			#extract stats from raxml output files
			
			treestats=os.popen('grep "Inference\[0\]" RAxML_info.SNPS_'+prefix).read()
			
			alpha=float(treestats.split(":")[2].split()[0])

			treestats=os.popen('grep "Likelihood" RAxML_info.SNPS_'+prefix).read()
			
			negloglike=float(treestats.strip().split(":")[1].replace("-",""))

			print "RAxML -log likelihood =", negloglike
		else:
			alpha=options.alpha
			negloglike="Unknown"
	
		#Read tree file
		
		print "Reading tree file"
		sys.stdout.flush()
		
		try:
			tree_string = open(options.tree).read()
			#if (options.runtree or options.tree=="") or iteration>1:
			#	os.system("mv "+options.tree+" "+options.tree+"iteration"+str(iteration))
				
		except IOError:
			DoError("Cannot open tree file "+options.tree)
		tree = Trees.Tree(tree_string, rooted=True)
		
		
		if options.outgroup!="" and options.outgroup!="None":
			print "Rooting tree on", options.outgroup
			sys.stdout.flush()
			tree.root_with_outgroup(outgroup=convertnameback[options.outgroup])
		elif options.outgroup!="None":
			print "Midpoint rooting tree"
			midpoint_root(tree)
	
		newtree=tree
		treestring=tree_to_string(tree, False, True, True, True,collapse=True, cutoff=1.0/(len(SNPlocations)+1))
		
		for name in sequencenames.keys():
			treestring=treestring.replace(name+":", sequencenames[name]+":")
		handle = open(prefix+"_iteration"+str(iteration)+".tre", "w")
		print >> handle, treestring+";"
		handle.close()
		
		#print treestring.replace(sequencenames[name]+":", name+":")
		
		
#		print oldtree
#		print newtree
	
		
		if iteration>1:
			os.system("cat "+prefix+"_iteration"+str(iteration)+".tre "+prefix+"_iteration"+str(iteration-1)+".tre | sed 's/:0.0;/;/g' > tmptrees.tre")
			
			os.system('~sh16/hashRF/hashrf tmptrees.tre 2 > rfout.tmp')
			try:
				rfout=open('rfout.tmp', 'rU').readlines()
				rf=float(rfout[11].split()[1])
			except IOError:
				rf="err"
			if options.geodesic:
				os.system('java -jar ~sh16/geodemaps/geodeMAPS.jar -o gdout.tmp tmptrees.tre >  /dev/null 2>&1')
				try:
					gdout=open('gdout.tmp', 'rU').readlines()
					gd=float(gdout[0].split()[2])
				except IOError:
					gd="err"
				print >> treedistances, ","+str(rf)+","+str(gd)+","+str(oldloglike)+","+str(negloglike)
				os.system("rm tmptrees.tre gdout.tmp rfout.tmp")
			else:
				print >> treedistances, ","+str(rf)+","+str(oldloglike)+","+str(negloglike)
				os.system("rm tmptrees.tre rfout.tmp")
				treedistances.flush()
			if int(rf)==0:
				break
	
		if oldtree.is_identical(newtree):
			break
		
		oldloglike=negloglike
		
		treestring=tree.to_string(False, True, True, True)
		for name in sequencenames.keys():
			treestring=treestring.replace(sequencenames[name]+":", name+":")
		handle = open(prefix+".tre", "w")
		print >> handle, treestring+";"
		handle.close()
			
		
		
			
		#create baseml control file for paml
		
		print "Running PAML to reconstruct ancestral states"
		sys.stdout.flush()
	
		create_baseml_control_file("AllSNPS_"+prefix+".phy", prefix+".tre", alpha)
		
		#run paml
		
		os.system("/nfs/users/nfs_m/mh10/software/paml41/bin/baseml > "+prefix+"temp.tmp")

		
		#remove spaces from rst alignment (necessary to allow easier reading of tree and ancestral sequences
		
		os.system("sed 's/node #//g' rst > rstnew")
		os.system('grep -a -v "^$" rstnew > rstnew2')
		
		#extract the tree with all nodes numbered from PAML rst output file
		
		print "Reading PAML tree"
		sys.stdout.flush()

		negloglikefile=os.popen("tail -n 1 rub")

		negloglike=negloglikefile.read().strip().split()[1]

		print "PAML -log likelihood =", negloglike
		
		pamltreefile=os.popen('grep -a -A 1 "tree with node labels for Rod Page\'s TreeView" rstnew | tail -n 1')
		
		pamltreestring=pamltreefile.read().strip()
		
		
		#read paml tree into memory (note that node ids will be stored as branchlengths)
		
		pamltree=Trees.Tree(pamltreestring, rooted=True)
		
		#convert branchlengths (node ids) into support attributes
		
		pamltree.branchlength2support()
		
		#get the root node number from the paml tree (I think this might always be 0)
		
		rootnode=pamltree.root
		
		#extract alignment PAML rst output file
		
		print "Reading ancestral sequences from PAML"
		sys.stdout.flush()
		
		pamlalignstats=os.popen('grep -a -A 1 "List of extant and reconstructed sequences" rstnew2 | tail -n 1')
		
		try:
			n=int(pamlalignstats.read().split()[0])
		except ValueError:
			DoError("PAML analysis failed")
		
		pamlalignfile=os.popen('grep -a -A '+str(n+1)+' "List of extant and reconstructed sequences" rstnew2 | tail -n '+str(n+1))
		
		pamlalignment = AlignIO.read(pamlalignfile, "phylip")
		
		#Now we have the ancestral state reconstructions and tree in memory we need to make the input for the recombination detection program
		
		#first create a dictionary of all paml snp sequences
		
		pamlsequences={}
		for record in pamlalignment:
			pamlsequences[record.id]=record.seq
		
				
		#then run the recombination detection script recursively across the branches of the tree, and create tab file of blocks found
		
		print "Identifying recombinations on each branch of the tree at a p-value of", options.cutoff, "using", options.permutations, "random permutations"
		sys.stdout.flush()
		
		tabout=open(prefix+"_rec.tab","w")
		
		binsnplist={}
		binsnppositions={}
		blocks={}
		
		snplocout=open(prefix+"_SNPS_per_branch.csv","w")
		snptabout=open(prefix+"_SNPS_per_branch.tab","w")
		total_number_of_nodes=count_nodes(tree)
		print >> snplocout, "SNP_location,Branch,ancestral_base,Daughter_base"
		nodecount=tree_recurse(rootnode,pamltree)
		print "100.00% complete"
		snplocout.close()
		snptabout.close()
		tabout.close()
		

		
			
		all_in1butnot2=0
		all_in2butnot1=0
		all_inboth=0
		
		for record in newalignment:
			oldseq=record.seq
			newseq=editablealignment[record.id]
			
			#print oldseq, newseq
			
			in1butnot2=0
			in2butnot1=0
			inboth=0
			
			for location in AllSNPlocations:
				if newseq[location]=="?" and oldseq[location]=="?":
					inboth+=1
				elif oldseq[location]=="?":
					in1butnot2+=1
				elif newseq[location]=="?":
					in2butnot1+=1
			
			#print in1butnot2, in2butnot1, inboth
			all_in1butnot2+=in1butnot2
			all_in2butnot1+=in2butnot1
			all_inboth+=inboth
			
		if options.maxiterations>1:
			print >> treedistances, str(iteration)+","+str(iteration+1)+","+str(all_in1butnot2)+","+str(all_in2butnot1)+","+str(all_inboth),
		
		
		print "\n",
		
		if iteration>=options.maxiterations or all_in1butnot2+all_in2butnot1==0:
			break

	#Just print some info about why the loop has ended

	if oldtree.is_identical(newtree):
		print "trees identical at iteration", iteration
	elif all_in1butnot2+all_in2butnot1==0:
		print "recombination blocks identical at iteration", iteration
	else:
		print "reached iteration", iteration
#		if iteration>1:
#			os.system("cat "+options.tree+"iteration"+str(iteration)+" "+options.tree+"iteration"+str(iteration-1)+" > tmptrees.tre")
#			os.system('~sh16/hashRF/hashrf tmptrees.tre 2 > rfout.tmp')
#			rfout=open('rfout.tmp', 'rU').readlines()
#			rf=float(rfout[11].split()[1])
#			os.system('java -jar ~sh16/geodemaps_v0.2/geodeMAPS.jar -v tmptrees.tre > gdout.tmp')
#			gdout=open('gdout.tmp', 'rU').readlines()
#			gd=float(gdout[0].split()[2])
#			print >> treestats, ",".join([(iteration-1), iteration, rf, gd])
#			os.system("rm tmptrees.tre gdout.tmp rfout.tmp")
		
	sys.stdout.flush()
	
	if options.maxiterations>1:
		treedistances.close()
	
	#print final tree
	
	treestring=tree.to_string(False, True, True, True)
	for name in sequencenames.keys():
		treestring=treestring.replace(name+":", sequencenames[name]+":")
	
	
#	handle = open(prefix+"_Final.tre", "w")
#	print >> handle, treestring+";"
#	handle.close()

	os.system("mv "+prefix+"_iteration"+str(iteration)+".tre "+prefix+"_Final.tre")
	os.system("~sh16/scripts/Genome_Diagram.py -q taxa -t "+prefix+"_Final.tre -o "+prefix+"_Final_recomb "+prefix+"_rec.tab")
	
	os.system("rm "+prefix+"temp.tmp baseml.ctl rst rst1 2base.t mlb lnf rub rstnew rstnew2 RAxML_*.SNPS_"+prefix+" "+prefix+".tre")
	
	
	if options.bootstrap:
		os.system("~sh16/scripts/run_RAxML.py -a "+prefix+"_iteration"+str(iteration)+".aln -o "+prefix+" -M 5 -q long -w")
	
	
