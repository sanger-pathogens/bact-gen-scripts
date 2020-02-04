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
sys.path.extend(map(os.path.abspath, ['/nfs/pathogen/sh16_scripts/modules/']))
from Si_nexus import *

import pylab
import numpy

#Requires Biopython and Pysam


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-a", "--alignment", action="store", dest="alignment", help="alignment file name", default="")
	parser.add_option("-d", "--detectiontype", action="store", dest="detectiontype", type="choice", choices=["hmm", "movingwindow"], help="recombination detection method to use (choose from 'hmm', 'movingwindow') [default=%default]", default="movingwindow")
	parser.add_option("-i", "--iterations", action="store", dest="maxiterations", type="int", help="maximum number of iterations to run [default=%default]", default=3)
	parser.add_option("-m", "--minsnps", action="store", dest="minsnps", type="int", help="minimum number of snps required on branch to run recombination detection algorithm [default=%default]", default=2)
	parser.add_option("-o", "--outgroup", action="store", dest="outgroup", help="outgroup", default="")
	parser.add_option("-p", "--prefix", action="store", dest="prefix", help="prefix for output files", default="")
	parser.add_option("-R", "--RAxML", action="store_true", dest="runtree", help="run phylogeny with RAxML [default=%default]", default=False)
	parser.add_option("-t", "--tree", action="store", dest="tree", help="starting tree (optional)", default="")
	
	parser.add_option("-A", "--alpha", action="store", dest="alpha", type="float", help="alpha parameter for gamma distribution (starting tree is specified) [default=%default]", default=1)
	parser.add_option("-P", "--PAML", action="store_false", dest="runpaml", help="do not run PAML to reconstruct ancestral sequences [default=run PAML]", default=True)
	parser.add_option("-r", "--reference", action="store", dest="reference", help="reference embl file to add to final diagram (optional)", default="")
	
	return parser.parse_args()


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
			
			if base!='-' and base not in foundbases:# 
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



#######################################################################
# Function to detect recombination regions using double moving window #
#######################################################################


def detect_recombination_using_moving_windows(binsnps, treeobject, node, daughter, nodenames, daughternames, downstreamtaxa):
	
#	print "Under construction"
#
#	print "\nCalculating chi-squared for snps..."
#	sys.stdout.flush()
	
	#Calculate Chi-square looking for snp clustering. Window size=1000
	window=100
	totallen=len(binsnps)
	X2=[]
	pvalues=[]
	currwindowcount=0
	currposn=0
	X2converter={}
	
	snpposns=[]
	gapposns=[]
	nongapposns=[]
	
	for x, binsnp in enumerate(binsnps):
		if binsnp==1:
			snpposns.append(len(nongapposns))
			nongapposns.append(x)
		elif binsnp==2:
			gapposns.append(x)
		elif binsnp==0:
			nongapposns.append(x)
	
	lennogaps=len(nongapposns)
	totalsnps=len(snpposns)
	expectedsnps=float(totalsnps)/lennogaps
	expectedsnps=expectedsnps*window
	expectednonsnps=window-expectedsnps
	windowcounts=[0]*(lennogaps)
	
	
	
	while (expectedsnps<10 or expectednonsnps<10) and window<(totallen/1000):
		#print "Expected frequencies too low for chi-squared test using a window size of "+str(window)+" (",expectedsnps,"and", expectednonsnps,"). Increasing window size to "+str(window*10)
		expectedsnps=expectedsnps*10
		expectednonsnps=expectednonsnps*10
		window=window*10
	
	
	if expectedsnps<1:
		expectedsnps=1
		expectednonsnps=window-expectedsnps
	
#	print "\nCalculating moving window snp counts..."
#	sys.stdout.flush()
	
	
	for position in snpposns:
		
		y=position-(window/2)
		curposn=y
		while (y+window)>curposn:
			if curposn<0:
				actualposn=lennogaps+curposn
			elif curposn>(lennogaps-1):
				actualposn=curposn-lennogaps
			else:
				actualposn=curposn
			
			windowcounts[actualposn]=windowcounts[actualposn]+1
			curposn=curposn+1
	
	
	
	
	
	
#	pylab.Figure()
#	n, bins, patches = pylab.hist(windowcounts, bins=20)
#	pylab.title("Moving window SNP number frequencies") 
#	pylab.xlabel("Number of SNPs") 
#	pylab.ylabel("Frequency") 
#	pylab.show()
#	
#	
	
	
	
	
#	for position in gapposns:
#		
#		y=position-(window/2)
#		curposn=y
#		while (y+window)>curposn:
#			if curposn<0:
#				actualposn=totallen+curposn
#			elif curposn>(totallen-1):
#				actualposn=curposn-totallen
#			else:
#				actualposn=curposn
#			
#			windowcounts[actualposn][1]=windowcounts[actualposn][1]+1
#			curposn=curposn+1

#	print "\nCalculating chi-squared values..."
#	sys.stdout.flush()	

	
	for x in windowcounts:
		if not X2converter.has_key(x):
			X2converter[x]=[ float(((abs(x-expectedsnps)-0.5)**2)/expectedsnps) + float(((abs((window-x)-expectednonsnps)-0.5)**2)/expectednonsnps),0]
			if x<expectedsnps:
				X2converter[x][1]=1
			else:
				X2converter[x][1]=(1-chi2.cdf(X2converter[x][0],1))
		
		X2.append(X2converter[x][0])
		pvalues.append(X2converter[x][1])
		
#	print "\nPrinting chi-squared files..."
#	sys.stdout.flush()
	
	#output=open('moving_window_snps'+str(window)+'.plot','w')
	#print >> output, "# BASE MWSNPs X2 p"
#	for x, pvalue in enumerate(pvalues):
#		print >> output, pvalue
#		print >> outputb, X2[x]
#		print >> outputc, windowcounts[x]

	inblock="n"
	blocks=[]
	blockstart=0
	
	x=0
#	while x<len(pvalues):
#	#for x, pvalue in enumerate(pvalues):
#		#if count>0:
#		#	print >> output, nongapposns[x]+1, windowcounts[x], X2[x], pvalues[x]
#		pvalue=pvalues[x]
#		
#		#print x, pvalue, inblock
#		
#		if pvalue<=0.05 and inblock=="n":
#			#blockstart=nongapposns[x]+1
#			inblock="y"
#			
#			lastvalue=pvalues[x-(window/2)]
#			#print pvalues[x-(window/2):x]
#			for y, newvalue in enumerate(pvalues[x-(window/2):x]):
#				if newvalue<lastvalue:
#					print x-(window/2)+y, pvalues[x-(window/2)+y], newvalue, x+y, pvalues[y+x], x, pvalues[x]
#					blockstart=nongapposns[x+y]+1
#					x=x+y
#					break
#				lastvalue=newvalue
#		elif pvalue>0.05 and inblock=="y":
#			inblock="n"
#			lastvalue=pvalues[x+(window/2)]
#			for y, newvalue in enumerate(pvalues[x+(window/2):x:-1]):
#				if newvalue<lastvalue:
#					print x+(window/2)-y, pvalues[x+(window/2)-y], newvalue, x-y, pvalues[x-y], x, pvalues[x]
#					if nongapposns[(x-y)-1]+1>blockstart:
#						blocks.append([blockstart,nongapposns[(x-y)-1]+1])
#					#x=x-y
#					break
#				lastvalue=newvalue
#			#blocks.append([blockstart,nongapposns[x-1]+1])



	while x<len(pvalues):
	#for x, pvalue in enumerate(pvalues):
		#if count>0:
		#	print >> output, nongapposns[x]+1, windowcounts[x], X2[x], pvalues[x]
		pvalue=pvalues[x]
		
		#print x, pvalue, inblock
		
		if pvalue<=0.05 and inblock=="n":
			#blockstart=nongapposns[x]+1
			
			loc=0
			
			#print binsnps[nongapposns[x]:nongapposns[x+(window/2)]]
			
			if x+(window/2)+1>len(nongapposns):
				position=-1
			else:
				position=x+(window/2)+1
			
			for y in binsnps[nongapposns[x]:nongapposns[position]]:
				if y==1:
					#print "start", x,x+loc,y,nongapposns[x], nongapposns[(x+loc)], loc, binsnps[nongapposns[x+loc]]
					blockstart=nongapposns[(x+loc)]+1
					x=loc+x
					inblock="y"
					break
				if y!=2:
					loc=loc+1
		elif pvalue>0.05 and inblock=="y":
			
			loc=0
			
			if x-(window/2)+1<0:
				position=0
			else:
				position=x-(window/2)+1
			
			#print binsnps[nongapposns[x]:nongapposns[x-(window/2)]:-1]
			for y in binsnps[nongapposns[x]:nongapposns[position]:-1]:
				if y==1:
					#print "end", x,x-loc,y,nongapposns[x], nongapposns[(x-loc)], loc, binsnps[nongapposns[x-loc]]
					if nongapposns[(x-loc)]+1>blockstart:
						blocks.append([blockstart,nongapposns[(x-loc)]+1])
						#x=x-loc
					inblock="n"
					break
				if y!=2:
					loc=loc+1
			
			if inblock=="y":
				print binsnps[nongapposns[x]:nongapposns[position]:-1]
				print "here"
				sys.exit()
			
			#blocks.append([blockstart,nongapposns[x-1]+1])
			
	
		x=x+1
		
	if inblock=="y":
		blocks.append([blockstart,nongapposns[-1]+1])
	#sys.exit()
	for block in blocks:
		downstreamnamelist=[]
		
		if node==0 and options.outgroup!="":
			
			editablealignment[options.outgroup]=editablealignment[options.outgroup][:block[0]]+"-"*(block[1]-block[0])+editablealignment[options.outgroup][block[1]:]
			
			print >> tabout, "FT   misc_feature    "+str(block[0])+".."+str(block[1])
			print >> tabout, "FT                   /colour=3"
			print >> tabout, 'FT                   /taxa="'+options.outgroup+'"'
			print >> tabout, 'FT                   /node="node_'+nodenames[1]+'->node_'+options.outgroup+'"'					
			
		else:
			for taxon in downstreamtaxa:
				downstreamnamelist.append(sequencenames[taxon.split("_")[1]])
				editablealignment[sequencenames[taxon.split("_")[1]]]=editablealignment[sequencenames[taxon.split("_")[1]]][:block[0]]+"-"*(block[1]-block[0])+editablealignment[sequencenames[taxon.split("_")[1]]][block[1]:]

	
			#This bit adds each block to the rec.tab output file
	
			if treeobject.is_internal(daughter):
				print >> tabout, "FT   misc_feature    "+str(block[0])+".."+str(block[1])
				print >> tabout, "FT                   /colour=2"
				print >> tabout, 'FT                   /taxa="'+', '.join(downstreamnamelist)+'"'
				print >> tabout, 'FT                   /node="node_'+nodenames[1]+'->node_'+daughternames[daughter][1]+'"'
						
			
			else:
				print >> tabout, "FT   misc_feature    "+str(block[0])+".."+str(block[1])
				print >> tabout, "FT                   /colour=4"
				print >> tabout, 'FT                   /taxa="'+sequencenames[daughternames[daughter][1]]+'"'
				print >> tabout, 'FT                   /node="node_'+nodenames[1]+'->'+sequencenames[daughternames[daughter][1]]+'"'


	
	return blocks
		
		
#	output.close()
#	outputb.close()
#	outputc.close()
	
	print "Done."
	sys.stdout.flush()
#	sys.exit()
	
	return

	

################################################################
# Function to recurse across tree and run recombination search #
################################################################


def tree_recurse(node,treeobject):
	
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
			tree_recurse(daughter,pamltree)
	
	
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
		binsnps=[0]*alignment.get_alignment_length()
		
		for x in gaplocations[nodenames[0]]:
			binsnps[x]=2
		
		for x in gaplocations[daughternames[daughter][0]]:
			binsnps[x]=2
		
		numsnps=0
		

		
		
		for x in range(0,pamlalignment.get_alignment_length()):
				
			
				
					
				if pamlsequences[nodenames[1]][x]=="?" or pamlsequences[daughternames[daughter][1]][x]=="?":
					binsnps[AllSNPlocations[x]]=2
				elif pamlsequences[nodenames[1]][x]!=pamlsequences[daughternames[daughter][1]][x]:
					binsnps[AllSNPlocations[x]]=1
					numsnps=numsnps+1
					
					#print SNP locations to snplocout file
					if treeobject.is_internal(daughter):
						print >> snplocout, str(AllSNPlocations[x]+1)+",node_"+str(nodenames[1])+"->node_"+str(daughternames[daughter][1])+","+pamlsequences[nodenames[1]][x]+","+pamlsequences[daughternames[daughter][1]][x]
						print >>snptabout, "FT   SNP             "+str(AllSNPlocations[x]+1)
						print >>snptabout, 'FT                   /node="node_'+str(nodenames[1])+'->node_'+str(daughternames[daughter][1])+'"'
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
						print >>snptabout, 'FT                   /node="node_'+str(nodenames[1])+'->'+str(sequencenames[daughternames[daughter][1]])+'"'
						downstreamtaxa=treeobject.get_taxa(daughter)
						downstreamnames=[]
						for downtax in downstreamtaxa:
							downstreamnames.append(sequencenames[downtax.split("_")[1]])
						print >>snptabout, 'FT                   /taxa="'+' '.join(treeobject.get_taxa(daughter))+'"'
						print >>snptabout, 'FT                   /SNP="'+pamlsequences[nodenames[1]][x]+'->'+pamlsequences[daughternames[daughter][1]][x]+'"'
						print >>snptabout, 'FT                   /colour=1'
		
		snpposns={}
		gapposns=[]
		nongapposns={}
		
		y=0
		z=0
		
		
		#create files to convert positions back once the 2s have been removed from binsnps
		
		for x, binsnp in enumerate(binsnps):
			if binsnp==1:
				snpposns[z]=x
				nongapposns[y]=x
				z=z+1
				y=y+1
			elif binsnp==2:
				gapposns.append(x)
			elif binsnp==0:
				nongapposns[y]=x
				y=y+1
		

		
		if numsnps>options.minsnps:
			#print
			downstreamtaxa=treeobject.get_taxa(daughter)
			
			#use moving window approach to detect recombinations
			
			if options.detectiontype=="movingwindow":
				
				blocks=detect_recombination_using_moving_windows(binsnps, treeobject, node, daughter, nodenames, daughternames, downstreamtaxa)
					

			else:
				binsnplist[str(node)+"_"+str(daughter)]=binsnps
				binsnppositions[str(node)+"_"+str(daughter)]=nongapposns
					




	
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
		
		if node==0 and options.outgroup!="":
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
		print "!!!Warning: Treefile provided and option to create tree with RAxML selected. Using user tree. RAxML will not be run. !!!"
	
	if options.prefix=="":
		prefix=options.alignment.split('.')[0]
	else:
		prefix=options.prefix
	
	#Read the alignment file
		
	print "Reading alignment file"
	sys.stdout.flush()
	try:
		alignment = AlignIO.read(open(options.alignment), "fasta")
	except ValueError:
		DoError("Cannot open alignment file "+options.alignment+". Is it in the correct format?")
		
	#Create a copy of the alignment that can be changed
	
	editablealignment={}
	for record in alignment:
		editablealignment[record.id]=record.seq
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
			newalignment.add_sequence(str(record.id),  str(editablealignment[record.id]).replace("N","-"))
			editablealignment[record.id]=record.seq

	
		#locate sites with a SNP
		if iteration==1:
			SNPlocations, gaplocations=Find_SNP_and_gap_locations(newalignment)
		else:
			SNPlocations=Find_SNP_locations(newalignment, AllSNPlocations)
			gaplocations=Allgaplocations.copy()
		
		
		#create alignment of just the SNP sites
		
		SNPalignment=Create_SNP_alignment(newalignment, SNPlocations)
		
		#print a phylip file of the SNP alignment
		
		sequencenames={}
		convertnameback={}
		seqnametoindex={}
		
		handle = open("SNPS_"+prefix+".phy", "w")
		print >> handle, len(SNPalignment), SNPalignment.get_alignment_length()
		count=1
		for record in SNPalignment:
	
			name="seq"+str(count)
			
			sequencenames[name]=record.id
			convertnameback[record.id]=name
			seqnametoindex[name]=count-1
			
			print >> handle, name+"  "+record.seq
	# add this back in to split the alignment into blocks of 60 bases
	#		for f in range(0,len(record.seq),60):
	#			if f+60< len(record.seq):
	#				print >> handle, record.seq[f:f+60]
	#			else:
	#				print >> handle, record.seq[f:]
				
			count=count+1
		handle.close()
		
		
		#the first time, make a copy of the SNP alignment file to be used in all future iterations of the PAML analyses
		
		if iteration==1:
			os.system("cp SNPS_"+prefix+".phy AllSNPS_"+prefix+".phy")
			AllSNPlocations=list(SNPlocations)
			Allgaplocations=gaplocations.copy()
		
		
		
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
		else:
			alpha=options.alpha
	
		#Read tree file
		
		print "Reading tree file"
		sys.stdout.flush()
		
		try:
			tree_string = open(options.tree).read()
		except IOError:
			DoError("Cannot open tree file "+options.tree)
		tree = Trees.Tree(tree_string, rooted=True)
		
		
		if options.outgroup!="":
			print "Rooting tree on", options.outgroup
			sys.stdout.flush()
			tree.root_with_outgroup(outgroup=convertnameback[options.outgroup])
		else:
			print "Midpoint rooting tree"
			midpoint_root(tree)
	
		newtree=tree
		
#		print oldtree
#		print newtree
	
		if oldtree.is_identical(newtree):
			break
	
		treestring=tree.to_string(False, True, True, True)
		for name in sequencenames.keys():
			treestring=treestring.replace(sequencenames[name]+":", name+":")
		
		handle = open(prefix+".tre", "w")
		print >> handle, treestring+";"
		handle.close()
			
		
		#If we have chosen to run paml
		
		if options.runpaml:
			
			#create baseml control file for paml
			
			print "Running PAML to reconstruct ancestral states"
			sys.stdout.flush()
		
			create_baseml_control_file("AllSNPS_"+prefix+".phy", prefix+".tre", alpha)
			
			#run paml
			
			os.system("/nfs/phd/mh10/software/paml41/bin/baseml > "+prefix+"temp.tmp")
		
		#remove spaces from rst alignment (necessary to allow easier reading of tree and ancestral sequences
		
		os.system("sed 's/node #//g' rst > rstnew")
		os.system('grep -v "^$" rstnew > rstnew2')
		
		#extract the tree with all nodes numbered from PAML rst output file
		
		print "Reading PAML tree"
		sys.stdout.flush()
		
		pamltreefile=os.popen('grep -A 1 "tree with node labels for Rod Page\'s TreeView" rstnew | tail -n 1')
		
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
		
		pamlalignstats=os.popen('grep -A 1 "List of extant and reconstructed sequences" rstnew2 | tail -n 1')
		
		try:
			n=int(pamlalignstats.read().split()[0])
		except ValueError:
			DoError("PAML analysis failed")
		
		pamlalignfile=os.popen('grep -A '+str(n+1)+' "List of extant and reconstructed sequences" rstnew2 | tail -n '+str(n+1))
		
		pamlalignment = AlignIO.read(pamlalignfile, "phylip")
		
		#Now we have the ancestral state reconstructions and tree in memory we need to make the input for the recombination detection program
		
		#first create a dictionary of all paml snp sequences
		
		pamlsequences={}
		for record in pamlalignment:
			pamlsequences[record.id]=record.seq
		
				
		#then run the recombination detection script recursively across the branches of the tree, and create tab file of blocks found
		
		print "Identifying recombinations on each branch of the tree"
		sys.stdout.flush()
		
		tabout=open(prefix+"_rec.tab","w")
		
		binsnplist={}
		binsnppositions={}
		blocks={}
		
		snplocout=open(prefix+"_SNPS_per_branch.csv","w")
		snptabout=open(prefix+"_SNPS_per_branch.tab","w")
		print >> snplocout, "SNP_location,Branch,ancestral_base,Daughter_base"
		tree_recurse(rootnode,pamltree)
		snplocout.close()
		snptabout.close()
		
		
		if options.detectiontype=="hmm":
			
			#concatenate all binsnplists to use as a training set for the hmm
			
			allbinsnps=[]
			allbinsnplist=[]
			for key in binsnppositions.keys():
				if key.split("_")[0]=="0":
					continue
				allbinsnps=allbinsnps+binsnplist[key]
				allbinsnplist.append([i for i in binsnplist[key] if i != 2])
			
			#count the total number of snps in allbinsnps to use for setting hmm parameters
			
			numsnps=allbinsnps.count(1)
			
			#remove all 2s (sites where the state is unknown) from allbinsnps
			
			allbinsnps = [i for i in allbinsnps if i != 2]

			#set the hmm parameters

			sigma = IntegerRange(0,2)
			
			#Discrete hmm with 2 emissions
			
			transitions = [[1-(float(len(binsnplist.keys()))/len(allbinsnps)), float(len(binsnplist.keys()))/len(allbinsnps)], [0.0001, 0.9999]]
			ebackbone=[1-(float(numsnps)/len(allbinsnps)),float(numsnps)/len(allbinsnps)]
			erecombination=[1-(float(numsnps*10)/len(allbinsnps)),float(numsnps*10)/len(allbinsnps)]
			emissions=[ebackbone, erecombination]
			pi=[0.9999,0.0001]
			
			#Discrete hmm with 3 emissions
			
#			transitions = [[1-(float(len(binsnplist.keys()))/len(binsnplist[key])), float(len(binsnplist.keys()))/len(binsnplist[key]), 0], [0.0001, 0.9998, 0.0001],[0, 0.1, 0.9]]
#			ebackbone=[1-(float(numsnps)/len(allbinsnps)),float(numsnps)/len(allbinsnps)]
#			erecombination=[1-(float(numsnps*10)/len(allbinsnps)),float(numsnps*10)/len(allbinsnps)]
#			ethirdstate=[1-(float(numsnps*5)/len(allbinsnps)),float(numsnps*5)/len(allbinsnps)]
#			emissions=[ebackbone, erecombination, ethirdstate]
#			pi=[0.9999,0.0001, 0.0]
			
			#create the hmm
			
			m = HMMFromMatrices(sigma, DiscreteDistribution(sigma), transitions, emissions, pi)
			
			print m
			
			initial_log = m.loglikelihood(EmissionSequence(sigma, allbinsnps))
			
			print "Untrained log-likelihood", initial_log

			#run a baum-welch optimisation on the training set
			print "Running Baum-Welch optimisation"
			sys.stdout.flush()
			
			m.baumWelch(EmissionSequence(sigma, allbinsnps), 100, 0.001)
			
			trained_log = m.loglikelihood(EmissionSequence(sigma, allbinsnps))
			
			print "Baum-Welch training improved the log-likelihood by ", initial_log - trained_log
			sys.stdout.flush()
			
			print m
			
			#for each node run a viterbi and create blocks for recombination regions
			
			print "Running viterbi decoding for each node on tree"
			blockdists={}
			for key in binsnplist.keys():
				
				binsnplist[key] = [i for i in binsnplist[key] if i != 2]
				
				#might need to create new transmission and emission probs for ref
				
				if key.split("_")[0]=="0":
					print "Reoptimising hmm for ref"
					continue
					numsnps=binsnplist[key].count(1)
					transitions = [[1-(float(len(binsnplist.keys()))/len(binsnplist[key])), float(len(binsnplist.keys()))/len(binsnplist[key])], [0.0001, 0.9999]]
					ebackbone=[1-(float(numsnps)/len(binsnplist[key])),float(numsnps)/len(binsnplist[key])]
					erecombination=[1-(float(numsnps*10)/len(binsnplist[key])),float(numsnps*10)/len(binsnplist[key])]
					emissions=[ebackbone, erecombination]
					pi=[0.9999,0.0001]
							
					m2=HMMFromMatrices(sigma, DiscreteDistribution(sigma), transitions, emissions, pi)
					m2.baumWelch(EmissionSequence(sigma, binsnplist[key]), 100, 0.001)
					v = m2.viterbi(EmissionSequence(sigma, binsnplist[key]))
				else:
					v = m.viterbi(EmissionSequence(sigma, binsnplist[key]))
				
				inblock="n"
				blocks[key]=[]
				blockstart=0
				currstate=0
				blockdists[key]=[]
				lastblockend=0
				firstblockstart=0
				totalblocklength=0
				
				for x, state in enumerate(v[0]):
						
					if state!=currstate and state>0 and inblock=="n":
						blockstart=binsnppositions[key][x]+1
						inblock="y"
						if lastblockend==0:
							firstblockstart=blockstart
						else:
							blockdists[key].append(blockstart-lastblockend)
					elif state!=currstate and inblock=="y":
						totalblocklength=totalblocklength+(binsnppositions[key][x-1]+1-blockstart)
						blocks[key].append([blockstart,binsnppositions[key][x-1]+1,currstate])
						inblock="n"
						if state!=0:
							blockstart=binsnppositions[key][x]+1
							inblock="y"
						else:
							lastblockend=binsnppositions[key][x-1]+1
					currstate=state
				if lastblockend!=0:
					blockdists[key].append((len(v[0])-lastblockend)+firstblockstart)
				
				create_tabfile_from_blocks(blocks[key], pamltree, key)

		
		tabout.close()
		
		
		
				
		
		
		
		print "\n",
		
		if iteration>=options.maxiterations:
			break

	#Just print some info about why the loop has ended

	if oldtree.is_identical(newtree):
		print "trees identical at iteration", iteration
	else:
		print "reached iteration", iteration
	sys.stdout.flush()
	
	#print final tree
	
	treestring=tree.to_string(False, True, True, True)
	for name in sequencenames.keys():
		treestring=treestring.replace(name+":", sequencenames[name]+":")

	handle = open(prefix+"_Final.tre", "w")
	print >> handle, treestring+";"
	handle.close()
	if options.reference=="":
		os.system("~sh16/biopython_tests/Genome_Diagram.py -q taxa -t "+prefix+"_Final.tre -o "+prefix+"_Final_recomb "+prefix+"_rec.tab")
	else:
		os.system("~sh16/biopython_tests/Genome_Diagram.py -q taxa -t "+prefix+"_Final.tre -o "+prefix+"_Final_recomb "+options.reference+" "+prefix+"_rec.tab")
	
	os.system("rm "+prefix+"temp.tmp baseml.ctl rst rst1 2base.t mlb lnf rub rstnew rstnew2 RAxML_*.SNPS_"+prefix)
