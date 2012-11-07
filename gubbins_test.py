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

#Requires Biopython and Pysam


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-a", "--alignment", action="store", dest="alignment", help="alignment file name", default="")
	parser.add_option("-d", "--detectiontype", action="store", dest="detectiontype", type="choice", choices=["hmm", "movingwindow"], help="recombination detection method to use (choose from 'hmm', 'movingwindow') [default=%default]", default="movingwindow")
	parser.add_option("-i", "--iterations", action="store", dest="maxiterations", type="int", help="maximum number of iterations to run [default=%default]", default=5)
	parser.add_option("-m", "--minsnps", action="store", dest="minsnps", type="int", help="minimum number of snps required on branch to run recombination detection algorithm [default=%default]", default=3)
	parser.add_option("-o", "--outgroup", action="store", dest="outgroup", help="outgroup", default="")
	parser.add_option("-p", "--prefix", action="store", dest="prefix", help="prefix for output files", default="")
	parser.add_option("-R", "--RAxML", action="store_true", dest="runtree", help="run phylogeny with RAxML [default=%default]", default=False)
	parser.add_option("-w", "--wholedata", action="store_true", dest = "wholedata", help="Rerun final phylogeny on whole data (rather than just SNP sites - better but slooooower", default=False)
	parser.add_option("-b", "--bootstrap", action="store_true", dest = "bootstrap", help="bootstrap final phylogeny", default=False)
	parser.add_option("-t", "--tree", action="store", dest="tree", help="starting tree (optional)", default="")
	parser.add_option("-v", "--pvalue", action="store", dest="pvalue_cutoff", type="float", help="[value cutoff for chi-squared tests [default=%default]", default=0.05)
	parser.add_option("-A", "--alpha", action="store", dest="alpha", type="float", help="alpha parameter for gamma distribution (starting tree is specified) [default=%default]", default=1)
	parser.add_option("-P", "--PAML", action="store_false", dest="runpaml", help="do not run PAML to reconstruct ancestral sequences [default=run PAML]", default=True)
	parser.add_option("-x", "--oldPAML", action="store_true", dest="usepreviouspamlrun", help="start from previous paml data", default=False)
	parser.add_option("-g", "--geodesic", action="store_true", dest="geodesic", help="calculate geodesic distance between trees", default=False)
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



#######################################################################
# Function to detect recombination regions using double moving window #
#######################################################################


def detect_recombination_using_moving_windows_old(binsnps, treeobject, node, daughter, nodenames, daughternames, downstreamtaxa):
	
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

	#if not treeobject.is_internal(daughter) and sequencenames[daughternames[daughter][1]]=="4232_8_2":
	if treeobject.is_internal(daughter):
		output=open(daughternames[daughter][1]+".plot","w")
	else:
		output=open(sequencenames[daughternames[daughter][1]]+".plot","w")
	print >> output, "# BASE VAL"
	y=0
	x=0
	for x, pvalue in enumerate(pvalues):
		while y<nongapposns[x]:
			print >> output, y+1, 0, 1
			y+=1
		print >> output, int(nongapposns[x])+1, pvalue, 0
		y+=1
	output.close()
	
	#print len(binsnps), len(pvalues), len(nongapposns)
	x=0
	y=0
	while x<len(pvalues):
	#for x, pvalue in enumerate(pvalues):
		#if count>0:
		#	print >> output, nongapposns[x]+1, windowcounts[x], X2[x], pvalues[x]
		pvalue=pvalues[x]
		
		#print x, pvalue, inblock
		
		if pvalue<=options.pvalue_cutoff and inblock=="n":
			#blockstart=nongapposns[x]+1
			
			loc=0
			
			#print binsnps[nongapposns[x]:nongapposns[x+(window/2)]]
			
#			if x-(window/2)+1<0:
#				startposition=0
#			else:
#				startposition=x-(window/2)+1
#			if x+(window/2)+1>len(nongapposns):
#				endposition=-1
#			else:
#				endposition=x+(window/2)+1
			
			if x+window>len(nongapposns):
				position=-1
			else:
				position=x+window
			
			
			for y in binsnps[nongapposns[x]:nongapposns[position]]:
				if y==1:
					#print "start", x,x+loc,y,nongapposns[x], nongapposns[(x+loc)], loc, binsnps[nongapposns[x+loc]]
					blockstart=nongapposns[(x+loc)]+1
					x=loc+x
					inblock="y"
					break
				if y!=2:
					loc=loc+1
		elif pvalue>options.pvalue_cutoff and inblock=="y":
			
			loc=0
			
#			if x-(window/2)+1<0:
#				endposition=0
#			else:
#				endposition=x-(window/2)+1
#			if x+(window/2)+1>len(nongapposns):
#				startposition=-1
#			else:
#				startposition=x+(window/2)+1
				
			if x-window<0:
				position=0
			else:
				position=x-window
			
#			if not treeobject.is_internal(daughter) and sequencenames[daughternames[daughter][1]]=="4232_8_2":
#				print binsnps[nongapposns[x]:nongapposns[position]:-1]
#			if blockstart==1142515:
#				print pvalue, pvalues[x-5:x+5]
#				print binsnps[nongapposns[x]:nongapposns[position]:-1]
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
#				print "here"
#				sys.exit()
				if nongapposns[(x-loc)]+1>blockstart:
					blocks.append([blockstart,nongapposns[(x-loc)]+1])
				inblock="n"			
			#blocks.append([blockstart,nongapposns[x-1]+1])
			
	
		x=x+1
		
	if inblock=="y":
		blocks.append([blockstart,nongapposns[-1]+1])
	#sys.exit()
	for block in blocks:
		downstreamnamelist=[]
		
		snpcount=0
		Ncount=0
		for x in binsnps[block[0]:block[1]]:
			if x==1:
				snpcount+=1
			elif x==2:
				Ncount+=1
		
		if snpcount<options.minsnps:
			continue
		
#		if daughternames[daughter][1]=="11":
#			print block, binsnps[block[0]-5:block[1]+5], binsnps[block[0]:block[1]]
		
		if node==0:
			
			editablealignment[options.outgroup]=editablealignment[options.outgroup][:block[0]]+"-"*(block[1]-block[0])+editablealignment[options.outgroup][block[1]:]
			
			print >> tabout, "FT   misc_feature    "+str(block[0])+".."+str(block[1])
			print >> tabout, "FT                   /colour=3"
			print >> tabout, 'FT                   /taxa="'+options.outgroup+'"'
			print >> tabout, 'FT                   /node="'+nodenames[1]+'->'+options.outgroup+'"'
			print >> tabout, 'FT                   /SNP_count='+str(snpcount)			
			
		else:
			for taxon in downstreamtaxa:
				downstreamnamelist.append(sequencenames[taxon.split("_")[1]])
				editablealignment[sequencenames[taxon.split("_")[1]]]=editablealignment[sequencenames[taxon.split("_")[1]]][:block[0]]+"-"*(block[1]-block[0])+editablealignment[sequencenames[taxon.split("_")[1]]][block[1]:]

	
			#This bit adds each block to the rec.tab output file
	
			if treeobject.is_internal(daughter):
				print >> tabout, "FT   misc_feature    "+str(block[0])+".."+str(block[1])
				print >> tabout, "FT                   /colour=2"
				print >> tabout, 'FT                   /taxa="'+', '.join(downstreamnamelist)+'"'
				print >> tabout, 'FT                   /node="'+nodenames[1]+'->'+daughternames[daughter][1]+'"'
				print >> tabout, 'FT                   /SNP_count='+str(snpcount)	
				print >> tabout, 'FT                   /N_count='+str(Ncount)
						
			
			else:
				print >> tabout, "FT   misc_feature    "+str(block[0])+".."+str(block[1])
				print >> tabout, "FT                   /colour=4"
				print >> tabout, 'FT                   /taxa="'+sequencenames[daughternames[daughter][1]]+'"'
				print >> tabout, 'FT                   /node="'+nodenames[1]+'->'+sequencenames[daughternames[daughter][1]]+'"'
				print >> tabout, 'FT                   /SNP_count='+str(snpcount)	
				print >> tabout, 'FT                   /N_count='+str(Ncount)	


	
	return blocks
		
		
#	output.close()
#	outputb.close()
#	outputc.close()
	
	print "Done."
	sys.stdout.flush()
#	sys.exit()
	
	return

def get_N_and_C(binsnps):
	N=0.0
	C=0.0
	for x in binsnps:
		if x==0:
			N+=1
		elif x==1:
			C+=1
			N+=1
	return N, C

def get_block_likelihood(start, end, binsnps, N, C):

	

	n=0.0
	c=0.0
	for x in binsnps[start:end+1]:
		if x==0:
			n+=1
		elif x==1:
			c+=1
			n+=1

	#print start, end, c, n, C, N
	part1=math.log((c/n),10)*c
	if n-c==0:
		part2=0
	else:
		part2=math.log((((n-c)/n)),10)*(n-c)
	if C-c==0:
		part3=0
	else:
		part3=math.log((((C-c)/(N-n))),10)*(C-c)
	if ((N-n)-(C-c))==0:
		part4=0
	else:
		part4=math.log(((((N-n)-(C-c))/(N-n))),10)*((N-n)-(C-c))
	
	likelihood=(part1+part2+part3+part4)*-1
	
	#print start, end, c, n, C, N, likelihood
	
	return likelihood





def reduce_factorial(l,i):
	
	factorial=math.log(1.0,10)

	for x in range(int(l)-(int(i)-1),int(l)+1):
		#print x, factorial
		factorial=factorial+math.log(x,10)
	
	
	return factorial




def get_blocks_from_windowcounts(binsnps, windowcounts, cutoff, nongapposns):

	tempblocks=[]
	blockstart=0
	inblock=False
	x=0
	y=0
	#print len(windowcounts), len(nongapposns)
	while x<len(windowcounts):
	#for x, pvalue in enumerate(pvalues):
		#if count>0:
		#	print >> output, nongapposns[x]+1, windowcounts[x], X2[x], pvalues[x]
		value=windowcounts[x]
		
		#print x, pvalue, inblock
		
		if value>cutoff and not inblock:
			
			blockstart=nongapposns[x]+1
			inblock=True
		elif value<=cutoff and inblock:
				
			tempblocks.append([0.0,blockstart,nongapposns[x]+1])
			inblock=False
	
		x=x+1
		
	if inblock:
		tempblocks.append([0.0,blockstart,nongapposns[-1]])
	
	#print tempblocks
	
	#print len(binsnps), nongapposns[-1], len(windowcounts)
	
	#Trim blocks using likelihood
	
	
	N,C=get_N_and_C(binsnps)
	
	newblocks=[]
	
	for block in tempblocks:
		
		
		#trim to first and last SNPs in block
		start=int(block[1])
		end=int(block[2])
		
		n,c=get_N_and_C(binsnps[start:end+1])
		if c<options.minsnps:
			continue


##		print start, int(block[1])
##		
##		print binsnps[int(block[1]):start: -1]
		
		
		if binsnps[start]!=1:
		
			for x, basetype in enumerate(binsnps[start: end]):
			
				if basetype==1:
					start=start+x
					break
					
		if binsnps[end]!=1:
			for x, basetype in enumerate(binsnps[end:start: -1]):
	
				if basetype==1:
					end=end-x
					break
		
		
		old_like=get_block_likelihood(start, end, binsnps, N, C)
		
		
		#print "initial likelihood=", old_like, "start=", start, "end=", end
		
		newstart=start
		oldstart=start
		for x, basetype in enumerate(binsnps[start+1: end]):
			
			if basetype==1:
				new_like=get_block_likelihood(x+start+1, end, binsnps, N, C)
				#print new_like
				if new_like>old_like:
					start=newstart
					break
				else:
					newstart=start+x+1
					old_like=new_like
		
		
		if start==oldstart:
			newstart=start
			for x, basetype in enumerate(binsnps[start-1:0:-1]):
				
				if basetype==1:
					new_like=get_block_likelihood((start-1)-x, end, binsnps, N, C)
					#print new_like
					if new_like>old_like:
						start=newstart
						break
					else:
						newstart=(start-1)-x
						old_like=new_like
		
		
		newend=end
		oldend=end
		for x, basetype in enumerate(binsnps[ end-1: start: -1]):
			if basetype==1:
				new_like=get_block_likelihood(start, (end-1)-x, binsnps, N, C)
				if new_like>old_like:
					end=newend
					break
				else:
					newend=(end-1)-x
					old_like=new_like
		
		if end==oldend:
			newend=end
			for x, basetype in enumerate(binsnps[ end+1:]):
				if basetype==1:
					new_like=get_block_likelihood(start, end+1+x, binsnps, N, C)
					if new_like>old_like:
						end=newend
						break
					else:
						newend=end+1+x
						old_like=new_like
		
		
		
		#print "new likelihood=", old_like, "start=", start, "end=", end
		block[0]=old_like
		block[1]=start
		block[2]=end
		snpcount=0
		for x, snp in enumerate(binsnps[start:end+1]):
			if snp==1:
				snpcount+=1
		if snpcount>=options.minsnps:
			newblocks.append([block[0], block[1], block[2]])
	
	return newblocks



################################################################
# Function to detect recombination regions using moving window #
################################################################


def detect_recombination_using_moving_windows(binsnps, treeobject, node, daughter, nodenames, daughternames, downstreamtaxa):
	

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
	totalsnps=float(len(snpposns))
	
	if totalsnps<options.minsnps:
		return []
	
	elif lennogaps<100:
		print "branch has fewer than 100 nongap sites. Results would be poor... skipping"
		return [] 
	
	windowcounts=[0]*(int(lennogaps))
	
	window=int(float(lennogaps)/(totalsnps/10))
	
	if window<100:
		window=100
	elif window>10000:
		window=10000
	
	if window>float(lennogaps)/10:
		window=int(float(lennogaps)/10)
	
	#print lennogaps, totalsnps, window
	
	#print nodenames, daughternames
	
	#print float(lennogaps)/(totalsnps), (totalsnps/float(lennogaps))
	
	threshold=1-(0.05/(float(lennogaps)/(window/10)))
	
	
	#print float(lennogaps)/window, threshold
	
	#sys.exit()
	
	for position in snpposns:
		
		y=position-(window/2)
		curposn=int(y)
		while (y+window)>curposn:
			if curposn<0:
				actualposn=lennogaps+curposn
			elif curposn>(lennogaps-1):
				actualposn=curposn-lennogaps
			else:
				actualposn=curposn
			try:
				windowcounts[actualposn]=windowcounts[actualposn]+1
			except IndexError:
				print "Error: window position is out of range of alignment length"
				print curposn, actualposn, position, len(windowcounts), lennogaps, window
				sys.exit()
			curposn=curposn+1

	
	
	added=True
	
	blocks=[]
	
	lastcutoff=-1
	
	while added and totalsnps>=options.minsnps:
	
	
		window=int(float(lennogaps)/(totalsnps/10))
		
		if window<100:
			window=100
		elif window>10000:
			window=10000
			
		if window>float(lennogaps)/10:
			window=int(float(lennogaps)/10)
		#print lennogaps, totalsnps, window
		
		#print float(lennogaps)/(totalsnps), (totalsnps/float(lennogaps))
		
		threshold=1-(0.05/(float(lennogaps)/(window/10)))
		cutoff=0
		pvalue=0.0
		while pvalue<=threshold:
			#print window, cutoff, pvalue, threshold, totalsnps, reduce_factorial(float(window),cutoff), reduce_factorial(window,cutoff), reduce_factorial(cutoff,cutoff)
			part1=reduce_factorial(window,cutoff)-reduce_factorial(cutoff,cutoff)
			part2=math.log((float(totalsnps)/lennogaps),10)*cutoff
			part3=math.log((1.0-(float(totalsnps)/lennogaps)),10)*(window-cutoff)
			
			logthing=part1 + part2 + part3
			
			pvalue+=10**logthing
			
			cutoff+=1
			
			
		cutoff-=1
		
		
		if cutoff!=lastcutoff:
		
			newblocks=get_blocks_from_windowcounts(binsnps, windowcounts, cutoff, nongapposns)
				
				
			
		else:
			newblocks=newblocks[1:]
		
		lastcutoff=cutoff
		
		
		added=False
		if len(newblocks)>0:
			newblocks.sort()			
			
			snpcount=0
			oldtotalsnps=totalsnps
			oldlennogaps=lennogaps
			
			for x, snp in enumerate(binsnps[newblocks[0][1]:newblocks[0][2]+1]):
				if snp==1:
					snpcount+=1
					binsnps[newblocks[0][1]+x]=0
					totalsnps-=1
				lennogaps-=1
			
			
			if snpcount>=options.minsnps:
				numgaps=0
				for x in binsnps[newblocks[0][1]:newblocks[0][2]+1]:
					if x==2:
						numgaps+=1
				
				blocklen=(((newblocks[0][2]+1)-newblocks[0][1])-numgaps)
				x=0
				pvalue=0.0
				while x<snpcount:
				
					part1=reduce_factorial(blocklen,x)-reduce_factorial(x,x)
					part2=math.log((float(oldtotalsnps)/oldlennogaps),10)*x
					part3=math.log((1.0-(float(oldtotalsnps)/oldlennogaps)),10)*(blocklen-x)
				
					logthing=part1 + part2 + part3
					
					pvalue+=10**logthing
					x+=1
				pvalue=1.0-round(pvalue,10)
				pvaluethreshold=(0.05/float(lennogaps))
				#print blocklen, pvaluethreshold, snpcount, pvalue, oldtotalsnps, 1-threshold
				if pvalue<pvaluethreshold:
					blocks.append([newblocks[0][1], newblocks[0][2], newblocks[0][0], snpcount, pvalue])
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
	
	return blocks
		









	

################################################################
# Function to recurse across tree and run recombination search #
################################################################


def tree_recurse(node,treeobject):
	
	if node==treeobject.root:
		isroot=True
		daughters=treeobject.node(node).get_succ()
		node=daughters[0]
		if treeobject.is_internal(node):
			nodenames=(str(int(treeobject.node(node).get_data().support)),str(int(treeobject.node(node).get_data().support)))
		else:
			nodenames=(sequencenames[treeobject.node(node).get_data().taxon.split("_")[1]],treeobject.node(node).get_data().taxon.split("_")[1])
		
		daughters=daughters[1:]
	else:
		isroot=False
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
	

	if isroot:
		tree_recurse(node,pamltree)
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
				
			
				
					
				if pamlsequences[nodenames[1]][x]=="?" or pamlsequences[daughternames[daughter][1]][x]=="?":# or pamlsequences[nodenames[1]][x]=="N" or pamlsequences[daughternames[daughter][1]][x]=="N":
					try:
						binsnps[AllSNPlocations[x]]=2
					except StandardError:
						print x, len(AllSNPlocations), len(binsnps)
					
				elif pamlsequences[nodenames[1]][x]!=pamlsequences[daughternames[daughter][1]][x] and pamlsequences[daughternames[daughter][1]][x]!="N" and pamlsequences[nodenames[1]][x]!="N":
					binsnps[AllSNPlocations[x]]=1
					numsnps=numsnps+1
					
					#print SNP locations to snplocout file
					if treeobject.is_internal(daughter):
						print >> snplocout, str(AllSNPlocations[x]+1)+",node_"+str(nodenames[1])+"->node_"+str(daughternames[daughter][1])+","+pamlsequences[nodenames[1]][x]+","+pamlsequences[daughternames[daughter][1]][x]
						print >>snptabout, "FT   SNP             "+str(AllSNPlocations[x]+1)
						print >>snptabout, 'FT                   /node="'+str(nodenames[1])+'->'+str(daughternames[daughter][1])+'"'
						downstreamtaxa=treeobject.get_taxa(daughter)
						downstreamnames=[]
						for downtax in downstreamtaxa:
							downstreamnames.append(sequencenames[downtax.split("_")[1]])
						print >>snptabout, 'FT                   /taxa="'+', '.join(downstreamnames)+'"'
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
						print >>snptabout, 'FT                   /taxa="'+', '.join(downstreamnames)+'"'
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






def add_node_names_to_tree(tree, pamltree):
	
	paml_node_names={}
	
	def get_paml_nodes(node):
		if pamltree.is_internal(node):
			nodename=str(int(pamltree.node(node).get_data().support))		
		else:
			return
		
		
		
		downstreamtaxa=pamltree.get_taxa(node)
		downstreamnamelist=[]
		for taxon in downstreamtaxa:
			downstreamnamelist.append(taxon.split("_")[1])
		
		print node, nodename, downstreamnamelist
		downstreamnamelist.sort()
		paml_node_names[' '.join(downstreamnamelist)]=nodename
		daughters=pamltree.node(node).get_succ()
		for daughter in daughters:
			get_paml_nodes(daughter)
	
	get_paml_nodes(pamltree.root)

	print paml_node_names
	
	
	def add_nodenames_to_tree(node):
		if not tree.is_internal(node):
			return		

		downstreamtaxa=tree.get_taxa(node)
		downstreamnamelist=[]
		for taxon in downstreamtaxa:
			downstreamnamelist.append(taxon)
		
		downstreamnamelist.sort()
		print node, downstreamnamelist, paml_node_names[' '.join(downstreamnamelist)]
		daughters=tree.node(node).get_succ()
		node_data=tree.node(node).get_data()
		print node_data.taxon
		node_data.taxon=paml_node_names[' '.join(downstreamnamelist)]
		print node_data.taxon
		tree.node(node).set_data(node_data)
		for daughter in daughters:
			add_nodenames_to_tree(daughter)
	
	add_nodenames_to_tree(tree.root)
	
	return tree

	
			
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
		prefix=options.alignment.split("/")[-1].split('.')[0]
	else:
		prefix=options.prefix.split("/")[-1]
	
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
		
		
		
#		
#		
#		#create newalignment from edited alignment for making the tree
#		
#		newalignment=Alignment(Gapped(IUPAC.unambiguous_dna))
#		for record in alignment:
#			newalignment.add_sequence(str(record.id),  str(editablealignment[record.id]))
#			editablealignment[record.id]=record.seq
#	
#		#locate sites with a SNP
#		
#		if iteration==1:
#			SNPlocations, gaplocations=Find_SNP_and_gap_locations(newalignment)
#		else:
#			SNPlocations=Find_SNP_locations(newalignment, AllSNPlocations)
#			gaplocations=Allgaplocations.copy()
#
#		
#		
#		#create alignment of just the SNP sites
#		if options.wholedata:
#			SNPalignment=newalignment
#		else:
#			SNPalignment=Create_SNP_alignment(newalignment, SNPlocations)
#		
#		#print a phylip file of the SNP alignment
#		
		sequencenames={}
		convertnameback={}
		seqnametoindex={}
#		
#		handle = open("SNPS_"+prefix+".phy", "w")
#		handleb = open(prefix+"_iteration"+str(iteration)+".aln", "w")	
		SNPalignment=read_alignment(prefix+"_iteration"+str(iteration)+".aln")
#		print >> handle, len(SNPalignment), SNPalignment.get_alignment_length()
		count=1
		print SNPalignment
		for record in SNPalignment:
	
			name="seq"+str(count)
			
			sequencenames[name]=record.id
			convertnameback[record.id]=name
			seqnametoindex[name]=count-1
			count+=1

		print sequencenames
#			
#			print >> handle, name+"  "+record.seq
#			print >> handleb, ">"+record.id
#			print >> handleb, record.seq
#	# add this back in to split the alignment into blocks of 60 bases
#	#		for f in range(0,len(record.seq),60):
#	#			if f+60< len(record.seq):
#	#				print >> handle, record.seq[f:f+60]
#	#			else:
#	#				print >> handle, record.seq[f:]
#				
#			count=count+1
#		handle.close()
#		handleb.close()
#		
#		
#		#the first time, make a copy of the SNP alignment file to be used in all future iterations of the PAML analyses
#		
#		if iteration==1:
#			os.system("cp SNPS_"+prefix+".phy AllSNPS_"+prefix+".phy")
#			AllSNPlocations=list(SNPlocations)
#			Allgaplocations=gaplocations.copy()
#			if options.maxiterations>1:
#				treedistances=open(prefix+"_treedistances.txt", "w")
#				if options.geodesic:
#					print >> treedistances, "Tree1,Tree2,Recomb SNPs in 1,Recomb SNPs in 2,Recomb SNPs in both,RF,GD,Tree1 RAxML -ll, Tree2 RAxML -ll"
#				else:
#					print >> treedistances, "Tree1,Tree2,Recomb SNPs in 1,Recomb SNPs in 2,Recomb SNPs in both,RF,Tree1 RAxML -ll, Tree2 RAxML -ll"
#		
#				
#		
#		
#		#run tree
#		
#		if (options.runtree or options.tree=="") or iteration>1:
#			
#			if os.path.isfile("RAxML_info.SNPS_"+prefix):
#				print "Removing old RAxML files"
#				os.system("rm RAxML_*.SNPS_"+prefix)
#				sys.stdout.flush()
#				
#			print "Running tree with RAxML"
#			sys.stdout.flush()
#			os.system("/software/pathogen/external/applications/RAxML/RAxML-7.0.4/raxmlHPC -f d -s SNPS_"+prefix+".phy -m GTRGAMMA -n SNPS_"+prefix+" > "+prefix+"temp.tmp")
#			
#			os.system("mv RAxML_result."+"SNPS_"+prefix+" "+prefix+"_Initial.tre")
#			options.tree=prefix+"_Initial.tre"
#			
#			#extract stats from raxml output files
#			
#			treestats=os.popen('grep "Inference\[0\]" RAxML_info.SNPS_'+prefix).read()
#			
#			alpha=float(treestats.split(":")[2].split()[0])
#
#			treestats=os.popen('grep "Likelihood" RAxML_info.SNPS_'+prefix).read()
#			
#			negloglike=float(treestats.strip().split(":")[1].replace("-",""))
#
#			print "RAxML -log likelihood =", negloglike
#		else:
#			alpha=options.alpha
#			negloglike="Unknown"
#	
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
			
			daughters=tree.node(tree.root).succ
			if len(daughters)>2:
				print "too many daughters in rooted tree"
				sys.exit()
			daughter1_data=tree.node(daughters[0]).get_data()
			daughter2_data=tree.node(daughters[1]).get_data()
			daughter1_data.branchlength=daughter1_data.branchlength+daughter2_data.branchlength
			daughter2_data.branchlength=0
			tree.node(daughters[0]).set_data(daughter1_data)
			tree.node(daughters[1]).set_data(daughter2_data)
			
		elif options.outgroup!="None":
			print "Midpoint rooting tree"
			midpoint_root(tree)
			daughters=tree.node(tree.root).succ
			if len(daughters)>2:
				print "too many daughters in rooted tree"
				sys.exit()
			daughter1_data=tree.node(daughters[0]).get_data()
			daughter2_data=tree.node(daughters[1]).get_data()
			daughter1_data.branchlength=daughter1_data.branchlength+daughter2_data.branchlength
			daughter2_data.branchlength=0
			tree.node(daughters[0]).set_data(daughter1_data)
			tree.node(daughters[1]).set_data(daughter2_data)
	
		newtree=tree
#		treestring=tree_to_string(tree, False, True, True, True,collapse=True, cutoff=1.0/(len(SNPlocations)+1))
#		
#		for name in sequencenames.keys():
#			treestring=treestring.replace(name+":", sequencenames[name]+":")
#		handle = open(prefix+"_iteration"+str(iteration)+".tre", "w")
#		print >> handle, treestring+";"
#		handle.close()
#		
#		#print treestring.replace(sequencenames[name]+":", name+":")
#		
#		
##		print oldtree
##		print newtree
#	
#		
#		if iteration>1:
#			os.system("cat "+prefix+"_iteration"+str(iteration)+".tre "+prefix+"_iteration"+str(iteration-1)+".tre | sed 's/:0.0;/;/g' > tmptrees.tre")
#			
#			os.system('~sh16/hashRF/hashrf tmptrees.tre 2 > rfout.tmp')
#			try:
#				rfout=open('rfout.tmp', 'rU').readlines()
#				rf=float(rfout[11].split()[1])
#			except IOError:
#				rf="err"
#			if options.geodesic:
#				os.system('java -jar ~sh16/geodemaps/geodeMAPS.jar -o gdout.tmp tmptrees.tre >  /dev/null 2>&1')
#				try:
#					gdout=open('gdout.tmp', 'rU').readlines()
#					gd=float(gdout[0].split()[2])
#				except IOError:
#					gd="err"
#				print >> treedistances, ","+str(rf)+","+str(gd)+","+str(oldloglike)+","+str(negloglike)
#				os.system("rm tmptrees.tre gdout.tmp rfout.tmp")
#			else:
#				print >> treedistances, ","+str(rf)+","+str(oldloglike)+","+str(negloglike)
#				os.system("rm tmptrees.tre rfout.tmp")
#				treedistances.flush()
#			if int(rf)==0:
#				break
#	
#		if oldtree.is_identical(newtree):
#			break
#		
#		oldloglike=negloglike
#		
#		treestring=tree.to_string(False, True, True, True)
#		for name in sequencenames.keys():
#			treestring=treestring.replace(sequencenames[name]+":", name+":")
#		handle = open(prefix+".tre", "w")
#		print >> handle, treestring+";"
#		handle.close()
#			
#		
#		#If we have chosen to run paml
#		
#		if options.runpaml and (not options.usepreviouspamlrun or iteration>1):
#			
#			#create baseml control file for paml
#			
#			print "Running PAML to reconstruct ancestral states"
#			sys.stdout.flush()
#		
#			create_baseml_control_file("AllSNPS_"+prefix+".phy", prefix+".tre", alpha)
#			
#			#run paml
#			
#			os.system("/nfs/users/nfs_m/mh10/software/paml41/bin/baseml > "+prefix+"temp.tmp")
#
#		elif options.usepreviouspamlrun:
#			print "Using previous paml run"
#		
#		#remove spaces from rst alignment (necessary to allow easier reading of tree and ancestral sequences
#		
#		os.system("sed 's/node #//g' rst > rstnew")
#		os.system('grep -a -v "^$" rstnew > rstnew2')
#		
#		#extract the tree with all nodes numbered from PAML rst output file
#		
#		print "Reading PAML tree"
#		sys.stdout.flush()
#
#		negloglikefile=os.popen("tail -n 1 rub")
#
#		negloglike=negloglikefile.read().strip().split()[1]
#
#		print "PAML -log likelihood =", negloglike
		
		pamltreefile=os.popen('grep -a -A 1 "tree with node labels for Rod Page\'s TreeView" rstnew | tail -n 1')
		
		pamltreestring=pamltreefile.read().strip()
		
		
		#read paml tree into memory (note that node ids will be stored as branchlengths)
		
		pamltree=Trees.Tree(pamltreestring, rooted=True)
		
		#convert branchlengths (node ids) into support attributes
		
		pamltree.branchlength2support()
		
		#get the root node number from the paml tree (I think this might always be 0)
		
		rootnode=pamltree.root
		
		
		tree.display()
		pamltree.display()
		
		tree=add_node_names_to_tree(tree, pamltree)
		
		tree.display()
		treestring=tree_to_string(tree, support_as_branchlengths=False,branchlengths_only=True,plain=False,plain_newick=False,ladderize=None, collapse=False, cutoff=0.0, treename=False, comments=False, node_names=True)
		for name in sequencenames.keys():
			treestring=treestring.replace(name+":", sequencenames[name]+":")
		print treestring
		sys.exit()
		
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
	if options.reference=="":
		os.system("~sh16/scripts/Genome_Diagram.py -q taxa -t "+prefix+"_Final.tre -o "+prefix+"_Final_recomb "+prefix+"_rec.tab")
	else:
		os.system("~sh16/scripts/Genome_Diagram.py -q taxa -t "+prefix+"_Final.tre -o "+prefix+"_Final_recomb "+options.reference+" "+prefix+"_rec.tab")
	
	#os.system("rm "+prefix+"temp.tmp baseml.ctl rst rst1 2base.t mlb lnf rub rstnew rstnew2 RAxML_*.SNPS_"+prefix+" "+prefix+".tre SNPS_"+prefix+".phy")
	
	
	if options.bootstrap:
		os.system("~sh16/scripts/run_RAxML.py -a "+prefix+"_iteration"+str(iteration)+".aln -o "+prefix+" -M 5 -q long -w")
	
	
