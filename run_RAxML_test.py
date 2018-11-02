#!/usr/bin/env python
# !/usr/bin/env python
import string
import os, sys
from Bio import Seq
from Bio import AlignIO
from Bio.Align import Generic
from Bio.Alphabet import IUPAC, Gapped
from optparse import OptionParser, OptionGroup
import glob
from random import *

sys.path.extend(map(os.path.abspath, ['/nfs/pathogen/sh16_scripts/modules/']))
from Si_general import *
from Si_SeqIO import *
from Si_SNPs_temp import *

from socket import gethostname



import time



####################
# Set some globals #
####################


RAXML_DIR="/software/pathogen/external/applications/RAxML/RAxML-7.0.4/"


gap_and_missing=set(["-", "N", "?"])
missing_set=([ "N", "?"])



##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	
	group = OptionGroup(parser, "General options")
	group.add_option("-a", "--alignment", action="store", dest="alignment", help="Input alignment file name", default="", metavar="FILE")
	group.add_option("-s", "--snpsonly", action="store_true", dest="SNPonly", help="Only analyse SNP sites", default=False)
	parser.add_option("-p", "--proportion", action="store", dest="proportion", help="maximum proportion of Ns to allow in a column for it to be included in the analysis (i.e. ignore any sites with > than this proportion of Ns) [default=%default]", default=0.05)
	group.add_option("-e", "--exclude", action="store", dest="exclude", help="Exclude sequences from analysis if they are less than INT% covered [Default= %default]", default=50, type="float", metavar="int")
	group.add_option("-t", "--type", action="store", dest="analysistype", help="Data type (DNA or protein). [Choices = DNA, protein] [Default = %default]", default="DNA", type="choice", choices=["DNA","protein"])
#	group.add_option("-d", "--dmodel", action="store", dest="dmodel", help="Model of evolution to use (for DNA anlysis). [Default= %default]", default="GTRGAMMA", type="choice", choices=["GTRGAMMA","GTR", "GTR", "GTR", "GTR"])
	group.add_option("-m", "--model", action="store", dest="model", help="Model of evolution to use (for protein analysis only). [Choices = DAYHOFF, DCMUT, JTT, MTREV, WAG, RTREV, CPREV, VT, BLOSUM62, MTMAM, LG] [Default = %default]", default="WAG", type="choice", choices=["DAYHOFF","DCMUT", "JTT", "MTREV", "WAG", "RTREV", "CPREV", "VT", "BLOSUM62", "MTMAM", "LG"])
	group.add_option("-v", "--asrv", action="store", dest="asrv", help="Method of correction for among site rate variation [Choices = GAMMA, CAT, CAT_GAMMA, MIX] [Default = %default]", default="GAMMA", type="choice", choices=["GAMMA","CAT", "CAT_GAMMA", "MIX"])
	group.add_option("-i", "--pinvar", action="store_true", dest="pinvar", help="Use correction for proportion of invariant sites", default=False)
	group.add_option("-F", "--frequencies", action="store_true", dest="f", help="Use empirical base frequencies (protein models only)", default=False)
	group.add_option("-N", "--number", action="store", dest="number", help="Number of alternative ML runs on distinct starting trees [Default = %default]", default=1, type="int", metavar="INT")
	group.add_option("-d", "--distance", action="store_true", dest="distance", help="Compute pairwise ML distances (GAMMA-based models of ASRV only). Will use a parsimony tree unless you specify a tree with the -T option", default=False)
	group.add_option("-g", "--constraint", action="store", dest="constraint", help="File name of a multifurcating constraint tree. Does not have to contain all taxa", default="", metavar="FILE")
	group.add_option("-T", "--tree", action="store", dest="tree", help="File name of a user-specified starting tree", default="", metavar="FILE")
	
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "bootstrap options")
	group.add_option("-b", "--bootstrap", action="store", dest="bootstrap", help="Number of bootstrap replicates [Default = %default]", default=100, type="int", metavar="INT")
	group.add_option("-f", "--fast", action="store_true", dest="fast", help="Use RAxML fast bootstrap method (this will all run on one node, and therefore may be slower)", default=False)
	parser.add_option_group(group)
	group = OptionGroup(parser, "LSF options")
	group.add_option("-q", "--queue", action="store", dest="queue", help="LSF queue [Choices = normal, long, basement, hugemem (farm only)] [Default = %default]", default="normal", type="choice", choices=["normal","long", "basement", "hugemem"])
	group.add_option("-M", "--memory", action="store", dest="mem", help="Amount of memory required for analysis (Gb). 0 means do not require a memory limit, 10 means ask for at least 10Gb. [Default= %default]", default=0, type="int")
	
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "Output options")
	group.add_option("-o", "--output", action="store", dest="suffix", help="Suffix for output files", default="")
	group.add_option("-w", "--overwrite", action="store_true", dest="overwrite", help="Force overwrite files", default=False)
	parser.add_option_group(group)
	
	
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


#	if options.processors>cpu_count():
#		DoError('You do not have '+str(options.processors)+' processors')
#
	if options.alignment=='':
		DoError('No alignment file selected')
	elif not os.path.isfile(options.alignment):
		DoError('Cannot find file '+options.alignment)
	
	if options.constraint!='' and not os.path.isfile(options.constraint):
		DoError('Cannot find file '+options.constraint)

	if options.exclude<0 or options.exclude>=100:
		DoError('exclude percentage must be >=0 and <100')
	
	if options.number<1 or options.number>=100:
		DoError('number of alternative ML runs must be >=1 and <100')
#	elif options.number>1:
#		options.bootstrap=0
	
	if options.tree!="" and options.constraint!="":
		DoError('Cannot run with a constraint tree and starting tree')
	if options.constraint!="" and not os.path.isfile(options.constraint):
		DoError('Cannot find file '+options.constraint)
	if options.tree!="" and not os.path.isfile(options.tree):
		DoError('Cannot find file '+options.tree)
	
	if options.suffix=="":
		DoError('No output suffix selected')
	
	if options.distance and options.asrv!="GAMMA":
		print "Only GAMMA-based asrv models can be used when calculating distance. Setting asrv to GAMMA"
		options.asrv="GAMMA"
	
	userinput=""
	filelist=glob.glob('RAxML_*.'+options.suffix)+glob.glob('RAxML_*.ml_'+options.suffix)+glob.glob('RAxML_*.boot_'+options.suffix)
	
	if options.overwrite and len(filelist)>0:#Should improve this to check for specific analyses
		for f in filelist:
			os.remove(f)
	elif len(filelist)>0:
		print '\nRAxML files with extension '+options.suffix+' already exist!'
		while userinput not in ['y','n']:
			userinput=raw_input('Overwrite? (y/n): ')
			userinput.lower()
		if userinput=='y':
			for f in filelist:
				os.remove(f)
	elif options.mem>30 or options.mem<0:
		DoError('Memory requirement (-M) must be between 0 and 30Gb')





################
# Main program #
################		

if __name__ == "__main__":

	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)
	
	#make random name for temporary files
	chars = string.ascii_letters + string.digits
	tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))
	
	
	alignment={}
	
	#if the user only wants to analyse SNP sites
	if options.SNPonly:
		print "Finding SNP sites"
		SNPsites=set([])
		seq=[]
		
		ok_bases=set(['A','G','C','T'])
		foundref=False
		inrefseq=False
		for line in open(options.alignment, "rU"):
			line=line.strip()
			if len(line)>0 and line[0]==">":
				seq=''.join(seq).upper()
				if not foundref:
					refseq=seq
					seq=[]
					foundref=True
					inrefseq=True
					continue
				if inrefseq:
					refseq=seq
					inrefseq=False
					seq=[]
					continue
				if len(seq)==len(refseq):
					for base in xrange(len(refseq)):
						if base in SNPsites:
							continue
						if refseq[base]!=seq[base] and refseq[base] in ok_bases and seq[base] in ok_bases:
							SNPsites.add(base)
				else:
					print len(seq), len(refseq)
					print "Sequences are not all the same length"
					sys.exit()
			
				seq=[]
			elif foundref:
				seq.append(''.join(line.split()))
		
		seq=''.join(seq).upper()
		if len(seq)==len(refseq):
			for base in xrange(len(refseq)):
				if base in SNPsites:
					continue
				if refseq[base]!=seq[base] and refseq[base] in ok_bases and seq[base] in ok_bases:
					SNPsites.add(base)
		
		print "Found", len(SNPsites), "SNP sites"
		
		sorted_snps=list(SNPsites)
		sorted_snps.sort()
		
		
		print "Extracting SNP sites"
		name=""
		seq=[]
		for line in open(options.alignment, "rU"):
			line=line.strip()
			if len(line)>0 and line[0]==">":
				if name!="":
					seq=''.join(seq).upper()
					snpseq=[]
					for x in sorted_snps:
						snpseq.append(seq[x])
	
					alignment[name]=''.join(snpseq)
					seq=[]
					
				name=line.split()[0][1:]
			elif name!="":
				seq.append(''.join(line.split()))

		if name!="":
			seq=''.join(seq).upper()
			snpseq=[]
			for x in sorted_snps:
				snpseq.append(seq[x])
			
			alignment[name]=''.join(snpseq)
			
	else:		
		print "Reading alignment"
		name=""
		seq=[]
		for line in open(options.alignment, "rU"):
			line=line.strip()
			if len(line)>0 and line[0]==">":
				if name!="":
					alignment[name]=''.join(seq).upper()
					seq=[]
					
				name=line.split()[0][1:]
			elif name!="":
				seq.append(''.join(line.split()))

		if name!="":
			
			alignment[name]=''.join(seq).upper()

	if name=="":
		print "Found no sequences"
		sys.exit()
	
	alnlen=len(alignment[name])
	
	print "Creating final alignment of", len(alignment), "taxa and",  alnlen, "sites"
	
	#Remove sequences with too many unknown entries
	print "Filtering sequences with less than", str(options.exclude)+"% coverage"
	
	seqnames=[]
	toremove=[]
	for seq in alignment:
		if len(alignment[seq].replace("-","").replace("N",""))<(float(options.exclude)/100)*alnlen:
			print "Excluding", seq, "with", str((float(len(alignment[seq].replace("-","").replace("N","")))/alnlen)*100)+"% coverage"
			toremove.append(seq)
		else:
			seqnames.append(seq)
	
	for seq in toremove:
		del alignment[seq]
	
	taxacount=len(alignment)
	
	if options.proportion>0 and options.proportion<1:
	
		toremove=set([])
		cutoff=2
		if float(taxacount)*options.proportion>cutoff:
			cutoff=float(taxacount)*options.proportion
		
		print "Filtering sites with greater than", cutoff, "Ns"
		
		for x in xrange(alnlen):
			ncount=0.0
			for seq in seqnames:
				if alignment[seq][x] in gap_and_missing:
					ncount+=1
			if ncount>cutoff:
				toremove.add(x)
		print len(toremove), "bases will be removed"
	else:
		toremove=set([])

	
	#remove any previous RAxML runs with this tempname (this should never happen)
	filelist=glob.glob('RAxML_*'+tmpname)
	for file in filelist:
		os.remove(file)
	

	
	#print alignment to phylip format
	
	print "Creating final alignment of", taxacount, "taxa and",  alnlen-len(toremove), "sites"
	
	output_handle=open(tmpname+".phy", "w")
	
	print >> output_handle, taxacount, alnlen-len(toremove)
	for seq in alignment:
		newseq=[]
		for x, base in enumerate(alignment[seq]):
			if not x in toremove:
				newseq.append(base)
		print >> output_handle, seq, ''.join(newseq)
	
	
	output_handle.close()
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	
	#MEM(AA+GAMMA)    = (n-2) * m * (80 * 8) bytes
	#MEM(AA+CAT)           = (n-2) * m * (20 * 8) bytes
	#MEM(DNA+GAMMA) = (n-2) * m * (16 * 8) bytes
	#MEM(DNA+CAT)        = (n-2) * m * (4  * 8)  bytes
	#1GB=1073741824 bytes
	#1MB=1048576 bytes

	
	
	
	bsubcommand=["bsub"]
	bsubcommand.append("-q "+options.queue)
	if options.mem>0:
		bsubcommand.append("-M "+str(options.mem)+'000000 -R \'select[mem>'+str(options.mem)+'000] rusage[mem='+str(options.mem)+'000]\'')
		
#	bsubcommand.append("-o out.%J_"+options.suffix)
#	bsubcommand.append("-e err.%J_"+options.suffix)
	
	bsub=" ".join(bsubcommand)
	
	
	if options.analysistype=="DNA":
		model="GTR"+options.asrv
		if options.pinvar:
			model=model+"I"
	elif options.analysistype=="protein":
		
		model="PROT"+options.asrv
		if options.pinvar:
			model=model+"I"
		if options.model=="LG":
			model=model+"WAG"
		else:
			model=model+options.model
		if options.f:
			model=model+"F"
		if options.model=="LG":
			model=model+" -P ~sh16/data/LG.dat"

	
	
	#If user only wants to calculate distances
	if options.distance:
		print "Running RAxML to calculate pairwise ML distances with "+model+" model of evolution..."
		if options.tree:
			print "Using user-specified tree:", options.tree
			model=model+" -t "+options.tree
		else:
			print "Using parsimony tree"
		
		os.system(bsub+' -J "'+tmpname+'_dist" RAxML -f x -m '+model+' -s '+tmpname+'.phy -n '+options.suffix)
		os.system('bsub -w \'ended('+tmpname+'_dist)\' rm -rf \'*'+tmpname+'*\'')
	#If they want to make a tree rather than distances
	else:
		if options.fast:
			#Run the ML+fast bootstrap tree over LSF
			print "Running RAxML phylogeny with "+model+" model of evolution and "+str(options.bootstrap)+" fast bootstrap replicates..."
			
			if options.constraint:
				mlmodel=model+" -g "+options.constraint
			elif options.tree:
				mlmodel=model+" -t "+options.tree
			else:
				mlmodel=model
			
			os.system(bsub+" RAxML -f a -x "+str(randrange(1,99999))+" -p "+str(randrange(1,99999))+" -# "+str(options.bootstrap)+" -m "+mlmodel+" -s "+tmpname+".phy -n "+options.suffix)
		else:
			#Run the ML tree over LSF
			print "Running RAxML phylogeny with "+model+" model of evolution..."
			
			if options.constraint:
				mlmodel=model+" -g "+options.constraint
			elif options.tree:
				mlmodel=model+" -t "+options.tree
			else:
				mlmodel=model
			
			sys.stdout.flush()
			if options.number>1:
				os.system(bsub+' -J "'+tmpname+'_ml" RAxML -f d -s '+tmpname+'.phy -N '+str(options.number)+' -m '+mlmodel+' -n ml_'+options.suffix)
			else:
				os.system(bsub+' -J "'+tmpname+'_ml" RAxML -f d -s '+tmpname+'.phy -m '+mlmodel+' -n ml_'+options.suffix)
			outputname="RAxML_result."+tmpname
			
		if not options.fast and options.bootstrap>0:
			
			print "Running "+str(options.bootstrap)+" bootstrap replicates..."
			
			numjobs=options.bootstrap
			
			#Run the bootstraps over LSF
			
			bootstrapstring="echo 'RAxML -f d -b "+str(randrange(1,9999))+"${LSB_JOBINDEX} -# 1 -m "+model+" -s "+tmpname+".phy -n boot_"+tmpname+"_${LSB_JOBINDEX}'"
			host=gethostname()
			if host[:4]=="pcs4":
				os.system(bootstrapstring+' | '+bsub+' -J "'+tmpname+'_boot[1-'+str(numjobs)+']%15"')
			else:
				os.system(bootstrapstring+' | '+bsub+' -J "'+tmpname+'_boot[1-'+str(numjobs)+']"')
				
			
			#When all bootstrap replicates are finished, cat them into a single file
			os.system('echo "cat RAxML_bootstrap.boot_'+tmpname+'* > RAxML_bootstrap.boot_'+options.suffix+'" | bsub -w \'ended('+tmpname+'_boot)\' -J "'+tmpname+'_cat"')
			
			#When the ml and bootstrap replicates are complete, put the bootstrap numbers onto the nodes of the ml tree
			if options.number>1:
				bsubcommand="echo \'x=$(grep \'Best\' RAxML_info.ml_"+options.suffix+" | awk \"{print \$6}\" |tr -d \':\' ) && cp RAxML_result.ml_"+options.suffix+".RUN.${x} RAxML_result.ml_"+options.suffix+" &&  RAxML -f b -t RAxML_result.ml_"+options.suffix+" -z RAxML_bootstrap.boot_"+options.suffix+" -s "+tmpname+".phy -m "+model+" -n "+options.suffix+"'"
				print bsubcommand
				os.system(bsubcommand+' | bsub -J "'+tmpname+'_join" -w \'ended('+tmpname+'_ml) && ended('+tmpname+'_cat)\'')
			else:
				os.system('bsub -M 2000000 -R \'select[mem>2000] rusage[mem=2000]\' -J "'+tmpname+'_join" -w \'ended('+tmpname+'_ml) && ended('+tmpname+'_cat)\' RAxML -f b -t RAxML_result.ml_'+options.suffix+' -z RAxML_bootstrap.boot_'+options.suffix+' -s '+tmpname+'.phy -m '+model+' -n '+options.suffix)
					
	
			#Clean up all temporary files created	
			os.system('bsub -w \'ended('+tmpname+'_join)\' rm -rf \'*'+tmpname+'*\'')
		elif not options.fast and options.bootstrap==0:
			#Clean up all temporary files created	
			os.system('bsub -w \'ended('+tmpname+'_ml)\' rm -rf \'*'+tmpname+'*\'')
		

