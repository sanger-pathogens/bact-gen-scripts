#!/usr/bin/env python


##################
# Import modules #
##################

import string, re
import os, sys, getopt, random, math, time
from random import *
from optparse import OptionParser, OptionGroup
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pysam
import datetime
import shlex, subprocess


####################
# Set some globals #
####################

SAMTOOLS_DIR="/nfs/users/nfs_s/sh16/samtools-0.1.18/"
BCFTOOLS_DIR="/nfs/users/nfs_s/sh16/samtools-0.1.18/bcftools/"
SSAHA_DIR="/nfs/users/nfs_s/sh16/ssaha2_v2.5.1_x86_64/"
BWA_DIR=""
SMALT_DIR="/software/pathogen/external/apps/usr/bin/"
MY_SCRIPTS_DIR="/nfs/users/nfs_s/sh16/scripts/"
BREAKDANCER_DIR="/nfs/users/nfs_s/sh16/breakdancer-1.1_2011_02_21/"
PINDEL_DIR="/nfs/users/nfs_s/sh16/pindel/trunk/"
DINDEL_DIR="/nfs/users/nfs_s/sh16/dindel/binaries/"
DINDEL_SCRIPTS_DIR="/nfs/users/nfs_s/sh16/dindel/dindel-1.01-python/"


##########################
# Error message function #
##########################

def DoError(errorstring):
	print "\nError:", errorstring
	print "\nFor help use -h or --help\n"
	sys.exit()

########################################################
# Print commands to bash script file with log checking #
########################################################

def print_command_to_bashfile(command, dir_name, file_name, handle, loglines, keep):
	
	command=' '.join(command.split())
	
	
	if keep:
		for logline in loglines:
			if logline=="":
				continue
			try:
				logcommand=(' '.join(logline.strip().split()[:-1])).strip()
			except StandardError:
				logcommand=''
			try:
				logreturn=logline.strip().split()[-1]
			except StandardError:
				logreturn=''
			if logcommand!='' and command==logcommand and logreturn=="0":
#				print "Done this bit:", command
				logout=open(dir_name+"/"+file_name+".log", "a")
				print >> logout, command, "0"
				logout.close()
				return True

#	print command
#	sys.exit()
	print >> handle, command
	print >> handle, "status=$?"
	print >> handle, "command='"+command+"'"
	print >> handle, "echo '"+command+"' $status >>",dir_name+"/"+file_name+".log"
	print >> handle, 'if [ $status -ne "0" ]'
	print >> handle, " then exit status"
	print >> handle, "fi"
	return False
	
##############################
# Get command line arguments #
##############################

def get_user_options():
	usage = "usage: %prog [options] <list of fastq/bam files to be mapped>"
	version="%prog 1.0. Written by Simon Harris, Wellcome Trust Sanger Institute, 2011"
	parser = OptionParser(usage=usage, version=version)
	
	#do not allow arguments to be interspersed. e.g. -a -b arg1 agr2. MUST be -a arg1 -b arg2.
	parser.disable_interspersed_args()
	
	group = OptionGroup(parser, "Required Options")
	group.add_option("-r", "--reference", action="store", dest="ref", help="reference DNA sequence (in fasta or multi-fasta format)", default="", metavar="FILE")
	parser.add_option_group(group)

	group = OptionGroup(parser, "Mapping Options")
	group.add_option("-s", "--studyid", action="store", dest="studyid", help="ID of study to analyse", default=0, metavar="INT", type="int")
	group.add_option("-C", "--QC", action="store_true", dest="passed", help="If using a study ID, only use manually passed runs", default=False)
	group.add_option("-p", "--program", action="store", type="choice", dest="program", choices=["bwa","ssaha", "smalt", "BWA","SSAHA", "SMALT"], help="Mapping program to use (choose from bwa, ssaha or smalt) [default= %default]", default="smalt")
	#group.add_option("-l", "--length", action="store", dest="readlength", help="Read length [default= %default]", default=54, type="int", metavar="INT")
	group.add_option("-S", "--single", action="store_false", dest="pairedend", help="reads are single ended (not paired)", default=True)
	group.add_option("-i", "--maxinsert", action="store", dest="maxinsertsize", help="maximum insert size (ssaha and smalt only) [default= %default]", default=600, type="int", metavar="INT")
	group.add_option("-j", "--mininsert", action="store", dest="mininsertsize", help="minimum insert size (ssaha and smalt only) [default= %default]", default=0, type="int", metavar="INT")
	group.add_option("-F", "--properpairfilter", action="store_false", dest="properpairfilter", help="Do not filter out reads that don't map in a proper pair from the bam file (after indel detection) [default = Filter]", default=True)

	
	
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "Pseudosequence creation options")
	group.add_option("-X", "--dna", action="store_false", dest="pseudosequence", help="Do not create pseudosequences", default=True)
	group.add_option("-x", "--noref", action="store_false", dest="incref", help="Do not include reference in pseudosequence alignment", default=True)
	group.add_option("-I", "--indels", action="store_false", dest="indels", help="Do not include small indels in pseudosequence alignment (i.e. alignment will be the same length as the reference)", default=True)
	group.add_option("-q", "--quality", action="store", type="int", dest="quality", help="Minimum base call quality to call a SNP (see samtools help for more information) [default= %default]", default=50, metavar="INT")
	group.add_option("-Q", "--mapq", action="store", type="int", dest="mapq", help="Minimum mapping quality to call a SNP (see samtools help for more information) [default= %default]", default=30, metavar="INT")	
	group.add_option("-d", "--depth", action="store", dest="depth", help="Minimum number of reads matching SNP [default= %default]", default=4, type="int")
	group.add_option("-D", "--stranddepth", action="store", dest="stranddepth", help="Minimum number of reads matching SNP per strand [default= %default]", default=2, type="int")
	group.add_option("-B", "--BAQ", action="store_false", dest="BAQ", help="Turn off samtools base alignment quality option (BAQ) ", default=True)
	#parser.add_option("-q", "--quality", action="store", dest="quality", help="Minimum base quality [default= %default]", default=120, type="int")
	#group.add_option("-S", "--RMS", action="store", dest="RMS", help="Minimum root mean squared mapping quality [default= %default]", default=25, type="int")
	#parser.add_option("-Q", "--strandquality", action="store", dest="strandquality", help="Minimum per strand base quality [default= %default]", default=60, type="int")
	group.add_option("-R", "--ratio", action="store", dest="ratio", help="SNP/Mapping quality ratio cutoff [default= %default]", default=0.8, type="float")
	
	#group.add_option("-S", "--SNPquality", action="store", type="int", dest="snpquality", help="Minimum site mapping quality for SNP calling [default= %default]", default=90, metavar="INT")
	#group.add_option("-R", "--ratio", action="store", type="float", dest="ratio", help="SNP/site mapping quality ratio cutoff [default= %default]", default=0.75, metavar="FLOAT")
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "Output options")
	group.add_option("-e", "--embl", action="store", dest="embl", help="reference annotation (in true embl or genbank format) [required for dN/dS calculations and CDS naming]", default="", metavar="FILE")
	group.add_option("-o", "--output", action="store", dest="output", help="output file prefix", default="")
	group.add_option("-f", "--force", action="store_true", dest="force", help="force overwrite of output files", default=False)
	group.add_option("-g", "--plots", action="store_true", dest="plots", help="create mapping plots", default=False)
	group.add_option("-t", "--tabfiles", action="store_true", dest="tabfile", help="Create tabfile of snps", default=False)
	group.add_option("-a", "--align", action="store_true", dest="alnfile", help="Create snp alignment file (in phylip format)", default=False)
	group.add_option("-P", "--phylogeny", action="store_true", dest="raxml", help="Run phylogeny with RAxML", default=False)
	group.add_option("-m", "--model", action="store", dest="model", help="Model of evolution to use. [Default= %default]", default="GTRGAMMA", type="choice", choices=["GTRGAMMA","GTRGAMMAI", "GTRCAT", "GTRMIX", "GTRMIXI"])
	group.add_option("-b", "--bootstrap", action="store", dest="bootstrap", help="Number of bootstrap replicates (0 = do not run bootstrap). [Default= %default]", default=100, type="int", metavar="int")
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "General usage options")
	#group.add_option("-I", "--interactive", action="store_true", dest="interactive", help="Enter interactive menu system", default=True)
	group.add_option("-k", "--keep", action="store_true", dest="keep", help="If old mapping files are present, do not rerun them", default=False)
	group.add_option("-L", "--LSF", action="store_false", dest="LSF", help="Do not use LSF to parallelise analyses", default=True)
	group.add_option("-U", "--queue", action="store", dest="LSFQ", help="LSF queue to submit to. [Default= %default]", default="normal", type="choice", choices=["normal","long", "basement", "hugemem"])
	group.add_option("-M", "--memory", action="store", dest="mem", help="Amount of memory required for analysis (Gb). [Default= %default]", default=2, type="int")#change to be able to specify amount of memory to use per analysis?
	group.add_option("-n", "--nodes", action="store", dest="nodes", help="Maximum number of jobs to run on nodes in parallel. [Default= %default]", default=30, type="int")#change to be able to specify amount of memory to use per analysis?
	group.add_option("-y", "--dirty", action="store_true", dest="dirty", help="Do not clean up temporary files. [Default= %default]", default=False)#change to be able to specify amount of memory to use per analysis?
	parser.add_option_group(group)
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):

	if options.output=='':
		options.output=options.ref.split("/")[-1].split(".")[0]+"_"+options.program


	while options.force==False and os.path.isfile(options.output+".aln"):
		outopt=""
		outopt=raw_input('\nOutput files with chosen prefix already exist.\n\nWould you like to overwrite (o), choose a new output file prefix (n) or quit (Q): ')
		if outopt=='Q':
			sys.exit()
		elif outopt=="o":
			options.force=True
		elif outopt=="n":
			options.output=raw_input('Enter a new output file prefix: ')
			
	
#	if not "-Q" in sys.argv and options.program in ["BWA", "bwa"]:
#		options.mapq=20
	
	if options.ref=='':
		DoError('No reference dna file (-r) selected!')
	elif args==[] and options.studyid==0:
		DoError('No input files selected!')
	elif options.maxinsertsize>10000 or options.maxinsertsize<0:
		DoError('Maximum insert size (-i) must be between 0 and 10,000')
	elif options.maxinsertsize<options.mininsertsize:
		DoError('Minimum insert size (-j) must be smaller than maximum insert size (-i). Currently -i='+str(options.maxinsertsize)+' and -j='+str(options.mininsertsize))
	elif options.mininsertsize>10000 or options.mininsertsize<0:
		DoError('Minimum insert size (-j) must be between 0 and 10,000')
	elif options.quality>99 or options.quality<1:
		DoError('Ssaha mapping quality score (-q) must be between 1 and 100')
	elif options.program in ["BWA", "bwa"] and (options.mapq>30 or options.mapq<0):
		DoError('Mapping quality score (-Q) must be between 0 and 30 for bwa')
	elif options.program not in ["BWA", "bwa"] and (options.mapq>60 or options.mapq<0):
		DoError('Mapping quality score (-Q) must be between 0 and 60 for '+options.program)
#	elif options.snpquality>100 or options.snpquality<1:
#		DoError('Minimum site mapping quality for SNP calling (-q) must be between 1 and 100!')
	elif options.ratio>1 or options.ratio<0:
		DoError('SNP/site mapping quality ratio cutoff (-R) must be between 0 and 1')
#	elif options.readlength>1000 or options.readlength<36:
#		DoError('Read length (-l) must be between 36 and 1000!')
	elif options.mem>30 or options.mem<0:
		DoError('Memory requirement (-M) must be between 0 and 30Gb')
	
	options.program=options.program.upper()
	
	return




#####################
# SNPanalysis class #
#####################

class SNPanalysis:
	def __init__(self, fastq='', name='', mapped={}, runssaha='n', CDSseq='', number=0, log=[], keep=False):
		self.fastq=fastq
		self.name=name
		self.dirname=''
		self.fastqdir=''
		self.number=number
		self.log=[]
		self.keep=False
	
	def runSsaha(self, ref, bashfile):	
		print "\nRunning Ssaha on "+self.name+'...',
		sys.stdout.flush()
		#meaninsert=((options.maxinsertsize-options.mininsertsize)/2)+options.mininsertsize
		#Ssaha commands.
		
		#single end mapping
		if options.pairedend==False:
			self.keep=print_command_to_bashfile(' '.join(map(str, [SSAHA_DIR+"ssaha2 -solexa -output sam_soft -outfile "+self.dirname+"tmp1.sam "+ref+" "+self.fastqdir+self.name+".fastq"])), self.dirname, self.name, bashfile, self.log, self.keep)
						
			
		#paired end
		else:
			#self.keep=print_command_to_bashfile(' '.join(map(str, [SSAHA_DIR+"ssaha2 -solexa -score "+str(options.quality)+" -skip 2 -diff 0 -kmer 13 -outfile "+self.dirname+"tmp.sam -multi 12345 -pair "+str(options.mininsertsize)+","+str(options.maxinsertsize)+" -output sam_soft "+ref+" "+self.fastqdir+self.name+"_1.fastq "+self.fastqdir+self.name+"_2.fastq"
			self.keep=print_command_to_bashfile(' '.join(map(str, [SSAHA_DIR+"ssaha2 -solexa -outfile "+self.dirname+"tmp1.sam -pair "+str(options.mininsertsize)+","+str(options.maxinsertsize)+" -output sam_soft "+ref+" "+self.fastqdir+self.name+"_1.fastq "+self.fastqdir+self.name+"_2.fastq"])), self.dirname, self.name, bashfile, self.log, self.keep)
		
		#self.keep=print_command_to_bashfile(' '.join(map(str, [SAMTOOLS_DIR+"samtools view -b -q "+str(options.mapq)+" -S",self.dirname+"tmp1.sam -t "+ref+".fai >", self.dirname+"tmp.bam"
		
		self.keep=print_command_to_bashfile(' '.join(map(str, [SAMTOOLS_DIR+"samtools view -b -S",self.dirname+"tmp1.sam -t "+ref+".fai >", self.dirname+"tmp1.bam"])), self.dirname, self.name, bashfile, self.log, self.keep)
		
		
		if options.pairedend:
			self.keep=print_command_to_bashfile(' '.join(map(str, [MY_SCRIPTS_DIR+"fix_circular_bams.py -b", self.dirname+"tmp1.bam -o", self.dirname+"tmp"])), self.dirname, self.name, bashfile, self.log, self.keep)
			self.keep=print_command_to_bashfile(' '.join(map(str, ["rm", self.dirname+"tmp1.bam"])), self.dirname, self.name, bashfile, self.log, self.keep)
		else:
			self.keep=print_command_to_bashfile(' '.join(map(str, ["mv", self.dirname+"tmp1.bam", self.dirname+"tmp.bam"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
		self.keep=print_command_to_bashfile(' '.join(map(str, [SAMTOOLS_DIR+"samtools view -H ",self.dirname+"tmp.bam >", self.dirname+"tmp2.sam"])), self.dirname, self.name, bashfile, self.log, self.keep)
#		
		self.keep=print_command_to_bashfile(' '.join(map(str, ["cat", self.dirname+"tmp2.sam", self.dirname+"tmp1.sam > ", self.dirname+"tmp.sam"])), self.dirname, self.name, bashfile, self.log, self.keep)
#		
		self.keep=print_command_to_bashfile(' '.join(map(str, ["rm", self.dirname+"tmp2.sam", self.dirname+"tmp1.sam"])), self.dirname, self.name, bashfile, self.log, self.keep)
		#self.keep=print_command_to_bashfile(' '.join(map(str, [SAMTOOLS_DIR+"samtools reheader", self.dirname+"tmp1.sam", self.dirname+"tmp.bam"
	
	
	def runBWA(self, ref, bashfile):	
		print "\nRunning BWA on "+self.name+'...',
		sys.stdout.flush()
		
		#Map the reads against the genome
		if options.pairedend:
			self.keep=print_command_to_bashfile(' '.join(map(str, [BWA_DIR+"bwa aln -q 15 ", options.ref, self.fastqdir+self.name+"_1.fastq >", self.dirname+"tmp.F.sai"])), self.dirname, self.name, bashfile, self.log, self.keep)
			self.keep=print_command_to_bashfile(' '.join(map(str, [BWA_DIR+"bwa aln -q 15 ", options.ref, self.fastqdir+self.name+"_2.fastq >", self.dirname+"tmp.R.sai"])), self.dirname, self.name, bashfile, self.log, self.keep)
			#Join both aligments
			self.keep=print_command_to_bashfile(' '.join(map(str, [BWA_DIR+"bwa sampe", ref,  self.dirname+"tmp.F.sai", self.dirname+"tmp.R.sai", self.fastqdir+self.name+"_1.fastq", self.fastqdir+self.name+"_2.fastq >", self.dirname+"tmp.sam"])), self.dirname, self.name, bashfile, self.log, self.keep)
		else:
			self.keep=print_command_to_bashfile(' '.join(map(str, [BWA_DIR+"bwa aln -q 15 ", options.ref, self.fastqdir+self.name+".fastq >", self.dirname+"tmp.F.sai"])), self.dirname, self.name, bashfile, self.log, self.keep)
			self.keep=print_command_to_bashfile(' '.join(map(str, [BWA_DIR+"bwa samse", ref,  self.dirname+"tmp.F.sai", self.fastqdir+self.name+".fastq >", self.dirname+"tmp.sam"])), self.dirname, self.name, bashfile, self.log, self.keep)
		
		
		#produce the BAM file
		#self.keep=print_command_to_bashfile(' '.join(map(str, [SAMTOOLS_DIR+"samtools view -b -q "+str(options.mapq)+" -S", self.dirname+"tmp.sam >", self.dirname+"tmp.bam"
		
		self.keep=print_command_to_bashfile(' '.join(map(str, [SAMTOOLS_DIR+"samtools view -b -S",self.dirname+"tmp1.sam -t "+ref+".fai >", self.dirname+"tmp1.bam"])), self.dirname, self.name, bashfile, self.log, self.keep)
		
		if options.pairedend:
			self.keep=print_command_to_bashfile(' '.join(map(str, [MY_SCRIPTS_DIR+"fix_circular_bams.py -b", self.dirname+"tmp1.bam -o", self.dirname+"tmp"])), self.dirname, self.name, bashfile, self.log, self.keep)
			self.keep=print_command_to_bashfile(' '.join(map(str, ["rm", self.dirname+"tmp1.bam"])), self.dirname, self.name, bashfile, self.log, self.keep)
		else:
			self.keep=print_command_to_bashfile(' '.join(map(str, ["mv", self.dirname+"tmp1.bam", self.dirname+"tmp.bam"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			
			
			
	def runSMALT(self, ref, bashfile):	
		print "\nRunning SMALT on "+self.name+'...',
		sys.stdout.flush()

		#Map the reads against the genome
		if options.pairedend:
			self.keep=print_command_to_bashfile(' '.join(map(str, [SMALT_DIR+"smalt map -x -i", options.maxinsertsize, " -j", options.mininsertsize, " -f samsoft -o "+self.dirname+"tmp1.sam", tmpname+".index", self.fastqdir+self.name+"_1.fastq", self.fastqdir+self.name+"_2.fastq"])), self.dirname, self.name, bashfile, self.log, self.keep)
		else:
			self.keep=print_command_to_bashfile(' '.join(map(str, [SMALT_DIR+"smalt map -x -f samsoft -o "+self.dirname+"tmp1.sam", tmpname+".index", self.fastqdir+self.name+".fastq"])), self.dirname, self.name, bashfile, self.log, self.keep)
		
		
		#produce the BAM file
		#self.keep=print_command_to_bashfile(' '.join(map(str, [SAMTOOLS_DIR+"samtools view -b -q "+str(options.mapq)+" -S", self.dirname+"tmp1.sam -t "+ref+".fai >", self.dirname+"tmp.bam"
		
		self.keep=print_command_to_bashfile(' '.join(map(str, [SAMTOOLS_DIR+"samtools view -b -S",self.dirname+"tmp1.sam -t "+ref+".fai >", self.dirname+"tmp1.bam"])), self.dirname, self.name, bashfile, self.log, self.keep)
		
		if options.pairedend:
			self.keep=print_command_to_bashfile(' '.join(map(str, [MY_SCRIPTS_DIR+"fix_circular_bams.py -b", self.dirname+"tmp1.bam -o", self.dirname+"tmp"])), self.dirname, self.name, bashfile, self.log, self.keep)
			self.keep=print_command_to_bashfile(' '.join(map(str, ["rm", self.dirname+"tmp1.bam"])), self.dirname, self.name, bashfile, self.log, self.keep)
		else:
			self.keep=print_command_to_bashfile(' '.join(map(str, ["mv", self.dirname+"tmp1.bam", self.dirname+"tmp.bam"])), self.dirname, self.name, bashfile, self.log, self.keep)
		
		#add header to sam file for making plots - no need. Can read bams too!
		self.keep=print_command_to_bashfile(' '.join(map(str, [SAMTOOLS_DIR+"samtools view -H ",self.dirname+"tmp.bam >", self.dirname+"tmp2.sam"])), self.dirname, self.name, bashfile, self.log, self.keep)
#		
		self.keep=print_command_to_bashfile(' '.join(map(str, ["cat", self.dirname+"tmp2.sam", self.dirname+"tmp1.sam > ", self.dirname+"tmp.sam"])), self.dirname, self.name, bashfile, self.log, self.keep)
#		
		self.keep=print_command_to_bashfile(' '.join(map(str, ["rm", self.dirname+"tmp2.sam", self.dirname+"tmp1.sam"])), self.dirname, self.name, bashfile, self.log, self.keep)
		#self.keep=print_command_to_bashfile(' '.join(map(str, [SAMTOOLS_DIR+"samtools reheader", self.dirname+"tmp1.sam", self.dirname+"tmp.bam"

		
	def makepileup_from_sam(self, ref, bashfile):
		
		if self.runmapping:
		
			#create plots from sam file
			if options.plots==True:
				self.keep=print_command_to_bashfile(' '.join(map(str, [MY_SCRIPTS_DIR+"odd_plots_from_sam.py", self.dirname+"tmp.sam", self.dirname+self.name+" >  /dev/null 2>&1"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			#remove the sam file as it is no longer needed
			self.keep=print_command_to_bashfile(' '.join(map(str, ["rm", self.dirname+"tmp.sam"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			#order the bam file
			self.keep=print_command_to_bashfile(' '.join(map(str, [SAMTOOLS_DIR+"samtools sort ", self.dirname+"tmp.bam", self.dirname+self.name])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			#index the bam file, to get the bai file.
			self.keep=print_command_to_bashfile(' '.join(map(str, [SAMTOOLS_DIR+"samtools index",  self.dirname+self.name+".bam"])), self.dirname, self.name, bashfile, self.log, self.keep)
		
		if self.findindels:
		
			self.keep=print_command_to_bashfile(' '.join(map(str, [BREAKDANCER_DIR+"perl/bam2cfg.pl",  self.dirname+self.name+".bam >", self.dirname+"tmp.config"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			self.keep=print_command_to_bashfile(' '.join(map(str, [BREAKDANCER_DIR+"cpp/breakdancer_max",  self.dirname+"tmp.config >", self.dirname+self.name+"_breakdancer.out"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			self.keep=print_command_to_bashfile(' '.join(map(str, [MY_SCRIPTS_DIR+"breakdancer_2_tab.py",  self.dirname+self.name+"_breakdancer.out", self.dirname+self.name+"_breakdancer.tab"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			self.keep=print_command_to_bashfile(' '.join(map(str, [SAMTOOLS_DIR+"samtools view ",  self.dirname+self.name+".bam |", PINDEL_DIR+"sam_2_pindel -", self.dirname+"tmp.pindelin 300 "+self.name+" 0"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			self.keep=print_command_to_bashfile(' '.join(map(str, [PINDEL_DIR+"pindel -x 6 -f "+ref+" -p", self.dirname+"tmp.pindelin -c ALL -o", self.dirname+"tmp.pindel.out", "-b", self.dirname+self.name+"_breakdancer.out"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			self.keep=print_command_to_bashfile(' '.join(map(str, ["cat", self.dirname+"tmp.pindel.out_[ID] >", self.dirname+"tmp.pindel.out"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			
			#NOTE: I edited this script so insertions in deletions are not output
			self.keep=print_command_to_bashfile(' '.join(map(str, [PINDEL_DIR+"vcfcreator -p", self.dirname+"tmp.pindel.out -r", ref, "-R", ref.split("/")[-1].split(".")[0], "-d", ''.join(map(str,[1979, 1, 30]))])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			
			self.keep=print_command_to_bashfile(' '.join(map(str, ['grep -A 1 "######" '+self.dirname+'tmp.pindel.out_* |grep -v "#####" |grep -v "\-\-" >', self.dirname+self.name+"_pindel.out"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			self.keep=print_command_to_bashfile(' '.join(map(str, [MY_SCRIPTS_DIR+"pindel_2_tab.py", self.dirname+self.name+"_pindel.out", self.dirname+self.name+"_pindel.tab"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			
			#NOTE: I edited this script to only output indels<10kb
			self.keep=print_command_to_bashfile(' '.join(map(str, [DINDEL_SCRIPTS_DIR+"convertVCFToDindel.py -i", self.dirname+"tmp.pindel.out.vcf", "-r", ref, "-o", self.dirname+"tmp.pindel.out.dindel.in"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			
			self.keep=print_command_to_bashfile(' '.join(map(str, [DINDEL_DIR+"dindel-1.01-linux-64bit --analysis realignCandidates --varFile ", self.dirname+"tmp.pindel.out.dindel.in", "--outputFile", self.dirname+"tmp.pindel.out.dindel.realigned.txt", "--ref", ref])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			
			self.keep=print_command_to_bashfile(' '.join(map(str, [DINDEL_DIR+"dindel-1.01-linux-64bit --analysis getCIGARindels --bamFile", self.dirname+self.name+".bam", "--ref", ref, "--outputFile", self.dirname+"tmp.dindel.out"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			
			self.keep=print_command_to_bashfile(' '.join(map(str, [DINDEL_SCRIPTS_DIR+"selectCandidates.py -i", self.dirname+"tmp.dindel.out.variants.txt -o", self.dirname+"tmp.dindel.out.variants_min2.txt"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			
			self.keep=print_command_to_bashfile(' '.join(map(str, ["cat", self.dirname+"tmp.dindel.out.variants_min2.txt", self.dirname+"tmp.pindel.out.dindel.realigned.txt.variants.txt >", self.dirname+"tmp.dindel_all.txt"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			
			self.keep=print_command_to_bashfile(' '.join(map(str, [DINDEL_SCRIPTS_DIR+"makeWindows.py --inputVarFile", self.dirname+"tmp.dindel_all.txt", "--windowFilePrefix", self.dirname+"tmp.dindel.windows", "--numWindowsPerFile 100000"])), self.dirname, self.name, bashfile, self.log, self.keep)
			

			self.keep=print_command_to_bashfile(' '.join(map(str, [DINDEL_DIR+"dindel-1.01-linux-64bit --analysis indels --doDiploid --bamFile", self.dirname+self.name+".bam", "--ref", ref, "--outputFile", self.dirname+"tmp.dindel.output", "--outputRealignedBAM --processRealignedBAM", MY_SCRIPTS_DIR+"dindeltobam.sh --varFile", self.dirname+"tmp.dindel.windows.1.txt", "--libFile", self.dirname+"tmp.dindel.out.libraries.txt", self.dirname+"tmp.dindel.windows.1.txt_stage2"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			
			self.keep=print_command_to_bashfile(' '.join(map(str, ["echo", self.dirname+"tmp.dindel.output.glf.txt >", self.dirname+"tmp.dindelfiles"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			self.keep=print_command_to_bashfile(' '.join(map(str, [DINDEL_SCRIPTS_DIR+"mergeOutputDiploid.py --inputFiles", self.dirname+"tmp.dindelfiles --outputFile", self.dirname+self.name+"_dindel.vcf --ref "+ref])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			self.keep=print_command_to_bashfile(' '.join(map(str, [MY_SCRIPTS_DIR+"dindel_2_tab.py ", self.dirname+self.name+"_dindel.vcf", self.dirname+self.name+"_dindel.tab"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			self.keep=print_command_to_bashfile(' '.join(map(str, [SAMTOOLS_DIR+"samtools view -b -S -T", ref, "-o", self.dirname+"tmp.dindel.bam", self.dirname+"tmp.dindel.output_realigned.sam"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			self.keep=print_command_to_bashfile(' '.join(map(str, [SAMTOOLS_DIR+"samtools sort", self.dirname+"tmp.dindel.bam", self.dirname+self.name+"_dindel"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			self.keep=print_command_to_bashfile(' '.join(map(str, [SAMTOOLS_DIR+"samtools index", self.dirname+self.name+"_dindel.bam"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			self.keep=print_command_to_bashfile(' '.join(map(str, [MY_SCRIPTS_DIR+"realign_dindel_bams.py", self.dirname+self.name+"_dindel.bam", self.dirname+self.name+"_dindel.vcf", self.dirname+self.name+".bam", self.dirname+"tmp.realigned.bam"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			self.keep=print_command_to_bashfile(' '.join(map(str, [SAMTOOLS_DIR+"samtools sort", self.dirname+"tmp.realigned.bam", self.dirname+self.name+"_realigned"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			self.keep=print_command_to_bashfile(' '.join(map(str, [SAMTOOLS_DIR+"samtools index", self.dirname+self.name+"_realigned.bam"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			self.keep=print_command_to_bashfile(' '.join(map(str, ["rm -f", self.dirname+"tmp.realigned.bam"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			bamfilename=self.dirname+self.name+"_realigned.bam"
		else:
			bamfilename=self.dirname+self.name+".bam"
		
			
			
		
		if self.makebcf:
			
			
			if options.properpairfilter:
				bamname=self.dirname+self.name+"_filtered.bam"
				self.keep=print_command_to_bashfile(' '.join(map(str, [SAMTOOLS_DIR+"samtools view -b -h -f 2 -o ", bamname, bamfilename])), self.dirname, self.name, bashfile, self.log, self.keep)
				
				self.keep=print_command_to_bashfile(' '.join(map(str, ["rm -f", bamfilename, bamfilename+".bai"])), self.dirname, self.name, bashfile, self.log, self.keep)
				
				self.keep=print_command_to_bashfile(' '.join(map(str, [SAMTOOLS_DIR+"samtools index", bamname])), self.dirname, self.name, bashfile, self.log, self.keep)
				
			else:
				bamname=bamfilename
		
		
			#produce the pileup file
			
			if options.BAQ:
				self.keep=print_command_to_bashfile(' '.join(map(str, [SAMTOOLS_DIR+"samtools mpileup -d 1000 -DSugBf ", ref, bamname, ">", self.dirname+"tmp.mpileup"])), self.dirname, self.name, bashfile, self.log, self.keep)
				
			else:
				self.keep=print_command_to_bashfile(' '.join(map(str, [SAMTOOLS_DIR+"samtools mpileup -d 1000 -DSugf ", ref, bamname, ">", self.dirname+"tmp.mpileup"])), self.dirname, self.name, bashfile, self.log, self.keep)
				
				
			
			self.keep=print_command_to_bashfile(' '.join(map(str, [BCFTOOLS_DIR+"bcftools view -bcg", self.dirname+"tmp.mpileup >", self.dirname+self.name+".bcf"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			
			self.keep=print_command_to_bashfile(' '.join(map(str, [BCFTOOLS_DIR+"bcftools index", self.dirname+self.name+".bcf"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			
			self.keep=print_command_to_bashfile(' '.join(map(str, [BCFTOOLS_DIR+"bcftools view -bcgv", self.dirname+"tmp.mpileup >", self.dirname+self.name+"_variant.bcf"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			
			self.keep=print_command_to_bashfile(' '.join(map(str, [BCFTOOLS_DIR+"bcftools index", self.dirname+self.name+"_variant.bcf"])), self.dirname, self.name, bashfile, self.log, self.keep)
			
				
		# clean up:
		if self.makebcf and not options.dirty:
			self.keep=print_command_to_bashfile(' '.join(map(str, ["rm", self.dirname+"tmp.*"])), self.dirname, self.name, bashfile, self.log, self.keep)
			

		
		
		if self.makebcf or self.makemfa:
		
			if options.pseudosequence:
			
#				self.keep=print_command_to_bashfile(' '.join(map(str, [MY_SCRIPTS_DIR+"add_dindel_indels_to_bcf.py -b ", self.dirname+self.name+".bcf", "-v ", self.dirname+self.name+"_dindel.vcf", "-r ", options.ratio, "-d ", options.depth, "-D ", options.stranddepth, "-q ", options.quality, "-m ", options.mapq, "-o", self.dirname])), self.dirname, self.name, bashfile, self.log, self.keep)
			
			
				self.keep=print_command_to_bashfile(' '.join(map(str, [MY_SCRIPTS_DIR+"bcf_2_pseudosequence_testing.py -b ", self.dirname+self.name+".bcf", "-B ", bamname, "-r ", options.ratio, "-d ", options.depth, "-D ", options.stranddepth, "-q ", options.quality, "-m ", options.mapq, "-o", self.dirname+self.name])), self.dirname, self.name, bashfile, self.log, self.keep)
				
				
		if options.plots:
			self.keep=print_command_to_bashfile(' '.join(map(str, ['gzip -f '+self.dirname+self.name+'*.plot'])), self.dirname, self.name, bashfile, self.log, self.keep)
			
		self.keep=print_command_to_bashfile(' '.join(map(str, [MY_SCRIPTS_DIR+'heterozygosity_plot.py -b', self.dirname+self.name+".bcf -o", self.dirname+self.name+"_contamination_plot.pdf", "-r", options.ratio, "-d", options.depth, "-D", options.stranddepth, "-q", options.quality])), self.dirname, self.name, bashfile, self.log, self.keep)
		
	
	
		
########
# Main #
########


if __name__ == "__main__":



	(options, args)=get_user_options()
	#print options, args
	check_input_validity(options, args)


	
	
	#make random name for files
	chars = string.ascii_letters + string.digits
	tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))
	
	pools=[]
	count=0
	poolsort=[]
	ziplist={}
	
	bamlist={}
	
	
	inputfiles=[]
	
	if options.studyid!=0:
		print "\nFinding fastq files in study "+str(options.studyid)+"..."
		sys.stdout.flush()
		pathfindstring="pathfind -t study -f fastq -id "+str(options.studyid)
		if options.passed:
			pathfindstring=pathfindstring+" -qc passed"
		pathfindarg = shlex.split(pathfindstring)
		returnval = subprocess.Popen(pathfindarg, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		inputfiles = returnval.stdout.read().replace(" ","").split("\n")
		for infile in inputfiles:
			if not os.path.isfile(infile):
				inputfiles.remove(infile)

		
	inputfiles=inputfiles+args
	
	if len(inputfiles)==0:
		DoError("No input files found. If you provided a studyid, please check it is correct")
	
	print '\nChecking input files...'
	sys.stdout.flush()
	
	
	for pool in inputfiles:
		
		if not os.path.isfile(pool):
			print "File "+pool+" not found! Skipping..."
			continue
		
		filetype='.'+pool.split('.')[-1]
		if filetype not in ['.fastq','.bam',".gz"]:
			print "WARNING: Input file name is not .fastq or .bam!"
		
		originalfastqdir=''
		if pool[-1]=='/':
			pool=pool[:-1]
		if len(pool.split('/'))>1:
			originalfastqdir='/'.join(pool.split('/')[:-1])+'/'
		
		if pool.split('.')[-1]=="gz" and pool.split('.')[-2]=="fastq":
			if not os.path.isdir(tmpname+"_unzipped"):
				os.system("mkdir "+tmpname+"_unzipped")
			#print "Unzipping ", pool
			
			#os.system("zcat "+pool+" > "+tmpname+"_unzipped/"+'.'.join(pool.split('/')[-1].split('.')[:-1]))
			#pool='.'.join(pool.split('.')[:-1])

			if options.pairedend==True:
				if not '.'.join(pool.replace("#","_").split('/')[-1].split('.')[:-2])[:-2] in ziplist:
					ziplist['.'.join(pool.replace("#","_").split('/')[-1].split('.')[:-2])[:-2]]=[pool]
				else:
					ziplist['.'.join(pool.replace("#","_").split('/')[-1].split('.')[:-2])[:-2]].append(pool)
			else:
				ziplist['.'.join(pool.replace("#","_").split('/')[-1].split('.')[:-2])]=[pool]
			originalpool='.'.join(pool.split('/')[-1].split('.')[:-2])
			pool=tmpname+"_unzipped/"+'.'.join(pool.replace("#","_").split('/')[-1].split('.')[:-1])
			#ziplist.append(pool)
			
		
		elif pool.split('.')[-1]=="bam":
			if not os.path.isdir(tmpname+"_unbammed"):
				os.system("mkdir "+tmpname+"_unbammed")
			bamlist['.'.join(pool.replace("#","_").split('/')[-1].split('.')[:-1])]=pool
			originalpool='.'.join(pool.split('/')[-1].split('.')[:-1])
			pool=tmpname+"_unbammed/"+'.'.join(pool.replace("#","_").split('/')[-1].split('.')[:-1])+".bam"
		
		else:
			originalpool='.'.join(pool.split('/')[-1].split('.')[:-1])
			
	
#		if not os.path.isfile(pool):
#			print "File "+pool+" not found! Skipping..."
#			continue
		fastqdir=''
		if pool[-1]=='/':
			pool=pool[:-1]
		if len(pool.split('/'))>1:
			fastqdir='/'.join(pool.split('/')[:-1])+'/'
		
		
		pool='.'.join(pool.split('/')[-1].split('.')[:-1])
		pool=pool.replace("_nonhuman","")
		
		if options.pairedend==True and filetype!=".bam":
			if pool[-2:]=='_1':
				pool=pool[:-2]
				originalpool=originalpool[:-2]
				if not os.path.isfile(originalfastqdir+originalpool+"_2.fastq") and not os.path.isfile(originalfastqdir+originalpool+"_2.fastq.gz"):
					print "File "+originalpool+"_2.fastq not found! Skipping..."
					continue
			elif pool[-2:]=='_2':
				pool=pool[:-2]
				originalpool=originalpool[:-2]
				if not os.path.isfile(originalfastqdir+originalpool+"_1.fastq") and not os.path.isfile(originalfastqdir+originalpool+"_1.fastq.gz"):
					print "File "+originalpool+"_1.fastq not found! Skipping..."
					continue
		
		name=pool
		
		pool=pool+"_"+options.program

		if pool in poolsort:
			continue
		print pool+'...',
		sys.stdout.flush()

		if not os.path.isdir(pool):
		 	print "pool "+pool+" not found! Creating...",
			os.system("mkdir "+pool)
						
		pools.append(SNPanalysis())
		pools[count].number=str(count+1)
		pools[count].dirname=pool+"/"
		pools[count].name=name
		pools[count].fastqdir=fastqdir+"/"
		pools[count].filetype=filetype
#		the rest can be globals
#		pools[count].quality=options.qualityoptions
#		pools[count].mininsertsize=options.mininsertsize
#		pools[count].maxinsertsize=options.maxinsertsize
#		pools[count].pairedend=options.pairedend
#		pools[count].readlength=options.readlength
		
		if options.keep and os.path.isfile(pools[count].dirname+pools[count].name+".bam"):
			
#			try:
#				samfile = pysam.Samfile( pools[count].dirname+pools[count].name+".bam", "rb" )
#				bcount=0
#				print samfile.references
#				print dir(samfile)
#				for read in samfile:
#					print read
#					bcount+=1
#				samfile.close
#				print bcount
#				pools[count].runmapping=False
#			except StandardError:
#				print "Not a bam file"
#				pools[count].runmapping=True
#			iterator = samfile
#			for read in iterator:
#				print read
#			bcount=0
#			for read in samfile:
#				print read
#				bcount+=1
#			print bcount
			pools[count].runmapping=False
		else:
			pools[count].runmapping=True
		
		
		
		pools[count].makebcf=True
		
		if options.pseudosequence:
			pools[count].makemfa=True
		else:
			pools[count].makemfa=False
		
		if options.indels:
			pools[count].findindels=True
		else:
			pools[count].findindels=False
		
		if options.keep:
			pools[count].keep=True
		else:
			pools[count].keep=False
		
		
		
		
		poolsort.append(pool)
		print 'ok'
		sys.stdout.flush()
		
		count=count+1
		
		
	if len(pools)==0:
		print "\nError: No valid input files!"
		sys.exit()

	
	
	seq_records=[]
	
	for seq_record in SeqIO.parse(open(options.ref), "fasta"):
		seq_records.append(seq_record)
		
	
	
	if len(seq_records)==0:
		DoError("Cannot open reference fasta file!")
	else:
	
		SeqIO.write(seq_records, open(options.ref,"w"), "fasta")
		if options.incref:
			concatenated_seq=""
			
			contigs={}
			contigorder=[]
			for record in seq_records:
				contigorder.append(record.id)
				contigs[record.id]=str(record.seq)
			
			keys=contigs.keys()
			keys.sort()
			for contig in contigorder:
				concatenated_seq=concatenated_seq+contigs[contig]
				
			my_seq_record = SeqRecord(Seq(concatenated_seq))
			my_seq_record.id=options.ref.split("/")[-1].split(".")[0]
			my_seq_record.description="Reference"
#			if options.pseudosequence:
#				SeqIO.write([my_seq_record], open(options.output+".aln","w"), "fasta")
		
		else:
			os.system("rm "+options.output+".aln")
	
	
	#Running Ssaha or bwa where required
	if options.program=='BWA':
		os.system("bwa index "+options.ref)
	elif options.program=='SMALT':
		os.system(SMALT_DIR+"smalt index -k 13 -s 1 "+tmpname+".index "+options.ref)
	#elif options.program=='ssaha':
	#	os.system("samtools faidx "+options.ref)
	os.system("samtools faidx "+options.ref)
		
	count=0
	
	for pool in pools:
		
		
		
		if options.keep and os.path.isfile(pool.dirname+pool.name+".log"):
			pool.log=open(pool.dirname+pool.name+".log", "rU").readlines()
		if os.path.isfile(pool.dirname+pool.name+".log"):
			os.remove(pool.dirname+pool.name+".log")
	
		count+=1
		bashfile=open(str(count)+tmpname+'_sbs.sh','w')
		
		if pool.runmapping and pool.name in ziplist:
			for filename in ziplist[pool.name]:
				pool.keep=print_command_to_bashfile(' '.join(map(str, ["zcat "+filename+" > "+pool.fastqdir+'.'.join(filename.replace("#","_").split('/')[-1].split('.')[:-1])])), pool.dirname, pool.name, bashfile, pool.log, pool.keep)
		
		
		elif pool.runmapping and pool.name in bamlist:
			if options.pairedend:
				pool.keep=print_command_to_bashfile(' '.join(map(str, [MY_SCRIPTS_DIR+'bam_filter.py -t all -b '+bamlist[pool.name]+' -o '+pool.fastqdir+pool.name])), pool.dirname, pool.name, bashfile, pool.log, pool.keep)
			else:
				pool.keep=print_command_to_bashfile(' '.join(map(str, [MY_SCRIPTS_DIR+'bam_filter.py -t all -f fastq -b '+bamlist[pool.name]+' -o '+pool.fastqdir+pool.name])), pool.dirname, pool.name, bashfile, pool.log, pool.keep)
		
		if pool.runmapping:
			if options.program=='BWA':
				pool.runBWA(options.ref, bashfile)
			elif options.program=='SMALT':
				pool.runSMALT(options.ref, bashfile)	
			elif options.program=='SSAHA':
				pool.runSsaha(options.ref, bashfile)
			
		 
		pool.makepileup_from_sam(options.ref, bashfile)
	
	
		bashfile.close()


		if options.LSF==False:
			os.system('bash '+str(count)+tmpname+'_sbs.sh')
			if not options.dirty:
				os.system('rm '+str(count)+tmpname+'_sbs.sh '+options.ref+'.pac '+options.ref+'.ann '+options.ref+'.rpac '+options.ref+'.amb '+options.ref+'.rbwt '+options.ref+'.bwt '+options.ref+'.sa '+options.ref+'.rsa')
	
	if options.pseudosequence:
		argstring=[]
		for pool in pools:	
			argstring.append(pool.dirname+pool.name+".mfa")
		if options.indels:
			joinstring=MY_SCRIPTS_DIR+"join_dna_files_with_indels_new.py -r "+options.ref+" -o "+options.output+".aln "+' '.join(argstring) #*_ssaha/*_test.mfa
		else:
			if not options.incref:
				joinstring="'cat "+' '.join(argstring)+" > "+options.output+".aln'"
			else:
				joinstring="'cat "+options.ref+" "+' '.join(argstring)+" > "+options.output+".aln'"
				
		summarystring=MY_SCRIPTS_DIR+'summarise_snps.py -g -w -r '+options.ref.split("/")[-1].split(".")[0]+' -o '+options.output+' -i '+options.output+'.aln'
	
		#print summarystring, options.embl
		
		if options.embl!="":
			summarystring=summarystring+" -e "+options.embl
		if options.alnfile:
			summarystring=summarystring+" -a"
		if options.tabfile:
			summarystring=summarystring+" -t"
		if options.raxml:
			summarystring=summarystring+" -p -l -b "+str(options.bootstrap)
		
	print
	if options.LSF==True:
		if count>0:
			if options.mem>0:
				memkb=str(options.mem*1000000)
				memmb=str(options.mem*1000)
				os.system('echo \'bash ${LSB_JOBINDEX}'+tmpname+'_sbs.sh\' | bsub -R \'select[mem>'+memmb+'] rusage[mem='+memmb+']\' -q '+options.LSFQ+'  -J'+tmpname+'_'+options.program+'"[1-'+str(count)+']%'+str(options.nodes)+'"  -M '+memkb+' -o '+tmpname+options.program+'-%I.out -e '+tmpname+options.program+'-%I.err')# run all ssaha jobs in job array
			else:
				os.system('echo \'bash ${LSB_JOBINDEX}'+tmpname+'_sbs.sh\' | bsub -q '+options.LSFQ+' -J'+tmpname+'_'+options.program+'"[1-'+str(count)+']%'+str(options.nodes)+'" -o '+tmpname+options.program+'-%I.out -e '+tmpname+options.program+'-%I.err')
			
			if not options.dirty:
					os.system('bsub -w \'ended('+tmpname+'_'+options.program+')\' rm -rf *'+tmpname+'_sbs.sh '+tmpname+'*.out '+tmpname+'*.err '+options.ref+'.pac '+options.ref+'.ann '+options.ref+'.rpac '+options.ref+'.amb '+options.ref+'.rbwt '+options.ref+'.bwt '+options.ref+'.sa '+options.ref+'.rsa '+tmpname+'*.sma '+tmpname+'*.smi '+tmpname+'_unzipped '+tmpname+'_unbammed')#when job array is all done delete
		
		
		if options.pseudosequence:
			if len(pools)>200:
				if count==0:
					os.system('echo '+joinstring+' | bsub -M 10000000 -q long -R \'select[mem>10000] rusage[mem=10000]\' -J'+tmpname+'_joining -o '+options.output+'_join.bsubout -e '+options.output+'_join.bsuberr')
				else:
					os.system('echo '+joinstring+' | bsub -M 10000000 -q long -R \'select[mem>10000] rusage[mem=10000]\' -J'+tmpname+'_joining -w \'ended('+tmpname+'_'+options.program+')\' -o '+options.output+'_join.bsubout -e '+options.output+'_join.bsuberr ')
				os.system('echo '+summarystring+' | bsub -M 16000000 -q long -R \'select[mem>16000] rusage[mem=16000]\' -w \'ended('+tmpname+'_joining)\' -o '+options.output+'_summary.bsubout -e '+options.output+'_summary.bsuberr')
			else:
				if count==0:
					os.system('echo '+joinstring+' | bsub -J'+tmpname+'_joining -o '+options.output+'_join.bsubout -e '+options.output+'_join.bsuberr')
				else:
					os.system('echo '+joinstring+' | bsub -J'+tmpname+'_joining -w \'ended('+tmpname+'_'+options.program+')\' -o '+options.output+'_join.bsubout -e '+options.output+'_join.bsuberr')
				os.system('echo '+summarystring+' | bsub -M 16000000 -q long -R \'select[mem>16000] rusage[mem=16000]\' -w \'ended('+tmpname+'_joining)\' -o '+options.output+'_sum.bsubout -e '+options.output+'_sum.bsuberr')
			
		#os.system('bsub -w \'ended('+tmpname+'_'+program+')\' rm *'+tmpname+'_sbs.sh; for f in tmpname*.err; do if test ! -s $f;then rm $f ${f%.err}.out;fi;done')
	
	elif options.pseudosequence:
		os.system(joinstring)
		os.system(summarystring)

		

			
			
			
			
