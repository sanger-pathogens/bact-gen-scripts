#!/usr/bin/env python


##################
# Import modules #
##################

import string, re
import os, sys, getopt, random, math, time, datetime
from random import *
from optparse import OptionParser, OptionGroup
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



####################
# Set some globals #
####################

SAMTOOLS_DIR=""
BCFTOOLS_DIR=""
SSAHA_DIR="/nfs/users/nfs_s/sh16/ssaha2_v2.5.1_x86_64/"
BWA_DIR=""

#SMALT_DIR="smalt"
MY_SCRIPTS_DIR="/nfs/users/nfs_s/sh16/scripts/"
GATK_LOC="/software/vertres/bin-external/GenomeAnalysisTK-1.5-9-ga05a7f2/GenomeAnalysisTK.jar"
JAVA_DIR="/software/jdk1.6.0_01/bin/"


##########################
# Error message function #
##########################

def DoError(errorstring):
	print "\nError:", errorstring
	print "\nFor help use -h or --help\n"
	sys.exit()


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
	group.add_option("-p", "--program", action="store", type="choice", dest="program", choices=["bwa","ssaha", "smalt", "BWA","SSAHA", "SMALT"], help="Mapping program to use (choose from bwa, ssaha or smalt) [default= %default]", default="smalt")
	group.add_option("-1", "--nomap", action="store_false", dest="domapping", help="Do not remap data - only available when input is bam (default is to map)", default=True)
	group.add_option("-v", "--smaltversion", action="store", type="choice", dest="version", choices=["latest","0.5.8", "0.6.3", "0.6.4", "0.7.4"], help="Version of SMALT to use (for backward compatibility). Choose from 0.5.8, 0.6.3, 0.6.4, 0.7.4 and latest (currently 0.7.4) [default= %default]", default="0.5.8")
	group.add_option("-H", "--human", action="store_true", dest="human", help="Mapping against human (or other large euk)", default=False)
	#group.add_option("-l", "--length", action="store", dest="readlength", help="Read length [default= %default]", default=54, type="int", metavar="INT")
	group.add_option("-s", "--single", action="store_false", dest="pairedend", help="reads are single ended (not paired)", default=True)
	group.add_option("-i", "--maxinsert", action="store", dest="maxinsertsize", help="maximum insert size (ssaha and smalt only) [default= %default]", default=1000, type="int", metavar="INT")
	group.add_option("-j", "--mininsert", action="store", dest="mininsertsize", help="minimum insert size (ssaha and smalt only) [default= %default]", default=50, type="int", metavar="INT")
	group.add_option("-S", "--ssahaquality", action="store", dest="ssahaquality", help="minimum ssaha quality score while mapping (ssaha only) [default= %default]", default=30, type="int", metavar="INT")
	group.add_option("-E", "--maprepeats", action="store_true", dest="maprepeats", help="randomly map repeats when using SMALT (default is to not map repeats)", default=False)
	group.add_option("-z", "--nomapid", action="store", dest="nomapid", help="Minimum identity threshold to report a mapping Specified as a positive integer or proportion of read length (smalt only) [default= %default]", default=0, type="float", metavar="float")
	group.add_option("-G", "--GATK", action="store_true", dest="GATK", help="Turn on GATK indel realignment (highly recommended). [Default= %default]", default=False)
	parser.add_option("-2", "--issequence", action="store", dest="ISfasta", help="Input fasta file of IS element file", default="")
	
	
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
	group.add_option("-c", "--circular", action="store_false", dest="circular", help="Contigs are not circular, so do not try to fix them", default=True)
	#parser.add_option("-q", "--quality", action="store", dest="quality", help="Minimum base quality [default= %default]", default=120, type="int")
	#group.add_option("-S", "--RMS", action="store", dest="RMS", help="Minimum root mean squared mapping quality [default= %default]", default=25, type="int")
	#parser.add_option("-Q", "--strandquality", action="store", dest="strandquality", help="Minimum per strand base quality [default= %default]", default=60, type="int")
	group.add_option("-R", "--ratio", action="store", dest="ratio", help="SNP/Mapping quality ratio cutoff [default= %default]", default=0.75, type="float")
	
	#group.add_option("-S", "--SNPquality", action="store", type="int", dest="snpquality", help="Minimum site mapping quality for SNP calling [default= %default]", default=90, metavar="INT")
	#group.add_option("-R", "--ratio", action="store", type="float", dest="ratio", help="SNP/site mapping quality ratio cutoff [default= %default]", default=0.75, metavar="FLOAT")
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "Output options")
	group.add_option("-e", "--embl", action="store", dest="embl", help="reference annotation (in true embl or genbank format) [required for dN/dS calculations and CDS naming]", default="", metavar="FILE")
	group.add_option("-o", "--output", action="store", dest="output", help="output file prefix", default="")
	group.add_option("-O", "--dir_output", action="store", dest="diroutput", help="output directory suffix", default="")
	group.add_option("-f", "--force", action="store_true", dest="force", help="force overwrite of output files", default=False)
	group.add_option("-F", "--filter_bam", action="store", dest="filter", help="filter or split bam file. Choose one of the following: 1) Include all reads, 2) Include only mapped reads 3) Include properly paired reads 4) Split into two files: mapped reads and unmapped reads, 5) Split into two files: properly paired and not properly paired [Default=%default]", type="choice", choices=["1", "2", "3", "4", "5"], default="1")
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
	group.add_option("-n", "--nodes", action="store", dest="nodes", help="Maximum number of jobs to run on nodes in parallel. [Default= %default]", default=20, type="int")#change to be able to specify amount of memory to use per analysis?
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
	elif args==[]:
		DoError('No input files selected!')
	elif options.maxinsertsize>100000 or options.maxinsertsize<10:
		DoError('Maximum insert size (-i) must be between 10 and 100,000')
	elif options.maxinsertsize<options.mininsertsize:
		DoError('Minimum insert size (-j) must be smaller than maximum insert size (-i). Currently -i='+str(options.maxinsertsize)+' and -j='+str(options.mininsertsize))
	elif options.mininsertsize>10000 or options.mininsertsize<10:
		DoError('Minimum insert size (-j) must be between 10 and 10,000')
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



####################
# Interactive Menu #
####################

def menusystem(options, args):
	
	run=False
	
	while run==False:
	
		os.system('clear')
		
		print "\nmultiple_mappings_to_bam.py: Written by Simon R Harris, Wellcome Trust Sanger Institute, UK. 2010"
		
		print "\nINPUT OPTIONS:"
		
		if options.ref=='':
			print "r: Reference dna sequence:\t\tNone selected (required)"
		else:
			print "r: Reference dna sequence:\t\t"+ref
			if len(args)==1:
				print len(args), "file to be mapped"
			else:
				print len(args), "files to be mapped"
				
			print "\nMAPPING OPTIONS:"
			print "p: Program:\t\t\t\t"+program
				
			if pairedend=='n':
				print "p: Use paired-end reads:\t\tno"
			else:
				print "p: Use paired-end reads:\t\tyes"
				if program=='maq':
					print "i: Maximum insert size:\t\t\t"+str(maxinsertsize)
					
				else:
					print "i: Maximum insert size:\t\t\t"+str(maxinsertsize)
					print "j: Minimum insert size:\t\t\t"+str(mininsertsize)
					
			if program=='maq':
				print "n: Maximum mismatches:\t\t\t"+str(mismatches)
				print "m: Maximum SNPs per read:\t\t"+str(maxsnps)
				print "d: Minimum mapping depth:\t\t"+str(mindepth)
				print "q: Minimum mapping quality:\t\t"+str(quality)
			elif program=='ssaha':
				print "q: Minimum mapping score:\t\t"+str(quality)
			
#			if velvet=='n':
#				print "v: Assemble non-mapping reads:\t\tno"
#			else:
#				print "v: Assemble non-mapping reads:\t\tyes"
			print "D: Reset to default mapping parameters"
			
		print "\nQ: QUIT"
		
		if ref=="":
			message="\nPlease select an option:"
			inputlist=['r', 'Q']
		else:
			message="\nPlease select an option or type y to run:"
			if program=='ssaha':
				inputlist=['r','q','v','y','t','Q','D','p','R','P']
			else:
				inputlist=['r','q','d','m','n','v','y','Q','D','p','P', 'R']
			if pairedend=='y':
				inputlist=inputlist+['i']
				if program=='ssaha':
					inputlist=inputlist+['j']
				
		ui=''
		while ui not in inputlist:
			ui=raw_input(message+' ')
	
		if ui=='y':
			os.system('clear')
			run=True
			
		elif ui=='r':
			ref=''
			while not os.path.isfile(ref):
				ref=raw_input('Enter reference file name including path or Q to go back to the menu: ')
				if ref=='Q':
					break
				elif not os.path.isfile(ref):
					print "File not found"
			options.ref=ref
			
		
		elif ui=='q':
			options.quality=0
			while options.quality > 100 or options.quality < 1:
				options.quality=int(raw_input('Enter minimum mapping quality (1-100): '))
		
#		elif ui=='v':
#			if velvet=='n':
#				velvet='y'
#			else:
#				velvet='n'
		
		elif ui=='o':
			outfile=''
			while outfile=='':
				outfile=raw_input('Enter prefix for output file names, or D to use the default: ')
				if outfile!='D' and outfile!='':
					outorig='n'
				elif outfile=='D':
					outorig='y'
		
		elif ui=='P':
			if program=='ssaha':
				program='bwa'
			elif program=='bwa':
				program='maq'
			else:
				program='ssaha'
		
		elif ui=='t':
			if rtype=='solexa':
				rtype='454'
			else:
				rtype='solexa'
		
		elif ui=='D':
			quality=30
			mindepth=5
			maxsnps=2
		
		elif ui=='p':
			if pairedend=='n':
				pairedend='y'
			else:
				pairedend='n'
		
		elif ui=='i':
			maxinsertsize=0
			while maxinsertsize > 10000 or maxinsertsize < 10 or maxinsertsize<=mininsertsize:
				maxinsertsize=int(raw_input('Enter maximum insert size (10-10,000). Must be more than min: '))
		
		elif ui=='j' and program=='ssaha':
			mininsertsize=0
			while mininsertsize > 10000 or mininsertsize < 10 or maxinsertsize<=mininsertsize:
				mininsertsize=int(raw_input('Enter minimum insert size (10-10,000). Must be less than max: '))
		
		elif ui=='R' and rtype=='solexa':
			readlength=0
			while readlength > 1000 or readlength < 10:
				readlength=int(raw_input('Enter read length (10-1000): '))
				
		elif ui=='Q':
			sys.exit()
	
#	if outorig=='y':
#		outfile=ref.split('.')[0]+"_q"+str(quality)+"_d"+str(mindepth)
		
	ref, inputdirs, quality, velvet, mindepth, maxsnps, pairedend, maxinsertsize, mininsertsize, mismatches, program, rtype, readlength, multiplexnamed, highmem=menusystem(ref, inputdirs, quality, velvet, mindepth, maxsnps, pairedend, maxinsertsize, mininsertsize, mismatches, program, rtype, readlength, multiplexnamed, highmem)

	return ref, inputdirs, quality, velvet, mindepth, maxsnps, pairedend, maxinsertsize, mininsertsize, mismatches, program, rtype, readlength, multiplexnamed, highmem






#####################
# SNPanalysis class #
#####################

class SNPanalysis:
	def __init__(self, fastq='', name='', mapped={}, runssaha='n', CDSseq='', number=0):
		self.fastq=fastq
		self.name=name
		self.runname=''
		self.fastqdir=''
		self.number=number
		self.pairedend=True
	
	def runSsaha(self, ref, bashfile):	
		print "\nRunning Ssaha on "+self.name+'...',
		sys.stdout.flush()
		#meaninsert=((options.maxinsertsize-options.mininsertsize)/2)+options.mininsertsize
		#Ssaha commands.
		
		#single end mapping
		if not self.pairedend:
			print >> bashfile, SSAHA_DIR+"ssaha2 -score "+str(options.ssahaquality)+" -kmer 13 -skip 2 -seeds 2 -score 12 -cmatch 9 -ckmer 6 -diff 0 -output sam_soft -outfile "+self.runname+"/tmp1.sam "+ref+" "+self.fastqdir+self.name+".fastq"
						
			
		#paired end
		else:
			#print >> bashfile, SSAHA_DIR+"ssaha2 -solexa -score "+str(options.quality)+" -skip 2 -diff 0 -kmer 13 -outfile "+self.runname+"/tmp.sam -multi 12345 -pair "+str(options.mininsertsize)+","+str(options.maxinsertsize)+" -output sam_soft "+ref+" "+self.fastqdir+self.name+"_1.fastq "+self.fastqdir+self.name+"_2.fastq"
			print >> bashfile, SSAHA_DIR+"ssaha2 -score "+str(options.ssahaquality)+" -kmer 13 -skip 2 -seeds 2 -score 12 -cmatch 9 -ckmer 6 -diff 0 -outfile "+self.runname+"/tmp1.sam -pair "+str(options.mininsertsize)+","+str(options.maxinsertsize)+" -output sam_soft "+ref+" "+self.fastqdir+self.name+"_1.fastq "+self.fastqdir+self.name+"_2.fastq"
		
		#print >> bashfile, SAMTOOLS_DIR+"samtools view -b -q "+str(options.mapq)+" -S",self.runname+"/tmp1.sam -t "+ref+".fai >", self.runname+"/tmp.bam"
		print >> bashfile, SAMTOOLS_DIR+"samtools view -b -S",self.runname+"/tmp1.sam -t "+ref+".fai >", self.runname+"/tmp1.bam"
		
		if self.pairedend and options.circular:
			print >> bashfile, MY_SCRIPTS_DIR+"fix_circular_bams.py -b", self.runname+"/tmp1.bam -o", self.runname+"/tmp"
			print >> bashfile, "rm", self.runname+"/tmp1.bam"
		else:
			print >> bashfile, "mv", self.runname+"/tmp1.bam", self.runname+"/tmp.bam"
			
		print >> bashfile, SAMTOOLS_DIR+"samtools view -H ",self.runname+"/tmp.bam >", self.runname+"/tmp2.sam"
#		
		print >> bashfile, "cat", self.runname+"/tmp2.sam", self.runname+"/tmp1.sam > ", self.runname+"/tmp.sam"
#		
		print >> bashfile, "rm", self.runname+"/tmp2.sam", self.runname+"/tmp1.sam"
		#print >> bashfile, SAMTOOLS_DIR+"samtools reheader", self.runname+"/tmp1.sam", self.runname+"/tmp.bam"
	
	
	def runBWA(self, ref, bashfile):	
		print "\nRunning BWA on "+self.name+'...',
		sys.stdout.flush()
		
		#Map the reads against the genome
		if self.pairedend:
			print >> bashfile, BWA_DIR+"bwa aln -q 15 ", options.ref, self.fastqdir+self.name+"_1.fastq >", self.runname+"/tmp.F.sai"
			print >> bashfile, BWA_DIR+"bwa aln -q 15 ", options.ref, self.fastqdir+self.name+"_2.fastq >", self.runname+"/tmp.R.sai"
			#Join both aligments
			print >> bashfile, BWA_DIR+"bwa sampe", ref,  self.runname+"/tmp.F.sai", self.runname+"/tmp.R.sai", self.fastqdir+self.name+"_1.fastq", self.fastqdir+self.name+"_2.fastq >", self.runname+"/tmp.sam"
		else:
			print >> bashfile, BWA_DIR+"bwa aln -q 15 ", options.ref, self.fastqdir+self.name+".fastq >", self.runname+"/tmp.F.sai"
			print >> bashfile, BWA_DIR+"bwa samse", ref,  self.runname+"/tmp.F.sai", self.fastqdir+self.name+".fastq >", self.runname+"/tmp.sam"
		
		
		#produce the BAM file
		#print >> bashfile, SAMTOOLS_DIR+"samtools view -b -q "+str(options.mapq)+" -S", self.runname+"/tmp.sam >", self.runname+"/tmp.bam"
		print >> bashfile, SAMTOOLS_DIR+"samtools view -b -S",self.runname+"/tmp.sam -t "+ref+".fai >", self.runname+"/tmp1.bam"
		if self.pairedend and options.circular:
			print >> bashfile, MY_SCRIPTS_DIR+"fix_circular_bams.py -b", self.runname+"/tmp1.bam -o", self.runname+"/tmp"
			print >> bashfile, "rm", self.runname+"/tmp1.bam"
		else:
			print >> bashfile, "mv", self.runname+"/tmp1.bam", self.runname+"/tmp.bam"
			
	def runSMALT(self, ref, bashfile):	
		
		if newsmalt:
			smaltoutput="bam"
		else:
			smaltoutput="sam"
		#Map the reads against the genome
		if self.domapping:
			print "\nRunning SMALT on "+self.name+'...',
			sys.stdout.flush()
			if self.pairedend:
				if options.maprepeats:
					print >> bashfile, SMALT_DIR+" map -y "+str(options.nomapid)+" -x -r 0 -i", options.maxinsertsize, " -j", options.mininsertsize, " -f "+smaltoutput+" -o "+self.runname+"/tmp1."+smaltoutput, tmpname+".index", self.fastqdir+self.name+"_1.fastq", self.fastqdir+self.name+"_2.fastq"
					cmdline="map -y "+str(options.nomapid)+" -x -r "+str(randrange(1,99999))+" -i", options.maxinsertsize, " -j", options.mininsertsize, " -f "+smaltoutput+"  -o "+self.runname+"/tmp1."+smaltoutput, tmpname+".index", self.fastqdir+self.name+"_1.fastq", self.fastqdir+self.name+"_2.fastq"
				else:
					if newsmalt:
						rbit=" -r -1"
					else:
						rbit=""
					print >> bashfile, SMALT_DIR+" map -y "+str(options.nomapid)+rbit+" -x -i", options.maxinsertsize, " -j", options.mininsertsize, " -f "+smaltoutput+" -o "+self.runname+"/tmp1."+smaltoutput, tmpname+".index", self.fastqdir+self.name+"_1.fastq", self.fastqdir+self.name+"_2.fastq"
					cmdline="map -y "+str(options.nomapid)+rbit+" -x -i", options.maxinsertsize, " -j", options.mininsertsize, " -f "+smaltoutput+" -o "+self.runname+"/tmp1."+smaltoutput, tmpname+".index", self.fastqdir+self.name+"_1.fastq", self.fastqdir+self.name+"_2.fastq"
			else:
				if options.maprepeats:
					print >> bashfile, SMALT_DIR+" map -y "+str(options.nomapid)+" -x -r 0 -f "+smaltoutput+" -o "+self.runname+"/tmp1."+smaltoutput, tmpname+".index", self.fastqdir+self.name+".fastq"
					cmdline="map -y "+str(options.nomapid)+" -x -r "+str(randrange(1,99999))+" -f "+smaltoutput+" -o "+self.runname+"/tmp1."+smaltoutput, tmpname+".index", self.fastqdir+self.name+".fastq"
				else:
					if newsmalt:
						rbit="-r -1"
					else:
						rbit=""
					print >> bashfile, SMALT_DIR+" map -y "+str(options.nomapid)+rbit+" -x -f "+smaltoutput+" -o "+self.runname+"/tmp1."+smaltoutput, tmpname+".index", self.fastqdir+self.name+".fastq"
					cmdline="map -y "+str(options.nomapid)+rbit+" -x -f "+smaltoutput+" -o "+self.runname+"/tmp1."+smaltoutput, tmpname+".index", self.fastqdir+self.name+".fastq"
			if not newsmalt:
				print >> bashfile, SAMTOOLS_DIR+"samtools view -b -S",self.runname+"/tmp1.sam -t "+ref+".fai >", self.runname+"/tmp1.bam"
				print >> bashfile, "rm", self.runname+"/tmp1.sam"
		else:
			print >> bashfile, "cp ", self.bam, self.runname+"/tmp1.bam"
		
		
		#produce the BAM file
		#print >> bashfile, SAMTOOLS_DIR+"samtools view -b -q "+str(options.mapq)+" -S", self.runname+"/tmp1.sam -t "+ref+".fai >", self.runname+"/tmp.bam"
		
		
		if self.pairedend and options.circular:
			print >> bashfile, MY_SCRIPTS_DIR+"fix_circular_bams.py -b", self.runname+"/tmp1.bam -o", self.runname+"/tmp"
			print >> bashfile, "rm", self.runname+"/tmp1.bam"
		else:
			print >> bashfile, "mv", self.runname+"/tmp1.bam", self.runname+"/tmp.bam"
		
		if options.GATK:
			print >> bashfile, "mv", self.runname+"/tmp.bam", self.runname+"/"+self.name+".bam"
			print >> bashfile, SAMTOOLS_DIR+"samtools view -H -o ", self.runname+"/tmphead.sam", self.runname+"/"+self.name+".bam"
			now = datetime.datetime.now()
			now = now.replace(microsecond=0)
			print >> bashfile, 'echo "@RG\tID:'+self.name+'\tCN:Sanger\tDT:'+now.isoformat()+'\tPG:SMALT\tPL:ILLUMINA\tSM:'+self.name+'" >>', self.runname+"/tmphead.sam"
			if self.domapping:
				print >> bashfile, "smaltversion=$( "+SMALT_DIR+" version  | grep Version | awk '{print $2}' )"
				print >> bashfile, 'echo "@PG\tID:SMALT\tPN:SMALT\tCL:'+' '.join(map(str,cmdline))+'\tVN:$smaltversion" >>', self.runname+'/tmphead.sam'
			print >> bashfile, SAMTOOLS_DIR+'samtools view -b -H -o', self.runname+'/tmphead.bam', self.runname+"/"+self.name+".bam"
			print >> bashfile, SAMTOOLS_DIR+'samtools merge -h ', self.runname+'/tmphead.sam -r', self.runname+"/tmp1.bam", self.runname+"/"+self.name+".bam", self.runname+'/tmphead.bam'
			print >> bashfile, "rm", self.runname+"/"+self.name+".bam"
			print >> bashfile, SAMTOOLS_DIR+'samtools sort', self.runname+"/tmp1.bam", self.runname+"/tmpsort"
			print >> bashfile, 'mv', self.runname+"/tmpsort.bam", self.runname+"/tmp1.bam"
			print >> bashfile, SAMTOOLS_DIR+'samtools index', self.runname+"/tmp1.bam"
			print >> bashfile, "cp", ref, self.runname+'/tmpref.fa'
			print >> bashfile, SAMTOOLS_DIR+'samtools faidx', self.runname+"/tmpref.fa"
			if options.mem>0:
				javamem=options.mem
			else:
				javamem=2
			print >> bashfile, JAVA_DIR+"java -Xmx"+str(javamem)+"g -jar", GATK_LOC, "-I", self.runname+"/tmp1.bam  -R", self.runname+"/tmpref.fa -T RealignerTargetCreator -o", self.runname+'/tmp.intervals'
			print >> bashfile, JAVA_DIR+"java -Xmx"+str(javamem)+"g -jar", GATK_LOC, "-I", self.runname+"/tmp1.bam  -R", self.runname+"/tmpref.fa -T IndelRealigner -targetIntervals", self.runname+'/tmp.intervals', "-o", self.runname+"/tmp.bam"
			print >> bashfile, "rm", self.runname+"/tmp1.bam", self.runname+"/tmp1.bam.bai",  self.runname+"/tmpref.*", self.runname+"/tmp.intervals", self.runname+"/tmphead.*"
		
		if options.plots:
			#add header to sam file for making plots - no need. Can read bams too!

			print >> bashfile, SAMTOOLS_DIR+"samtools view -h ",self.runname+"/tmp.bam >", self.runname+"/tmp.sam"
		else:
			print >> bashfile, "rm", self.runname+"/tmp1.sam"
		
			#print >> bashfile, "cat", self.runname+"/tmp2.sam", self.runname+"/tmp1.sam > ", self.runname+"/tmp.sam"
			#		
			#print >> bashfile, "rm", self.runname+"/tmp2.sam", self.runname+"/tmp1.sam"
			#print >> bashfile, SAMTOOLS_DIR+"samtools reheader", self.runname+"/tmp1.sam", self.runname+"/tmp.bam"

		
	def makepileup_from_sam(self, ref, bashfile):
	
		#create plots from sam file
		if options.plots:
			print >> bashfile, MY_SCRIPTS_DIR+"odd_plots_from_sam.py", self.runname+"/tmp.sam", self.runname+"/"+self.name+" >  /dev/null 2>&1"
			print >> bashfile, "rm", self.runname+"/tmp1.sam"
		
		#remove the sam file as it is no longer needed
		print >> bashfile, "rm", self.runname+"/tmp.sam"
		
		#order the bam file
		print >> bashfile, SAMTOOLS_DIR+"samtools sort ", self.runname+"/tmp.bam", self.runname+"/tmp1"
		
		#filter the bam file if requested
		if options.filter=="1":
			print >> bashfile, "mv", self.runname+"/tmp1.bam", self.runname+"/"+self.name+".bam"
		elif options.filter=="2":
			print >> bashfile, SAMTOOLS_DIR+"samtools view -F 4 -b -o", self.runname+"/"+self.name+".bam", self.runname+"/tmp1.bam"
		elif options.filter=="3":
			print >> bashfile, SAMTOOLS_DIR+"samtools view -f 2 -b -o", self.runname+"/"+self.name+".bam", self.runname+"/tmp1.bam"
		elif options.filter=="4":
			print >> bashfile, SAMTOOLS_DIR+"samtools view -F 4 -b -o", self.runname+"/"+self.name+".bam", self.runname+"/tmp1.bam"
			print >> bashfile, SAMTOOLS_DIR+"samtools view -f 4 -b -o", self.runname+"/"+self.name+"_unmapped.bam", self.runname+"/tmp1.bam"
		elif options.filter=="5":
			print >> bashfile, SAMTOOLS_DIR+"samtools view -f 2 -b -o", self.runname+"/"+self.name+".bam", self.runname+"/tmp1.bam"
			print >> bashfile, SAMTOOLS_DIR+"samtools view -F 2 -b -o", self.runname+"/"+self.name+"_unpaired.bam", self.runname+"/tmp1.bam"
		
		
		print >> bashfile, "rm", self.runname+"/tmp.bam", self.runname+"/tmp1.bam"	
		
		#index the bam file, to get the bai file.
		print >> bashfile, SAMTOOLS_DIR+"samtools index",  self.runname+"/"+self.name+".bam"
		 
		
		#produce the pileup file
		
		if options.BAQ:
			print >> bashfile, SAMTOOLS_DIR+"samtools mpileup -d 1000 -m", options.depth, " -DSugBf ", ref, self.runname+"/"+self.name+".bam >", self.runname+"/tmp.mpileup"
		else:
			print >> bashfile, SAMTOOLS_DIR+"samtools mpileup -d 1000 -m", options.depth, " -DSugf ", ref, self.runname+"/"+self.name+".bam >", self.runname+"/tmp.mpileup"
			
		
		print >> bashfile, BCFTOOLS_DIR+"bcftools view -bcg", self.runname+"/tmp.mpileup >", self.runname+"/"+self.name+".bcf"
		
		print >> bashfile, BCFTOOLS_DIR+"bcftools index", self.runname+"/"+self.name+".bcf"
		
		print >> bashfile, BCFTOOLS_DIR+"bcftools view -bcgv", self.runname+"/tmp.mpileup >", self.runname+"/"+self.name+"_variant.bcf"
		
		print >> bashfile, BCFTOOLS_DIR+"bcftools index", self.runname+"/"+self.name+"_variant.bcf"
		
		print >> bashfile, "rm", self.runname+"/tmp.mpileup"
		
		if options.indels:
			print >> bashfile, "mv ", self.runname+"/"+self.name+".bam", self.runname+"/tmp.bam"
			print >> bashfile, "mv ", self.runname+"/"+self.name+".bam.bai", self.runname+"/tmp.bam.bai"
			if options.ISfasta!="":
				print >> bashfile, MY_SCRIPTS_DIR+"realign_indels_test.py -b", self.runname+"/tmp.bam -B", self.runname+"/"+self.name+"_variant.bcf -r", ref, "-p", options.ratio, "-d", options.depth, "-o", self.runname+"/"+self.name, "-D", self.runname, "-I", options.ISfasta
			else:
				print >> bashfile, MY_SCRIPTS_DIR+"realign_indels_test.py -b", self.runname+"/tmp.bam -B", self.runname+"/"+self.name+"_variant.bcf -r", ref, "-p", options.ratio, "-d", options.depth, "-o", self.runname+"/"+self.name, "-D", self.runname
			
			if options.BAQ:
				print >> bashfile, SAMTOOLS_DIR+"samtools mpileup -d 1000 -m", options.depth, " -DSugBf ", ref, self.runname+"/"+self.name+".bam >", self.runname+"/tmp.mpileup"
			else:
				print >> bashfile, SAMTOOLS_DIR+"samtools mpileup -d 1000 -m", options.depth, " -DSugf ", ref, self.runname+"/"+self.name+".bam >", self.runname+"/tmp.mpileup"
				
			
			print >> bashfile, BCFTOOLS_DIR+"bcftools view -bcg", self.runname+"/tmp.mpileup >", self.runname+"/"+self.name+".bcf"
			
			print >> bashfile, BCFTOOLS_DIR+"bcftools index", self.runname+"/"+self.name+".bcf"
			
			print >> bashfile, BCFTOOLS_DIR+"bcftools view -bcgv", self.runname+"/tmp.mpileup >", self.runname+"/"+self.name+"_variant.bcf"
			
			print >> bashfile, BCFTOOLS_DIR+"bcftools index", self.runname+"/"+self.name+"_variant.bcf"
			
		# clean up:
		if not options.dirty:
			print >> bashfile, "rm", self.runname+"/tmp.*"

		#produce pseudosequence if requested
		#if options.pseudosequence==True:
		#print MY_SCRIPTS_DIR+"samtools_pileup_2_pseudosequence.py -p", self.runname+"/"+self.name+".pileup", "-b", self.runname+"/"+self.name+".bam", "-r", options.ratio, "-q", options.snpquality, "-o", self.runname+"/"+self.name
		if options.pseudosequence:
			print >> bashfile, MY_SCRIPTS_DIR+"bcf_2_pseudosequence_new_noindels.py -b ", self.runname+"/"+self.name+".bcf", "-B ", self.runname+"/"+self.name+".bam", "-r ", options.ratio, "-d ", options.depth, "-D ", options.stranddepth, "-q ", options.quality, "-m ", options.mapq, "-o", self.runname+"/"+self.name
			#print >> bashfile, 'cat '+self.runname+"/"+self.name+'.dna >> '+options.output+".aln"
			#print >> bashfile, 'gzip -f '+self.runname+"/"+self.name+'.dna '+self.runname+"/"+self.name+'.mfa '+self.runname+"/"+self.name+'.pileup '+self.runname+"/"+self.name+'*.plot'
		if options.plots:
			print >> bashfile, 'gzip -f '+self.runname+"/"+self.name+'*.plot'
		if not options.LSF:
			print >> bashfile, MY_SCRIPTS_DIR+'heterozygosity_plot.py -b', self.runname+"/"+self.name+".bcf -o", self.runname+"/"+self.name+"_contamination_plot.pdf", "-r", options.ratio, "-d", options.depth, "-D", options.stranddepth, "-q", options.quality
		
	
	
		
########
# Main #
########


if __name__ == "__main__":



	(options, args)=get_user_options()
	#print options, args
	check_input_validity(options, args)
	
	if options.version=="latest" or options.version=="0.7.4":
		SMALT_DIR="/nfs/users/nfs_s/sh16/smalt-0.7.4/smalt_x86_64"
		newsmalt=True
	elif options.version=="0.6.4":
		SMALT_DIR="/nfs/users/nfs_s/sh16/smalt-0.6.4/smalt_x86_64"
		newsmalt=False
	elif options.version=="0.6.3":
		SMALT_DIR="/nfs/users/nfs_s/sh16/smalt-0.6.3/smalt_x86_64"
		newsmalt=False
	elif options.version=="0.5.8":
		SMALT_DIR="/nfs/users/nfs_s/sh16/smalt-0.5.8/smalt_x86_64"
		newsmalt=False
	else:
		print "Unknown smalt version"
		sys.exit()


	print '\nChecking input files...'
	sys.stdout.flush()
	
	#make random name for files
	chars = string.ascii_letters + string.digits
	tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))
	
	pools=[]
	count=0
	poolsort=[]
	ziplist={}
	
	bamlist={}
	
	for pool in args:
		
		pairedend=options.pairedend
		if not os.path.isfile(pool):
			print "File "+pool+" not found! Skipping..."
			continue
		
		filetype='.'+pool.split('.')[-1]
		if filetype not in ['.fastq','.bam',".gz"]:
			print "WARNING: Input file name is not .fastq or .bam!"
		
		bam=""
		originalfastqdir=''
		if pool[-1]=='/':
			pool=pool[:-1]
		if len(pool.split('/'))>1:
			originalfastqdir='/'.join(pool.split('/')[:-1])+'/'
		if "_nonhuman" in pool:
			nonhumanpool="_nonhuman"
			pool=pool.replace("_nonhuman","")
		else:
			nonhumanpool=""
		
		if pool.split('.')[-1]=="gz" and pool.split('.')[-2]=="fastq":
			if not os.path.isdir(tmpname+"_unzipped"):
				os.system("mkdir "+tmpname+"_unzipped")
			#print "Unzipping ", pool
			
			#os.system("zcat "+pool+" > "+tmpname+"_unzipped/"+'.'.join(pool.split('/')[-1].split('.')[:-1]))
			#pool='.'.join(pool.split('.')[:-1])

			if options.pairedend:
				if pool.split('.')[-3][-2:] in ["_1", "_2"]:
					if not '.'.join(pool.split('/')[-1].split('.')[:-2])[:-2] in ziplist:
						ziplist['.'.join(pool.split('/')[-1].split('.')[:-2])[:-2]]=[pool]
					else:
						ziplist['.'.join(pool.split('/')[-1].split('.')[:-2])[:-2]].append(pool)
				else:
					ziplist['.'.join(pool.split('/')[-1].split('.')[:-2])]=[pool]
			else:
				ziplist['.'.join(pool.split('/')[-1].split('.')[:-2])]=[pool]
			
			pool=tmpname+"_unzipped/"+'.'.join(pool.split('/')[-1].split('.')[:-1])
			
			#ziplist.append(pool)
			
			
			
		
		elif pool.split('.')[-1]=="bam":
			
			if options.domapping:
				if not os.path.isdir(tmpname+"_unbammed"):
					os.system("mkdir "+tmpname+"_unbammed")
				bamlist['.'.join(pool.replace("#","_").split('/')[-1].split('.')[:-1])]=pool
				
				pool=tmpname+"_unbammed/"+'.'.join(pool.replace("#","_").split('/')[-1].split('.')[:-1])+".bam"
			else:
				bam=pool
			
	
#		if not os.path.isfile(pool):
#			print "File "+pool+" not found! Skipping..."
#			continue
		fastqdir=''
		if pool[-1]=='/':
			pool=pool[:-1]
		if len(pool.split('/'))>1:
			fastqdir='/'.join(pool.split('/')[:-1])+'/'
		
		
		pool='.'.join(pool.split('/')[-1].split('.')[:-1])
		

		if options.pairedend and filetype!=".bam":
			if pool[-2:]=='_1':
				
				if not os.path.isfile(originalfastqdir+pool[:-2]+"_2"+nonhumanpool+".fastq") and not os.path.isfile(originalfastqdir+pool[:-2]+"_2"+nonhumanpool+".fastq.gz"):
					print "File "+pool+"_2.fastq not found! Treating "+pool+" as unpaired..."
					pairedend=False
				else:
					pool=pool[:-2]
			elif pool[-2:]=='_2':
				if not os.path.isfile(originalfastqdir+pool[:-2]+"_1"+nonhumanpool+".fastq") and not os.path.isfile(originalfastqdir+pool[:-2]+"_1"+nonhumanpool+".fastq.gz"):
					print "File "+pool+"_1.fastq not found! Treating "+pool+" as unpaired..."
					pairedend=False
				else:
					pool=pool[:-2]
			else:
				print "Not a typical paired-end name format! Treating "+pool+" as unpaired..."
				pairedend=False
			
		
		name=pool

		if options.diroutput=="":
			pool=pool+"_"+options.program
		else:
			pool=pool+"_"+options.diroutput

		if pool in poolsort:
			continue
		print pool+'...',
		sys.stdout.flush()

		if not os.path.isdir(pool):
		 	print "pool "+pool+" not found! Creating...",
			os.system("mkdir "+pool)
						
		pools.append(SNPanalysis())
		pools[count].number=str(count+1)
		pools[count].runname=pool
		pools[count].name=name
		pools[count].fastqdir=fastqdir
		pools[count].filetype=filetype
		pools[count].pairedend=pairedend
		if bam!="":
			pools[count].domapping=False
			pools[count].bam=bam
			print "As you chose not to remap, using", bam, "as mapped data"
		else:
			pools[count].domapping=True
#		the rest can be globals
#		pools[count].quality=options.qualityoptions
#		pools[count].mininsertsize=options.mininsertsize
#		pools[count].maxinsertsize=options.maxinsertsize
#		pools[count].pairedend=options.pairedend
#		pools[count].readlength=options.readlength
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
		
	
	if not options.human:
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
		if options.human:
			os.system(SMALT_DIR+" index -k 20 -s 13 "+tmpname+".index "+options.ref)
		else:
			os.system(SMALT_DIR+" index -k 13 -s 1 "+tmpname+".index "+options.ref)
	#elif options.program=='ssaha':
	#	os.system("samtools faidx "+options.ref)
	os.system("samtools faidx "+options.ref)
		
	count=0

	for pool in pools:
		
		if options.keep and options.pseudosequence and os.path.isfile(pool.runname+"/"+pool.name+".bam") and os.path.isfile(pool.runname+"/"+pool.name+".bcf") and os.path.isfile(pool.runname+"/"+pool.name+".mfa"):# and os.path.isfile(pool.runname+"/"+pool.name+"_indels.txt"):
			print "\n"+pool.runname, "contains a .bam, .bcf and .mfa file already. As you selected the keep option, it will not be mapped again.",
		elif options.keep and os.path.isfile(pool.runname+"/"+pool.name+".bam") and os.path.isfile(pool.runname+"/"+pool.name+".bcf"):# and os.path.isfile(pool.runname+"/"+pool.name+"_indels.txt"):
			print "\n"+pool.runname, "contains a .bam, and .bcf file already. As you selected the keep option, it will not be mapped again.",
		else:
			count+=1
			bashfile=open(str(count)+tmpname+'_sbs.sh','w')

			if pool.name in ziplist:
				for filename in ziplist[pool.name]:
					print >> bashfile, "zcat "+filename+" > "+pool.fastqdir+'.'.join(filename.split('/')[-1].split('.')[:-1])
			
			
			elif pool.name in bamlist:
				if options.pairedend:
					print >> bashfile, MY_SCRIPTS_DIR+'bam_filter.py -t all -b '+bamlist[pool.name]+' -o '+pool.fastqdir+pool.name
				else:
					print >> bashfile, MY_SCRIPTS_DIR+'bam_filter.py -t all -f fastq -b '+bamlist[pool.name]+' -o '+pool.fastqdir+pool.name
			
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
			argstring.append(pool.runname+"/"+pool.name+".mfa")
		if options.indels:
			joinstring=MY_SCRIPTS_DIR+"join_dna_files_with_indels_newb.py -r "+options.ref+" -o "+options.output+".aln "+' '.join(argstring) #*_ssaha/*_test.mfa
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
				os.system('echo \'bash ${LSB_JOBINDEX}'+tmpname+'_sbs.sh\' | bsub -R \'select[mem>'+memmb+'] rusage[mem='+memmb+']\' -q '+options.LSFQ+'  -J'+tmpname+'_'+options.program+'"[1-'+str(count)+']%'+str(options.nodes)+'"  -M '+memkb+' -o '+tmpname+options.program+'-%I.out -e '+tmpname+options.program+'-%I.err > '+tmpname+'jobid')# run all ssaha jobs in job array . add this to exclude a node: -R \'hname!=pcs4k\'
			else:
				os.system('echo \'bash ${LSB_JOBINDEX}'+tmpname+'_sbs.sh\' | bsub -q '+options.LSFQ+' -J'+tmpname+'_'+options.program+'"[1-'+str(count)+']%'+str(options.nodes)+'"" -o '+tmpname+options.program+'-%I.out -e '+tmpname+options.program+'-%I.err > '+tmpname+'jobid')
			
			jobid=open(tmpname+'jobid', "rU").read().split(">")[0].split("<")[1]
			os.system('bsub -w \'ended('+tmpname+'_'+options.program+')\' \"bacct -e -l '+jobid+' > '+options.output+"_mapping_job_bacct.txt \"")
			os.system('bsub -w \'ended('+tmpname+'_'+options.program+')\' \"bhist -l '+jobid+' > '+options.output+"_mapping_job_bhist.txt \"")
			
			if not options.dirty:
					os.system('bsub -w \'ended('+tmpname+'_'+options.program+')\' rm -rf *'+tmpname+'_sbs.sh '+tmpname+'*.out '+tmpname+'*.err '+options.ref+'.pac '+options.ref+'.ann '+options.ref+'.rpac '+options.ref+'.amb '+options.ref+'.rbwt '+options.ref+'.bwt '+options.ref+'.sa '+options.ref+'.rsa '+tmpname+'*.sma '+tmpname+'*.smi '+tmpname+'_unzipped '+tmpname+'_unbammed '+tmpname+'jobid')#when job array is all done delete
		
		
		if options.pseudosequence:
			if len(pools)>200:
				if count==0:
					os.system('echo '+joinstring+' | bsub -M 10000000 -q long -R \'select[mem>10000] rusage[mem=10000]\' -J'+tmpname+'_joining -o '+options.output+'_join.out -e '+options.output+'_join.err')
				else:
					os.system('echo '+joinstring+' | bsub -M 10000000 -q long -R \'select[mem>10000] rusage[mem=10000]\' -J'+tmpname+'_joining -w \'ended('+tmpname+'_'+options.program+')\' -o '+options.output+'_join.out -e '+options.output+'_join.err ')
				os.system('echo '+summarystring+' | bsub -M 16000000 -q long -R \'select[mem>16000] rusage[mem=16000]\' -w \'ended('+tmpname+'_joining)\' -o '+options.output+'_sum.out -e '+options.output+'_sum.err')
			else:
				if count==0:
					os.system('echo '+joinstring+' | bsub -M 2000000 -R \'select[mem>2000] rusage[mem=2000]\' -J'+tmpname+'_joining -o '+options.output+'_join.out -e '+options.output+'_join.err')
				else:
					os.system('echo '+joinstring+' | bsub -M 2000000 -R \'select[mem>2000] rusage[mem=2000]\' -J'+tmpname+'_joining -w \'ended('+tmpname+'_'+options.program+')\' -o '+options.output+'_join.out -e '+options.output+'_join.err ')
				os.system('echo '+summarystring+' | bsub -M 16000000 -q long -R \'select[mem>16000] rusage[mem=16000]\' -w \'ended('+tmpname+'_joining)\' -o '+options.output+'_sum.out -e '+options.output+'_sum.err')
			
		#os.system('bsub -w \'ended('+tmpname+'_'+program+')\' rm *'+tmpname+'_sbs.sh; for f in tmpname*.err; do if test ! -s $f;then rm $f ${f%.err}.out;fi;done')
	
	elif options.pseudosequence:
		os.system(joinstring)
		os.system(summarystring)

		

			
			
