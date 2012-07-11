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



####################
# Set some globals #
####################

SAMTOOLS_DIR=""
SSAHA_DIR=""
BWA_DIR=""
MY_SCRIPTS_DIR=""


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
	usage = "usage: %prog [options] <list of fastq files>"
	version="%prog 1.0. Written by Simon Harris, Wellcome Trust Sanger Institute, 2010"
	parser = OptionParser(usage=usage, version=version)
	
	#do not allow arguments to be interspersed. e.g. -a -b arg1 agr2. MUST be -a arg1 -b arg2.
	parser.disable_interspersed_args()
	
	group = OptionGroup(parser, "Required Options")
	group.add_option("-r", "--reference", action="store", dest="ref", help="reference DNA sequence (in fasta or multi-fasta format)", default="", metavar="FILE")
	parser.add_option_group(group)

	group = OptionGroup(parser, "Mapping Options")
	group.add_option("-p", "--program", action="store", type="choice", dest="program", choices=["bwa","ssaha"], help="Mapping program to use (choose from bwa or ssaha) [default= %default]", default="bwa")
	group.add_option("-l", "--length", action="store", dest="readlength", help="Read length [default= %default]", default=54, type="int", metavar="INT")
	group.add_option("-s", "--single", action="store_false", dest="pairedend", help="reads are single ended (not paired)", default=True)
	group.add_option("-i", "--maxinsert", action="store", dest="maxinsertsize", help="maximum insert size [default= %default]", default=400, type="int", metavar="INT")
	group.add_option("-j", "--mininsert", action="store", dest="mininsertsize", help="minimum insert size [default= %default]", default=100, type="int", metavar="INT")
	group.add_option("-q", "--quality", action="store", type="int", dest="quality", help="Minimum mapping quality [default= %default]", default=30, metavar="INT")
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "Pseudosequence creation options")
	group.add_option("-d", "--dna", action="store_true", dest="pseudosequence", help="Create pseudosequences", default=False)
	group.add_option("-S", "--SNPquality", action="store", type="int", dest="snpquality", help="Minimum site mapping quality for SNP calling [default= %default]", default=30, metavar="INT")
	group.add_option("-R", "--ratio", action="store", type="float", dest="ratio", help="SNP/site mapping quality ratio cutoff [default= %default]", default=0.75, metavar="FLOAT")
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "Output options")
	group.add_option("-o", "--output", action="store", dest="output", help="output file prefix", default="")
	group.add_option("-a", "--align", action="store_true", dest="alnfile", help="Create snp alignment file (in phylip format)", default=False)
	group.add_option("-P", "--phylogeny", action="store_true", dest="raxml", help="Run phylogeny with RAxML", default=False)
	group.add_option("-m", "--model", action="store", dest="model", help="Model of evolution to use. [Default= %default]", default="GTRGAMMA", type="choice", choices=["GTRGAMMA","GTRGAMMAI", "GTRCAT", "GTRMIX", "GTRMIXI"])
	group.add_option("-b", "--bootstrap", action="store", dest="bootstrap", help="Number of bootstrap replicates (0 = do not run bootstrap). [Default= %default]", default=100, type="int", metavar="int")
	parser.add_option_group(group)
	
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):

	if options.output=='':
		options.output=options.ref.split("/")[-1].split(".")[0]+"_"+options.program


	while os.path.isfile(options.output+".aln"):
		outopt=""
		outopt=raw_input('\nOutput files with chosen prefix already exist.\n\nWould you like to overwrite (o), choose a new output file prefix (n) or quit (Q): ')
		if outopt=='Q':
			sys.exit()
		elif outopt=="o":
			break
		elif outopt=="n":
			options.output=raw_input('Enter a new output file prefix: ')
			

	if options.ref=='':
		DoError('No reference dna file (-r) selected!')
	elif args==[]:
		DoError('No input files selected!')
	elif options.maxinsertsize>10000 or options.maxinsertsize<10:
		DoError('Maximum insert size (-i) must be between 10 and 10,000!')
	elif options.maxinsertsize<options.mininsertsize:
		DoError('Minimum insert size (-j) must be smaller than maximum insert size (-i)! Currently -i='+str(options.maxinsertsize)+' and -j='+str(options.mininsertsize))
	elif options.mininsertsize>10000 or options.mininsertsize<10:
		DoError('Minimum insert size (-j) must be between 10 and 10,000!')
	elif options.quality>100 or options.quality<1:
		DoError('Mapping quality score (-q) must be between 1 and 100!')
	elif options.snpquality>100 or options.snpquality<1:
		DoError('Minimum site mapping quality for SNP calling (-q) must be between 1 and 100!')
	elif options.ratio>1 or options.ratio<0:
		DoError('SNP/site mapping quality ratio cutoff (-R) must be between 0 and 1!')
	elif options.readlength>1000 or options.readlength<36:
		DoError('Read length (-l) must be between 36 and 1000!')
	
	if options.raxml==True and options.pseudosequence==False:
		options.pseudosequence=True
		
	return




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
	
	def runSsaha(self, ref):	
		print "\nRunning Ssaha on "+self.name+'...',
		sys.stdout.flush()
		
		#single end mapping
		if options.pairedend==False:
			os.system(SSAHA_DIR+"ssaha2 -"+self.rtype+" -score "+str(options.quality)+" -skip 2 -diff 0 -output sam -outfile "+self.runname+"/tmp.sam "+ref+" "+self.fastqdir+self.name+".fastq")
						
			
		#paired end
		else:
			os.system(SSAHA_DIR+"ssaha2 -solexa -outfile "+self.runname+"/tmp1.sam -pair "+str(options.mininsertsize)+","+str(options.maxinsertsize)+" -output sam_soft "+ref+" "+self.fastqdir+self.name+"_1.fastq "+self.fastqdir+self.name+"_2.fastq")

		os.system(SAMTOOLS_DIR+"samtools faidx "+ref )
		
		os.system(SAMTOOLS_DIR+"samtools view -b -S "+self.runname+"/tmp1.sam -t "+ref+".fai > "+self.runname+"/tmp.bam")
			
		os.system(SAMTOOLS_DIR+"samtools view -H "+self.runname+"/tmp.bam > "+self.runname+"/tmp2.sam")
		
		os.system("cat "+self.runname+"/tmp2.sam "+self.runname+"/tmp1.sam > "+self.runname+"/tmp.sam")
		
		os.system("rm "+self.runname+"/tmp2.sam "+self.runname+"/tmp1.sam")
		
	
	
	def runBWA(self, ref):	
		print "\nRunning BWA on "+self.name+'...',
		sys.stdout.flush()
		
		#Map the reads against the genome
		os.system(BWA_DIR+"bwa aln -q 15 "+options.ref+" "+self.fastqdir+self.name+"_1.fastq > "+self.runname+"/tmp.F.sai")
		os.system(BWA_DIR+"bwa aln -q 15 "+options.ref+" "+self.fastqdir+self.name+"_2.fastq > "+self.runname+"/tmp.R.sai")
		
		#Join both aligments
		os.system(BWA_DIR+"bwa sampe "+ref+" "+self.runname+"/tmp.F.sai "+self.runname+"/tmp.R.sai "+self.fastqdir+self.name+"_1.fastq "+self.fastqdir+self.name+"_2.fastq > "+self.runname+"/tmp.sam")
		
		#produce the BAM file
		os.system(SAMTOOLS_DIR+"samtools view -b -S "+self.runname+"/tmp.sam > "+self.runname+"/tmp.bam")
		
		
	def makepileup_from_sam(self, ref):
		
		#remove the sam file as it is no longer needed
		os.system("rm "+self.runname+"/tmp.sam")
		
		#order the bam file
		os.system(SAMTOOLS_DIR+"samtools sort "+self.runname+"/tmp.bam "+self.runname+"/"+self.name)
		
		#index the bam file, to get the bai file.
		os.system(SAMTOOLS_DIR+"samtools index "+self.runname+"/"+self.name+".bam")
		
		#produce the pileup file
		os.system(SAMTOOLS_DIR+"samtools pileup -c -f "+ref+" "+self.runname+"/"+self.name+".bam > "+self.runname+"/"+self.name+".pileup")
		
		# clean up:
		os.system("rm "+self.runname+"/tmp*")

		#produce pseudosequence if requested
		if options.pseudosequence==True:
			os.system(MY_SCRIPTS_DIR+"samtools_pileup_2_pseudosequence.py -p "+self.runname+"/"+self.name+".pileup "+"-b "+self.runname+"/"+self.name+".bam -r "+str(options.ratio)+" -q "+str(options.snpquality)+" -o "+self.runname+"/"+self.name)
			os.system('cat '+self.runname+"/"+self.name+'.dna >> '+options.output+".aln")
		
########
# Main #
########


if __name__ == "__main__":



	(options, args)=get_user_options()
	#print options, args
	check_input_validity(options, args)


	print '\nChecking input files...'
	sys.stdout.flush()

	
	pools=[]
	count=0
	poolsort=[]
	for pool in args:

		if not os.path.isfile(pool):
			print "File "+pool+" not found! Skipping..."
			continue
		fastqdir=''
		if pool[-1]=='/':
			pool=pool[:-1]
		if len(pool.split('/'))>1:
			fastqdir='/'.join(pool.split('/')[:-1])+'/'
		filetype='.'+pool.split('.')[-1]
		
		if filetype!='.fastq':
			print "WARNING: Input file name is not .fastq!"
		
		pool='.'.join(pool.split('/')[-1].split('.')[:-1])
		pool=pool.replace("_nonhuman","")

		
		if options.pairedend==True:
			if pool[-2:]=='_1':
				pool=pool[:-2]
				if not os.path.isfile(fastqdir+pool+"_2.fastq"):
					print "File "+pool+"_2.fastq not found! Skipping..."
					continue
			elif pool[-2:]=='_2':
				pool=pool[:-2]
				if not os.path.isfile(fastqdir+pool+"_1.fastq"):
					print "File "+pool+"_1.fastq not found! Skipping..."
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
		pools[count].runname=pool
		pools[count].name=name
		pools[count].fastqdir=fastqdir
		pools[count].filetype=filetype
		
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
		concatenated_seq=""
		for record in seq_records:
			concatenated_seq=concatenated_seq+str(record.seq)
			
		my_seq_record = SeqRecord(Seq(concatenated_seq))
		my_seq_record.id=options.ref.split("/")[-1].split(".")[0]
		my_seq_record.description="Reference"
		
		SeqIO.write([my_seq_record], open(options.output+".aln","w"), "fasta")
		
		
	
	
	#Running Ssaha or bwa where required
	if options.program=='bwa':
		os.system("bwa index "+options.ref)
	elif options.program=='ssaha':
		os.system("samtools faidx "+options.ref)
		
	
		
	for pool in pools:
		if options.program=='bwa':
			pool.runBWA(options.ref)
			
		elif options.program=='ssaha':
			pool.runSsaha(options.ref)
		
			
		pool.makepileup_from_sam(options.ref)
		
			
	
	summarystring=MY_SCRIPTS_DIR+'summarise_snps.py -r '+options.ref.split("/")[-1].split(".")[0]+' -i '+options.output+'.aln'

	
	if options.alnfile:
		summarystring=summarystring+" -a"
	if options.raxml:
		summarystring=summarystring+" -p -m "+options.model+" -b "+str(options.bootstrap)
		
	
	os.system(summarystring)

		

			
			
			
			
