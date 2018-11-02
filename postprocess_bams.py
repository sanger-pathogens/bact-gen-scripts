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
#SMALT_DIR="/nfs/users/nfs_s/sh16/smalt-0.5.5/"
#SMALT_DIR=""
MY_SCRIPTS_DIR="/nfs/pathogen/sh16_scripts/"
BREAKDANCER_DIR="/nfs/users/nfs_s/sh16/breakdancer-1.1_2011_02_21/"
PINDEL_DIR="/nfs/users/nfs_s/sh16/pindel/trunk/"
DINDEL_DIR="/nfs/users/nfs_s/sh16/dindel/binaries/"
DINDEL_SCRIPTS_DIR="/nfs/users/nfs_s/sh16/dindel/dindel-1.01-python/"


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
	usage = "usage: %prog [options] <list of bam files>"
	version="%prog 1.0. Written by Simon Harris, Wellcome Trust Sanger Institute, 2011"
	parser = OptionParser(usage=usage, version=version)
	
	#do not allow arguments to be interspersed. e.g. -a -b arg1 agr2. MUST be -a arg1 -b arg2.
	parser.disable_interspersed_args()
	
	group = OptionGroup(parser, "Required Options")
	group.add_option("-r", "--reference", action="store", dest="ref", help="Reference DNA sequence (in fasta or multi-fasta format)", default="", metavar="FILE")
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "I/O Options")
	group.add_option("-s", "--studyid", action="store", dest="studyid", help="ID of study to analyse", default=0, metavar="INT", type="int")
	group.add_option("-p", "--passed", action="store_true", dest="passed", help="If using a study ID, only use manually passed runs", default=False)
	group.add_option("-e", "--embl", action="store", dest="embl", help="Reference annotation (in true embl or genbank format) [required for dN/dS calculations and CDS naming]", default="", metavar="FILE")
	group.add_option("-o", "--output", action="store", dest="output", help="Output file prefix/suffix", default="")
	group.add_option("-b", "--bcf", action="store_false", dest="bcf", help="Turn off bcf file creation", default=True)
	group.add_option("-I", "--indels", action="store_false", dest="indels", help="Turn off identification of potential indels with breakdancer, pindel and dindel", default=True)
	group.add_option("-X", "--dna", action="store_false", dest="pseudosequence", help="Turn off pseudosequence creation", default=True)
	group.add_option("-x", "--noref", action="store_false", dest="incref", help="Do not include reference in pseudosequence alignment", default=True)
	group.add_option("-g", "--plots", action="store_true", dest="plots", help="Create mapping plots", default=False)
	group.add_option("-t", "--tabfiles", action="store_true", dest="tabfile", help="Create tabfile of snps", default=False)
	group.add_option("-f", "--force", action="store_true", dest="force", help="Force overwrite of output files", default=False)
	parser.add_option_group(group)

	group = OptionGroup(parser, "Filtering Options")
	#group.add_option("-l", "--length", action="store", dest="readlength", help="Read length [default= %default]", default=54, type="int", metavar="INT")
	#group.add_option("-s", "--single", action="store_false", dest="pairedend", help="reads are single ended (not paired)", default=True)
	group.add_option("-i", "--maxinsert", action="store", dest="maxinsertsize", help="Maximum insert size (ssaha and smalt only) [default= %default]", default=600, type="int", metavar="INT")
	group.add_option("-j", "--mininsert", action="store", dest="mininsertsize", help="Minimum insert size (ssaha and smalt only) [default= %default]", default=0, type="int", metavar="INT")
	group.add_option("-F", "--properpairfilter", action="store_false", dest="properpairfilter", help="Do not filter out reads that don't map in a proper pair from the bam file (after indel detection) [default = Filter]", default=True)
	group.add_option("-q", "--quality", action="store", type="int", dest="quality", help="Minimum base call quality to call a SNP (see samtools help for more information) [default= %default]", default=50, metavar="INT")
	group.add_option("-Q", "--mapq", action="store", type="int", dest="mapq", help="Minimum mapping quality to call a SNP (see samtools help for more information) [default= %default]", default=30, metavar="INT")	
	group.add_option("-d", "--depth", action="store", dest="depth", help="Minimum number of reads matching SNP [default= %default]", default=4, type="int")
	group.add_option("-D", "--stranddepth", action="store", dest="stranddepth", help="Minimum number of reads matching SNP per strand [default= %default]", default=2, type="int")
	group.add_option("-B", "--BAQ", action="store_false", dest="BAQ", help="Turn off samtools base alignment quality option (BAQ) ", default=True)
	#parser.add_option("-q", "--quality", action="store", dest="quality", help="Minimum base quality [default= %default]", default=120, type="int")
	#group.add_option("-S", "--RMS", action="store", dest="RMS", help="Minimum root mean squared mapping quality [default= %default]", default=25, type="int")
	#parser.add_option("-Q", "--strandquality", action="store", dest="strandquality", help="Minimum per strand base quality [default= %default]", default=60, type="int")
	group.add_option("-R", "--ratio", action="store", dest="ratio", help="SNP/Mapping quality ratio cutoff [default= %default]", default=0.8, type="float")

	
	group = OptionGroup(parser, "General usage options")
	#group.add_option("-I", "--interactive", action="store_true", dest="interactive", help="Enter interactive menu system", default=True)
	group.add_option("-k", "--keep", action="store_true", dest="keep", help="If old output files are present, do not rerun them", default=False)
	group.add_option("-L", "--LSF", action="store_false", dest="LSF", help="Do not use LSF to parallelise analyses", default=True)
	group.add_option("-U", "--queue", action="store", dest="LSFQ", help="LSF queue to submit to. [Default= %default]", default="normal", type="choice", choices=["normal","long", "basement", "hugemem"])
	group.add_option("-M", "--memory", action="store", dest="mem", help="Amount of memory required for analysis (Gb). [Default= %default]", default=2, type="int")#change to be able to specify amount of memory to use per analysis?
	group.add_option("-y", "--dirty", action="store_true", dest="dirty", help="Do not clean up temporary files. [Default= %default]", default=False)#change to be able to specify amount of memory to use per analysis?
	parser.add_option_group(group)
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):

	if options.output=='':
		options.output=options.ref.split("/")[-1].split(".")[0]+"_bam_analysis"


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
	elif options.mapq>60 or options.mapq<0:
		DoError('Mapping quality score (-Q) must be between 0 and 60')
#	elif options.snpquality>100 or options.snpquality<1:
#		DoError('Minimum site mapping quality for SNP calling (-q) must be between 1 and 100!')
	elif options.ratio>1 or options.ratio<0:
		DoError('SNP/site mapping quality ratio cutoff (-R) must be between 0 and 1')
#	elif options.readlength>1000 or options.readlength<36:
#		DoError('Read length (-l) must be between 36 and 1000!')
	elif options.mem>30 or options.mem<0:
		DoError('Memory requirement (-M) must be between 0 and 30Gb')

	
	return




#####################
# SNPanalysis class #
#####################

class bamanalysis:
	def __init__(self, bamfile='', name='', mapped={}, runssaha='n', CDSseq='', number=0, log=[], keep=False):
		self.bamfile=bamfile
		self.name=name
		self.dirname=''
		self.originaldir=''
		self.number=number
		self.makebcf=False
		self.makemfa=False
		self.findindels=False
		self.log=[]
		self.keep=False
		
	def makepileup_from_sam(self, ref, bashfile):
		
		
		if not os.path.isfile(self.originaldir+self.name+".bam.bai"):
			print >> bashfile, SAMTOOLS_DOR+"samtools index", self.originaldir+self.name+".bam.bai"
		
		bamfilename=self.originaldir+self.name+".bam"
		
			
			
		
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
	
	bamfiles=[]
	count=0
	bamlist=[]
	inputbamfiles=[]
	
	if options.studyid!=0:
		print "\nFinding bam files in study "+str(options.studyid)+"..."
		sys.stdout.flush()
		pathfindstring="pathfind -t study -f fastq -id "+str(options.studyid)
		if options.passed:
			pathfindstring=pathfindstring+" -qc passed"
		pathfindarg = shlex.split(pathfindstring)
		returnval = subprocess.Popen(pathfindarg, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		inputbamfiles = returnval.stdout.read().split("\n")
		print inputbamfiles
		sys.exit()
		
	inputbamfiles=inputbamfiles+args
	
	print '\nChecking input files...'
	sys.stdout.flush()
	for bamfile in inputbamfiles:
		
		if not os.path.isfile(bamfile):
			print "File "+bamfile+" not found! Skipping..."
			continue
		
		filetype=bamfile.split('.')[-1]
		if filetype != "bam":
			print "WARNING: Input file name does not end in .bam!"
			continue
		
		originaldir=''
		if bamfile[-1]=='/':
			bamfile=bamfile[:-1]
		if len(bamfile.split('/'))>1:
			originaldir='/'.join(bamfile.split('/')[:-1])+'/'
		
		
		name='.'.join(bamfile.split('/')[-1].split('.')[:-1])
		
		dirname=name+"_"+options.output

		if bamfile in bamlist:
			continue
		print bamfile+'...',
		sys.stdout.flush()

		if not os.path.isdir(dirname):
		 	print "directory "+dirname+" not found. Creating...",
			os.system("mkdir "+dirname)
						
		bamfiles.append(bamanalysis())
		bamfiles[count].number=str(count+1)
		bamfiles[count].bamfile=bamfile
		bamfiles[count].name=name
		originaldir
		bamfiles[count].originaldir=originaldir
		bamfiles[count].dirname=dirname+"/"

	

		if options.bcf:
			bamfiles[count].makebcf=True
		else:
			bamfiles[count].makebcf=False
		
		if options.pseudosequence:
			bamfiles[count].makemfa=True
		else:
			bamfiles[count].makemfa=False
		
		if options.indels:
			bamfiles[count].findindels=True
		else:
			bamfiles[count].findindels=False
		
		if options.keep:
			bamfiles[count].keep=True
		else:
			bamfiles[count].keep=False
		
		
		bamlist.append(bamfile)
		print 'ok'
		sys.stdout.flush()
		
		count=count+1
		
		
	if len(bamlist)==0:
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

		
		else:
			os.system("rm "+options.output+".aln")
	

	os.system(SAMTOOLS_DIR+"samtools faidx "+options.ref)
		
	count=0

	for bamfile in bamfiles:
			
		if options.keep and os.path.isfile(bamfile.dirname+"/"+bamfile.name+".log"):
			bamfile.log=open(bamfile.dirname+"/"+bamfile.name+".log", "rU").readlines()
		if os.path.isfile(bamfile.dirname+"/"+bamfile.name+".log"):
			os.remove(bamfile.dirname+"/"+bamfile.name+".log")
		
		count+=1
		bashfile=open(str(count)+tmpname+'_sbs.sh','w')
			
		 
		bamfile.makepileup_from_sam(options.ref, bashfile)
	
	
		bashfile.close()
		if bamfile.keep:
			os.system("rm -f "+str(count)+tmpname+'_sbs.sh')
			count-=1
			continue


		if options.LSF==False:
			os.system('bash '+str(count)+tmpname+'_sbs.sh')
			if not options.dirty:
				os.system('rm '+str(count)+tmpname+'_sbs.sh '+options.ref+'.pac '+options.ref+'.ann '+options.ref+'.rpac '+options.ref+'.amb '+options.ref+'.rbwt '+options.ref+'.bwt '+options.ref+'.sa '+options.ref+'.rsa')
	
	if options.pseudosequence:
		argstring=[]
		for bamfile in bamfiles:	
			argstring.append(bamfile.dirname+"/"+bamfile.name+".mfa")
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
#		if options.alnfile:
#			summarystring=summarystring+" -a"
		if options.tabfile:
			summarystring=summarystring+" -t"
#		if options.raxml:
#			summarystring=summarystring+" -p -l -b "+str(options.bootstrap)
		
	print
	if options.LSF==True:
		if count>0:
			if options.mem>0:
				memkb=str(options.mem*1000000)
				memmb=str(options.mem*1000)
				os.system('echo \'bash ${LSB_JOBINDEX}'+tmpname+'_sbs.sh\' | bsub -R \'select[mem>'+memmb+'] rusage[mem='+memmb+']\' -q '+options.LSFQ+'  -J'+tmpname+'_'+'"[1-'+str(count)+']"  -M '+memkb+' -o '+tmpname+'-%I.out -e '+tmpname+'-%I.err')# run all ssaha jobs in job array
			else:
				os.system('echo \'bash ${LSB_JOBINDEX}'+tmpname+'_sbs.sh\' | bsub -q '+options.LSFQ+' -J'+tmpname+'_'+'"[1-'+str(count)+']" -o '+tmpname+'-%I.out -e '+tmpname+'-%I.err')
			
			if not options.dirty:
					os.system('bsub -w \'ended('+tmpname+'_'+')\' rm -rf *'+tmpname+'_sbs.sh '+tmpname+'*.out '+tmpname+'*.err '+options.ref+'.pac '+options.ref+'.ann '+options.ref+'.rpac '+options.ref+'.amb '+options.ref+'.rbwt '+options.ref+'.bwt '+options.ref+'.sa '+options.ref+'.rsa '+tmpname+'*.sma '+tmpname+'*.smi '+tmpname+'_unzipped '+tmpname+'_unbammed')#when job array is all done delete
		
		
		if options.pseudosequence:
			if len(bamfiles)>200:
				if count==0:
					os.system('echo '+joinstring+' | bsub -M 10000000 -q long -R \'select[mem>10000] rusage[mem=10000]\' -J'+tmpname+'_joining -o '+options.output+'_join.bsubout -e '+options.output+'_join.bsuberr')
				else:
					os.system('echo '+joinstring+' | bsub -M 10000000 -q long -R \'select[mem>10000] rusage[mem=10000]\' -J'+tmpname+'_joining -w \'ended('+tmpname+'_'+')\' -o '+options.output+'_join.bsubout -e '+options.output+'_join.bsuberr ')
				os.system('echo '+summarystring+' | bsub -M 16000000 -q long -R \'select[mem>16000] rusage[mem=16000]\' -w \'ended('+tmpname+'_joining)\' -o '+options.output+'_summary.bsubout -e '+options.output+'_summary.bsuberr')
			else:
				if count==0:
					os.system('echo '+joinstring+' | bsub -J'+tmpname+'_joining -o '+options.output+'_join.bsubout -e '+options.output+'_join.bsuberr')
				else:
					os.system('echo '+joinstring+' | bsub -J'+tmpname+'_joining -w \'ended('+tmpname+'_'+')\' -o '+options.output+'_join.bsubout -e '+options.output+'_join.bsuberr')
				os.system('echo '+summarystring+' | bsub -M 16000000 -q long -R \'select[mem>16000] rusage[mem=16000]\' -w \'ended('+tmpname+'_joining)\' -o '+options.output+'_sum.bsubout -e '+options.output+'_sum.bsuberr')
			
		#os.system('bsub -w \'ended('+tmpname+'_'+program+')\' rm *'+tmpname+'_sbs.sh; for f in tmpname*.err; do if test ! -s $f;then rm $f ${f%.err}.out;fi;done')
	
	elif options.pseudosequence:
		os.system(joinstring)
		os.system(summarystring)

		

			
			
			
			
