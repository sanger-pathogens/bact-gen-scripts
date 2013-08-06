#!/usr/bin/env python

##################
# Import modules #
##################

import string, re
import os, sys, random, math, time
from optparse import OptionParser, OptionGroup
from random import *
from Bio import SeqIO
from Bio.Seq import Seq
from shutil import copyfile
import mimetypes
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
import farm

SMALT_LOC="/software/pathogen/external/apps/usr/bin/smalt"
SAMTOOLS_LOC="/software/pathogen/external/apps/usr/bin/samtools"
AGA_DIR="/nfs/users/nfs_s/sh16/scripts/Aga/"


##########################
# Error message function #
##########################

def DoError(errorstring):
	print "\nError:", errorstring
	print "\nFor help use -h or --help\n"
	sys.exit()


###########################################
# Function to report errors in bash files #
###########################################

def add_bash_error_function(file_handle):
	
	print >> file_handle, "function error_exit"
	print >> file_handle, "{"
	print >> file_handle, '	echo "$1" 1>&2'
	print >> file_handle, "	exit 1"
	print >> file_handle, "}"
	
	return


####################################################
# Define a class to contain fastq file information #
####################################################

class fastq:
	def __init__(self, file_name):
			
		self.file_name=file_name
		self.is_gzipped=False
		self.is_read1=False
		self.is_read2=False
		self.is_paired=False
		self.name=""
		self.base_name=""
		self.read_ok=False
		self.absolute_path=""
		self.extensions=[]
		self.run=False
		self.lane=False
		self.tag=False
		self.sanger_name=False
			
		if not os.path.isfile(self.file_name):
			print "Cannot find file", self.file_name
			print "File excluded"
			return
		
		self.absolute_path=os.path.abspath(self.file_name)
		self.name=os.path.basename(self.file_name)
		self.extensions=self.name.split(".")[1:]
		self.base_name=self.name.split(".")[0]
		
		if not "fastq" in self.extensions:
			print self.file_name, "does not appear to be a fastq file (does not include the extension fastq)"
			print "File excluded"
			return
		
		if mimetypes.guess_type(self.file_name)[1]=="gzip":
			self.is_gzipped=True
		
		
		if re.search(r"^\d+_\d+#\d+_[12]$", self.base_name) is not None:
			self.sanger_name=True
			self.is_paired=True
			split_name=self.base_name.replace("#","_").split("_")
			self.run=split_name[0]
			self.lane=split_name[1]
			self.tag=split_name[2]
			if split_name[3]=="1":
				self.is_read1=True
			elif split_name[3]=="2":
				self.is_read2=True
			self.base_name="_".join(self.base_name.split("_")[:-1])
		
		elif re.search(r"^\d+_\d+#\d+$", self.base_name) is not None:
			self.sanger_name=True
			self.is_paired=False
			split_name=self.base_name.replace("#","_").split("_")
			self.run=split_name[0]
			self.lane=split_name[1]
			self.tag=split_name[2]
			if split_name[3]=="1":
				self.is_read1=True
			else:
				self.is_read2=True
		
		elif re.search(r"^\d+_\d+_[12]$", self.base_name) is not None:
			self.sanger_name=True
			self.is_paired=True
			split_name=self.base_name.replace("#","_").split("_")
			self.run=split_name[0]
			self.lane=split_name[1]
			if split_name[2]=="1":
				self.is_read1=True
			elif split_name[3]=="2":
				self.is_read2=True
			self.base_name="_".join(self.base_name.split("_")[:-1])
		
		elif re.search(r"^\d+_\d+$", self.base_name) is not None:
			self.sanger_name=True
			self.is_paired=False
			split_name=self.base_name.replace("#","_").split("_")
			self.run=split_name[0]
			self.lane=split_name[1]
		
		else:
			split_name=self.base_name.split("_")
			
			if len(split_name)>1:
				if split_name[-1]=="1":
					self.is_read1=True
					self.is_paired=True
					self.base_name="_".join(self.base_name.split("_")[:-1])
				elif split_name[-1]=="2":
					self.is_read2=True
					self.is_paired=True
					self.base_name="_".join(self.base_name.split("_")[:-1])

		self.read_ok=True
		return
	
	def describe(self):
		print
		print "Information for fastq file", self.name+":"
		print
		print "\tAbsolute path:", self.absolute_path
		print "\tBase name:", self.base_name
		if self.sanger_name:
			print "\tRun:", self.run
			print "\tLane:", self.lane
			if self.tag:
				print "\tTag:", self.tag
		if self.is_paired:
			if self.is_read1:
				print "\tRead is first of a pair"
			elif self.is_read2:
				print "\tRead is second of a pair"
		else:
			print "\tRead is single ended"
		if self.is_gzipped:
			print "\tFile is gzipped"
		
		print
		
		



##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	
	group = OptionGroup(parser, "General options")
	group.add_option("-r", "--reference", action="store", dest="ref", help="Reference fasta", default="", metavar="FILE")
	group.add_option("-t", "--tab", action="store", dest="tab", help="Reference tab file of MGEs", default="", metavar="FILE")
	group.add_option("-o", "--output_prefix", action="store", dest="prefix", help="Output prefix", default="")
	group.add_option("-i", "--maxinsert", action="store", dest="maxinsertsize", help="maximum insert size [default= %default]", default=1000, type="int", metavar="INT")
	group.add_option("-j", "--mininsert", action="store", dest="mininsertsize", help="minimum insert size [default= %default]", default=50, type="int", metavar="INT")
	group.add_option("-z", "--nomapid", action="store", dest="nomapid", help="Minimum identity threshold to report a mapping Specified as a positive integer or proportion of read length [default= %default]", default=0, type="float", metavar="float")
	parser.add_option_group(group)
#	group.add_option("-d", "--contaminant_database", action="store", dest="contaminants", help="Name file containing contaminant accession numbers", default=False, metavar="FILE")
#	group.add_option("-H", "--human", action="store_true", dest="human", help="Blast primers against human genome", default=False)
#	parser.add_option_group(group)
#if len(sys.argv)!=3:
#	print "Usage: create_pan_genome.py <ssaha_folders>"
#	sys.exit()


	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


#	if options.processors>cpu_count():
#		DoError('You do not have '+str(options.processors)+' processors')
#
	if options.ref=='':
		DoError('No reference file selected')
	elif not os.path.isfile(options.ref):
		DoError('Cannot find file '+options.ref)
	
	
	return
	

################
# Main program #
################		

if __name__ == "__main__":
	
	starttime=time.clock()

	#Get command line arguments

	(options, args) = main()
	
	#Do some checking of the input files
	
	check_input_validity(options, args)
	
	#Create a temporary name
	
	chars = string.ascii_letters + string.digits
	tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))
	tmpname='Sitmp'
	
	#Create a temporary directory to run the analyses in
	
	if not os.path.exists(tmpname):
		os.makedirs(tmpname)
	
	#Check the arguments are all real fastq files and parse them in a sensible way
	
	fastqs={}
	for arg in args:
		arg_fastq=fastq(arg)
		if arg_fastq.read_ok:
			if not arg_fastq.base_name in fastqs:
				fastqs[arg_fastq.base_name]={"read1": False, "read2": False}
			if arg_fastq.is_read1:
				fastqs[arg_fastq.base_name]["read1"]=arg_fastq
			elif arg_fastq.is_read2:
				fastqs[arg_fastq.base_name]["read2"]=arg_fastq
			elif not arg_fastq.is_paired:
				fastqs[arg_fastq.base_name]["read1"]=arg_fastq
			else:
				print "Cannot understand file name of", arg+". Skipping..."
		else:
			print "Cannot understand file name of", arg+". Skipping..."
		
		
	#Create a bsub job to index the reference file ready for mapping
	
	smalt_index_command=SMALT_LOC+" index -k 13 -s 1 "+tmpname+"/"+tmpname+".index "+options.ref
	
	job1 = farm.Bsub(tmpname+"/SMALT_index_bsub.out", tmpname+"/SMALT_index_bsub.err", tmpname+"_smalt_index", "normal", 1, smalt_index_command)
	job1_id = job1.run()
	
	
	#for each fastq file (pair), create a bash script to run SMALT, extract unmapping reads and assemble them
	
	assembly_list=[]
	
	for i, fastq in enumerate(fastqs):
		if options.nomapid:
			smalt_map_command=' '.join(map(str,[SMALT_LOC+" map -y "+str(options.nomapid)+" -x -r 0 -i", options.maxinsertsize, " -j", options.mininsertsize, " -f "+smaltoutput+" -o "+tmpname+"/"+tmpname+"_"+str(i+1)+".bam", tmpname+"/"+tmpname+".index"]))
		else:
			smalt_map_command=' '.join(map(str,[SMALT_LOC+" map -x -r 0 -i", options.maxinsertsize, " -j", options.mininsertsize, " -f bam -o "+tmpname+"/"+tmpname+"_"+str(i+1)+".bam", tmpname+"/"+tmpname+".index"]))
		
		smalt_run_file=open(tmpname+"/"+tmpname+"_smalt_"+str(i+1)+".sh", "w")
		add_bash_error_function(smalt_run_file)
		
		name=""
		
		if fastqs[fastq]["read1"] and fastqs[fastq]["read2"]:
#			print >> smalt_run_file, smalt_map_command, fastqs[fastq]["read1"].absolute_path, fastqs[fastq]["read2"].absolute_path, '|| error_exit "SMALT command failed! Aborting"'
			name=fastqs[fastq]["read1"].base_name
		elif fastqs[fastq]["read1"]:
#			print >> smalt_run_file, smalt_map_command, fastqs[fastq]["read1"].absolute_path, '|| error_exit "SMALT command failed! Aborting"'
			name=fastqs[fastq]["read1"].base_name
		elif fastqs[fastq]["read2"]:
#			print >> smalt_run_file, smalt_map_command, fastqs[fastq]["read2"].absolute_path, '|| error_exit "SMALT command failed! Aborting"'
			name=fastqs[fastq]["read2"].base_name
#		
#		print >> smalt_run_file, SAMTOOLS_LOC+" sort "+tmpname+"/"+tmpname+"_"+str(i+1)+".bam "+tmpname+"/"+name+' || error_exit "samtools sort command failed! Aborting"'
#		print >> smalt_run_file, SAMTOOLS_LOC+" index "+tmpname+"/"+name+".bam"+' || error_exit "samtools index command failed! Aborting"'
#		print >> smalt_run_file, "rm -f "+tmpname+"/"+tmpname+"_"+str(i+1)+'.bam || error_exit "rm command failed! Aborting"'
		
		if fastqs[fastq]["read1"] and fastqs[fastq]["read2"]:
			print >> smalt_run_file, AGA_DIR+"bam_filter.py -b "+tmpname+"/"+name+".bam -o "+tmpname+"/"+name+"_unmapped -f pairedfastq -t atleastoneunmapped"+' || error_exit "bam_filter command failed! Aborting"'
			print >> smalt_run_file, AGA_DIR+"velvet_assembly.sh -f "+tmpname+"/"+name+"_unmapped_1.fastq -r "+tmpname+"/"+name+"_unmapped_2.fastq -s "+tmpname+"/"+name+"_shuffled.fastq -n -p"+' || error_exit "velvet assembly command failed! Aborting"'
		else:
			print >> smalt_run_file, AGA_DIR+"bam_filter.py -b "+tmpname+"/"+name+".bam -o "+tmpname+"/"+name+"_unmapped -f fastq -t atleastoneunmapped"+' || error_exit "bam_filter command failed! Aborting"'
			print >> smalt_run_file, AGA_DIR+"velvet_assembly.sh -f "+tmpname+"/"+name+"_unmapped.fastq -s "+tmpname+"/"+name+"_shuffled.fastq -n"+' || error_exit "velvet assembly command failed! Aborting"'
		
		print >> smalt_run_file, "mv "+tmpname+"/"+name+"_shuffled_velvet/contigs.fa "+tmpname+"/"+name+"_assembled.fasta"+' || error_exit "mv command failed! Aborting"'
		print >> smalt_run_file, "rm -rf "+tmpname+"/"+name+"_shuffled_velvet "+tmpname+"/"+name+"_shuffled.fastq "+tmpname+"/"+name+'_shuffled_velvet.log || error_exit "rm command failed! Aborting"'
		assembly_list.append(tmpname+"/"+name+"_assembled.fasta")
		
		
		smalt_run_file.close()
	
	#create a bsub array to run the smalt scripts
	
	smalt_run_file_command="bash "+tmpname+"/"+tmpname+"_smalt_INDEX.sh"
	
	job2 = farm.Bsub(tmpname+"/SMALT_bsub.out", tmpname+"/SMALT_bsub.err", tmpname+"_smalt", "normal", 2, smalt_run_file_command, start=1, end=i+1)
	#add a dependency so the mappings only run once the indexing is complete
	job2.add_dependency(job1_id) 
	job2_id = job2.run()
	
	
	
	#Pull out the core and accessory regions of the reference genome using the tab file (if provided)
	
	if options.prefix=="":
		options.prefix="Reference"
	
	
	acc_file=options.prefix+"_acc.mfa"
	core_file=options.prefix+"_core.mfa"
	
	
	if options.tab!="":
		try:
			tablines=open(options.tab,"rU").readlines()
		except StandardError:
			print "Failed to read tab file"
			sys.exit()
		ftloc=-1
		note=""
		regions=[]
		for line in tablines:
#			         111111111122222222223333333333
#			123456789012345678901234567890123456789
#			FT   misc_feature    2825393..2826315
			line=line.strip()
			if len(line)>22 and line[:2]=="FT" and len(line[2:21].strip())>0:
				try: ftloc=map(int,(line[21:].split("..")))
				except StandardError:
					print line
					sys.exit()
				note=""
			elif len(line)>22 and line[:2]=="FT" and len(line[2:21].strip())==0 and ftloc!=-1:
				if line[21:].split("=")[0]=="/note":
					note='='.join(line[22:].split("=")[1:])
					if len(ftloc)==1:
						regions.append([ftloc[0], ftloc[0], note])
					else:
						regions.append([ftloc[0], ftloc[1], note])
						
					ftloc=-1

		
		refseq=''.join(open(options.ref, "rU").read().split('\n')[1:])
		
		cur=0
		count=1
		core_count=1
		acc_count=1
		core_out=open(core_file, "w")
		acc_out=open(acc_file, "w")
		
		for region in regions:
			if region[0]-cur>1:
				print "ref", cur, "..",region[0]-1
				print >> core_out, ">Ref_"+str(count), cur, "..",region[0]-1
				print >> core_out, refseq[cur+1:region[0]]
				core_count+=1
				count+=1
				
			print "acc", region[0], "..", region[1], region[2]
			print >> acc_out, ">Ref_"+str(count), region[0], "..",region[1]-1, region[2]
			print >> acc_out, refseq[region[0]:region[1]]
			cur=region[1]
			acc_count+=1
			count+=1
		if cur<len(refseq):
			print >> core_out, ">Ref_"+str(count), cur, "..",len(refseq)
			print >> core_out, refseq[cur:]
		core_out.close()
		acc_out.close()
		
	else:
		copyfile(options.ref, core_file)
		open(acc_file, 'w').close()
			
	
	#Create a bsub job to join all of the assemblies and prepare everything for blasting against the reference core regions
	
	cleanup1_file=open(tmpname+"/"+tmpname+"_cleanup1.sh", "w")
	
	add_bash_error_function(cleanup1_file)
	
	print >> cleanup1_file, "cat "+' '.join(assembly_list)+' '+acc_file+' > '+tmpname+'/'+tmpname+'_all_seqs.fasta || error_exit "cat command failed! Aborting"'
	print >> cleanup1_file, AGA_DIR+"extract_contig_lengths.py "+tmpname+'/'+tmpname+'_all_seqs.fasta > '+tmpname+'/'+tmpname+'_seq_lengths.txt || error_exit "extracting of contig lengths failed! Aborting"'
	print >> cleanup1_file, "formatdb -p F -i "+core_file+' || error_exit "reference formatdb failed! Aborting"'

	cleanup1_file.close()
	
	#Run the bsub command once the SMALT mappingsa re all complete
	
	cleanup1_run_command="bash "+tmpname+"/"+tmpname+"_cleanup1.sh"
	
	job3 = farm.Bsub(tmpname+"/cleanup1_bsub.out", tmpname+"/cleanup1_bsub.err", tmpname+"_cleanup1", "normal", 1, cleanup1_run_command)
	job3.add_dependency(job2_id)
	job3_id = job3.run()
	
	
	#For each assembly, create a bash script to blast the assembly against the reference core and remove contigs that match (should this script also remove any overlapping contig ends that match the core?)
	
	for i, assembly in enumerate(assembly_list):
		core_blast_file=open(tmpname+"/"+tmpname+"_core_blast_"+str(i+1)+".sh", "w")
	
		add_bash_error_function(core_blast_file)
		
		print >> core_blast_file, "blastall -p blastn -i "+assembly+" -m 8 -e 1e-10 -o "+assembly+".core.blast -d "+core_file+' -F F || error_exit "core blast failed! Aborting"'
		print >> core_blast_file, AGA_DIR+"filter_ref_contig_blast.py "+assembly+'core.blast '+tmpname+'/'+tmpname+'_seq_lengths.txt '+assembly+' > '+tmpname+'/'+assembly+'.noncore.fasta || error_exit "filtering of reference blast hits failed! Aborting"'
	
		core_blast_file.close()
	
	core_blast_run_command="bash "+tmpname+"/"+tmpname+"_core_blast_INDEX.sh"
	
	job4 = farm.Bsub(tmpname+"/core_blast_bsub.out", tmpname+"/core_blast_bsub.err", tmpname+"_core_blast", "normal", 1, core_blast_run_command, start=1, end=i+1)
	job4.add_dependency(job3_id) 
	job4_id = job4.run()
