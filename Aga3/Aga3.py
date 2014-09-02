#!/usr/bin/env python

##################
# Import modules #
##################

import string, re
import os, sys, math, time
from optparse import OptionParser, OptionGroup
from random import *
from Bio import SeqIO
from Bio.Seq import Seq
from shutil import copyfile
import mimetypes
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
import farm
import subprocess
import gzip

SMALT_LOC="/software/pathogen/external/apps/usr/bin/smalt"
SAMTOOLS_LOC="/software/pathogen/external/apps/usr/bin/samtools"
AGA_DIR="/nfs/users/nfs_s/sh16/scripts/Aga3/"


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


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options] <list of fasta assemblies>"
	parser = OptionParser(usage=usage)
	
	parser.add_option("-r", "--reference", action="store", dest="ref", help="Reference sequence file. NOTE: Currently this must be a single DNA sequence, not a multifasta.", default="", metavar="FILE")
	parser.add_option("-t", "--tab", action="store", dest="tab", help="Reference tab file of MGEs", default="", metavar="FILE")
	parser.add_option("-o", "--output_prefix", action="store", dest="prefix", help="Output prefix", default="")
	parser.add_option("-T", "--tree", action="store", dest="tree", help="Tree file to allow ordering by clade (only works if taxon names are the same in the assemblies and tree)", default="")
	parser.add_option("-M", "--midpoint", action="store_true", dest="midpoint", help="Midpoint root tree[default= %default]", default=False)
	parser.add_option("-L", "--ladderise", action="store", choices=['right', 'left'], dest="ladderise", help="ladderise tree (choose from right or left) [default= %default]", type="choice", default=None)
	

	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


	if options.ref=='':
		DoError('No reference file selected')
	elif not os.path.isfile(options.ref):
		DoError('Cannot find file '+options.ref)
		
	if len(args)==0:
		DoError("No assemblies specified")
	
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
	
	
	
	#Pull out the core and accessory regions of the reference genome using the tab file (if provided)
	
	if options.prefix=="":
		options.prefix="Reference"
	
	
	acc_file=options.prefix+"_acc.mfa"
	core_file=options.prefix+"_core.mfa"
	
	
	if options.tab!="":
		print "Extracting core and accessory regions for reference"
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
				try: ftloc=map(int,(line[21:].replace("complement(","").replace(")","").split("..")))
				except StandardError:
					print "Failed to read tab file"
					print line[21:].split("..")
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
				print >> core_out, ">Ref_"+str(count)+"#"+str(cur)+".."+str(region[0]-1)
				print >> core_out, refseq[cur+1:region[0]]
				core_count+=1
				count+=1
				
			print "acc", region[0], "..", region[1], region[2]
			print >> acc_out, ">Ref_"+str(count)+"#"+str(region[0])+".."+str(region[1]-1), region[2]
			print >> acc_out, refseq[region[0]:region[1]]
			cur=region[1]
			acc_count+=1
			count+=1
		if cur<len(refseq):
			print "ref", cur, "..",len(refseq)
			print >> core_out, ">Ref_"+str(count)+"#"+str(cur)+".."+str(len(refseq))
			print >> core_out, refseq[cur:]
		core_out.close()
		acc_out.close()
		
	else:
		copyfile(options.ref, core_file)
		open(acc_file, 'w').close()
	
	
	assembly_list=[]
	#Create a bsub job to map all of the assemblies against the reference core and pull out unmatched regions
	print "Running jobs to map assemblies to the reference core using mummer"
	for i, assembly in enumerate(args):
		name=assembly.split("/")[-1].split(".")[0]
		mummer_file=open(tmpname+"/"+tmpname+"_mummer_"+str(i+1)+".sh", "w")
		print >> mummer_file, AGA_DIR+"align_contigs_with_mummer.py -r "+core_file+" -q "+assembly+" -o "+tmpname+"/"+name+' || error_exit "mummer script failed! Aborting"'
		assembly_list.append(tmpname+"/"+name+"_novel.fasta")

		mummer_file.close()
	
	#Run the bsub command
	
	mummer_run_file_command="bash "+tmpname+"/"+tmpname+"_mummer_INDEX.sh"
	
	job1 = farm.Bsub(tmpname+"/mummer_bsub.out", tmpname+"/mummer_bsub.err", tmpname+"_mummer", "normal", 1, mummer_run_file_command, start=1, end=i+1)
	job1_id = job1.run()
	
	
	#Create a bsub job to join all of the noncore contigs and prepare everything for blasting against themselves
	print "Running job to combine unmatched regions"
	cluster_noncore=open(tmpname+"/"+tmpname+"_cluster_noncore.sh", "w")
	
	add_bash_error_function(cluster_noncore)
	assembly_list.append(acc_file)
	
	treebit=""
	if options.tree!="":
		treebit=treebit+"-t "+options.tree
	if options.midpoint:
		treebit=treebit+" -m"
	if options.ladderise!=None:
		treebit=treebit+" -l "+options.ladderise
	

	print >> cluster_noncore, AGA_DIR+"filter_contigs_with_mummer_sequential.py -o "+tmpname+'/'+tmpname+" "+acc_file+" "+tmpname+' || error_exit "Filtering noncore with mummer failed! Aborting"'
	print >> cluster_noncore, AGA_DIR+"reorder_contigs.py "+treebit+" -o "+options.prefix+"_accessory.fasta -c "+tmpname+'/'+tmpname+'_filtered.fasta || error_exit "Sorting filtered noncore failed! Aborting"'

	cluster_noncore.close()
	
	#Run the bsub command once the SMALT mappings are all complete
	
	cluster_noncore_run_command="bash "+tmpname+"/"+tmpname+"_cluster_noncore.sh"
	
	job2 = farm.Bsub(tmpname+"/cluster_noncore_bsub.out", tmpname+"/cluster_noncore_bsub.err", tmpname+"_cluster_noncore", "normal", 10, cluster_noncore_run_command)
	job2.add_dependency(job1_id)
	job2_id = job2.run()
	
	
	#Create a bsub job to join all of the noncore contigs and prepare everything for blasting against themselves
	
	create_pan_genome=open(tmpname+"/"+tmpname+"_create_pan_genome.sh", "w")
	
	add_bash_error_function(create_pan_genome)
	
	print >> create_pan_genome, "cat "+core_file+' '+options.prefix+"_accessory.fasta > "+options.prefix+'_pan.fasta || error_exit "cat command failed! Aborting"'
	print >> create_pan_genome, AGA_DIR+"better_blast.py -q "+options.prefix+'_pan.fasta -s '+options.prefix+'_pan.fasta -o '+options.prefix+"_self_blast"

	create_pan_genome.close()
	
	#Run the bsub command to create the final pan and accessory genome files
	
	create_pan_genome_command="bash "+tmpname+"/"+tmpname+"_create_pan_genome.sh"
	print "Running job to create accessory and pan genomes"
	job3 = farm.Bsub(tmpname+"/create_pan_genome_bsub.out", tmpname+"/create_pan_genome_bsub.err", tmpname+"_create_pan_genome", "normal", 1, create_pan_genome_command)
	job3.add_dependency(job2_id)
	job3_id = job3.run()
	
	
	#Run the bsub command to map assemblies back against the pan genome to create tab files
	
	print "Running jobs to map assemblies to the pan genome core using mummer"
	for i, assembly in enumerate(args):
		name=assembly.split("/")[-1].split(".")[0]
		mummer_file=open(tmpname+"/"+tmpname+"_mummerb_"+str(i+1)+".sh", "w")
		print >> mummer_file, AGA_DIR+"align_contigs_with_mummer.py -t -s -q "+options.prefix+"_pan.fasta -r "+assembly+" -o "+name+'_pan_mapping || error_exit "mummer script failed! Aborting"'
		

		mummer_file.close()
	
	#Run the bsub command
	
	mummer_run_file_command="bash "+tmpname+"/"+tmpname+"_mummerb_INDEX.sh"
	
	job4 = farm.Bsub(tmpname+"/mummerb_bsub.out", tmpname+"/mummerb_bsub.err", tmpname+"_mummerb", "normal", 1, mummer_run_file_command, start=1, end=i+1)
	job4.add_dependency(job3_id)
	job4_id = job4.run()
	

	
	
	

#
#Make a plot of the accessory genome presence/absence
#~/scripts/iCANDY.py -n -d area -y 100 *.plot Accessory_genome.fasta
#
#
#identify each gene in accessory with prodigal
#prodigal -a accessory_prodigal.faa -d accessory_prodigal.fasta -o accessory_prodigal.tab < accessory_oneseq.fasta
#run these through orthomcl, but with e=1e-50 and inflation of 2
#extract each cluster into a multifasta
#align each of these
#use the alignments as input to psi-blast vs nr database
	
	
	
