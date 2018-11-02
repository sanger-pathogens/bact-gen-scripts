#!/usr/bin/env python

#Compares mummer snp output files.  NOTE: must be from the same reference sequence



##################
# Import modules #
##################

import string, re
import os, sys, random, math, time
from optparse import OptionParser, OptionGroup
from random import *

##########################
# Error message function #
##########################

def DoError(errorstring):
	print "\nError:", errorstring
	print "\nFor help use -h or --help\n"
	sys.exit()


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	
	group = OptionGroup(parser, "General options")
	group.add_option("-r", "--reference", action="store", dest="ref", help="Reference fasta", default="", metavar="FILE")
	group.add_option("-t", "--input_type", action="store", dest="input_type", help="Input type (choice of fastq or assembly) [default= %default]", default="fastq", type="choice", choices=["fastq","assembly"])
	group.add_option("-s", "--single", action="store_true", dest="single", help="Fastq files are not paired", default=False)
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

	chars = string.ascii_letters + string.digits
	
	tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))
	
	pwd=os.getcwd()
	
	assemblies=[]
	
	if options.input_type=="fastq":
		fastqprefixlist={}
		if not options.single:
			for fastqfile in args:
				if not os.path.isfile(fastqfile):
					print "Cannot find file", fastqfile, "skipping..."
					continue
				if not '.'.join(fastqfile.split(".")[:-1]).split("_")[-1] in ["1","2"]:
					print "file", fastqfile, "does not have _1 or _2 suffix... are your reads paired?"
				elif "_".join('.'.join(fastqfile.split(".")[:-1]).split("_")[:-1]) not in fastqprefixlist:
					fastqprefixlist["_".join('.'.join(fastqfile.split(".")[:-1]).split("_")[:-1])]=1
				else:
					fastqprefixlist["_".join('.'.join(fastqfile.split(".")[:-1]).split("_")[:-1])]+=1
		
		
		
			for fastqfile in fastqprefixlist:
				
				os.system("/nfs/pathogen/sh16_scripts/iterative_assembler.py -f "+fastqfile+"_1.fastq -r "+fastqfile+"_2.fastq -o "+tmpname+"_"+fastqfile+".fasta")
				assemblies.append(tmpname+"_"+fastqfile+".fasta")
				
		else:
			print "single end reads not available yet"
			sys.exit()
			for fastqfile in fastqprefixlist:
				os.system("/nfs/pathogen/sh16_scripts/iterative_assembler.py -f "+fastqfile)
	
	else:
		assemblies=args
	
	#Make pseudofastq for reference
	os.system("/nfs/pathogen/sh16_scripts/fasta2fastq_shredder.py "+options.ref+" "+tmpname+" 76 3 c 200")
	
	foldernames=[]
	for assembly in assemblies:
		if len(assembly.split('.'))>1:
			foldername='.'.join(assembly.split('.')[:-1])+"_dir"
		else:
			foldername=assembly+"_dir"
		foldernames.append(foldername)
		os.system("mkdir "+foldername)
		os.system("mv "+assembly+" "+foldername)
		os.chdir(pwd+"/"+foldername)
	
		os.system("/nfs/pathogen/sh16_scripts/multiple_mappings_to_bam.py -r "+assembly+" -M -f -p bwa "+pwd+"/"+tmpname+"_[12].fastq >  /dev/null 2>&1")
	
		os.chdir(pwd)
	
	
	
	
	
	
	sys.exit()
	
	reffile=tmpname+"_ref.mfa"
	
	
	
	folder=sys.argv[1]
	if os.path.isfile(folder+"/"+folder.replace("_bwa","")+"_unmapped.fastq.gz"):
		os.system("gunzip -f "+folder+"/"+folder.replace("_bwa","")+"_unmapped.fastq.gz")
	
	#outfile=open(tmpname+".fastq","w")
	#if os.path.isfile(folder+"/"+folder.replace("_bwa","")+"_unmapped.fastq"):
	#	lines=open(folder+"/"+folder.replace("_bwa","")+"_unmapped.fastq", "rU").readlines()
	#else:
	#	print "Cannot find "+folder+"/"+folder.replace("_bwa","")+"_unmapped.fastq"
	#	sys.exit()
	#x=0
	#while x<len(lines):
	#	newlines=[lines[x].strip()]
	#	x+=1
	#	for y in range(0,3):
	#		newlines.append(lines[x].strip())
	#		x+=1
	#	
	#	if x==len(lines) or lines[x].strip()[:-2]!=newlines[0].strip()[:-2]:
	#		continue
	#		
	#	for y in range(0,4):
	#		newlines.append(lines[x].strip())
	#		x+=1
	#	
	#	for outline in newlines:
	#		print >> outfile, outline
	#
	#outfile.close()
	
	depths=[]
	for folder in sys.argv[1:]:
		if os.path.isfile(folder+"/"+folder.replace("_bwa","")+"_Mapped.plot.gz"):
			os.system("gunzip -f "+folder+"/"+folder.replace("_bwa","")+"_Mapped.plot.gz")
			
		if os.path.isfile(folder+"/"+folder.replace("_bwa","")+"_Mapped.plot"):
			lines=open(folder+"/"+folder.replace("_bwa","")+"_Mapped.plot", "rU").readlines()
			
			exp=0.0
			for line in lines:
				exp+=int(line.strip())
			
			exp=exp/len(lines)
			depths.append([exp, folder])
			
			os.system("gzip "+folder+"/"+folder.replace("_bwa","")+"_Mapped.plot")
			
		else:
			exp=15
			depths.append([0, folder])
	
	depths.sort()
	depths.reverse()
	
	print depths[0][1], depths[0][0]
	
	folder=depths[0][1]
	exp=depths[0][0]
	if exp >40:
		exp=40
	elif exp<30:
		exp=30
	os.system('/nfs/pathogen/sh16_scripts/velvet_assembly.sh -n -e '+str(int(exp/2))+' -o "-min_contig_lgth 500" -p -f '+folder+"/"+folder.replace("_bwa","")+"_unmapped.fastq")
	os.system("gzip "+folder+"/"+folder.replace("_bwa","")+"_unmapped.fastq")
	os.system("mv "+folder+"/"+folder.replace("_bwa","")+"_unmapped_velvet "+tmpname+"_velvet")
	os.system("cp "+tmpname+"_velvet/contigs.fa "+reffile)
	
	for strains in depths[2:]:
	
		exp=strains[0]
		folder=strains[1]
		
		print folder, exp
		
		if exp >40:
			exp=40
		elif exp<30:
			exp=30
		
		if os.path.isfile(folder+"/"+folder.replace("_bwa","")+"_unmapped.fastq.gz"):
			os.system("gunzip -f "+folder+"/"+folder.replace("_bwa","")+"_unmapped.fastq.gz")
		
		outfile1=open(tmpname+"_1.fastq","w")
		outfile2=open(tmpname+"_2.fastq","w")
		
		if os.path.isfile(folder+"/"+folder.replace("_bwa","")+"_unmapped.fastq"):
			lines=open(folder+"/"+folder.replace("_bwa","")+"_unmapped.fastq", "rU").readlines()
		else:
			print "Cannot find "+folder+"/"+folder.replace("_bwa","")+"_unmapped.fastq"
			sys.exit()
		x=0
		while x<len(lines):
			newlines1=[lines[x].strip()]
			x+=1
			for y in range(0,3):
				newlines1.append(lines[x].strip())
				x+=1
			
			
			
			if x==len(lines) or lines[x].strip()[:-2]!=newlines1[0].strip()[:-2]:
				continue
			
			newlines2=[lines[x].strip()]
			x+=1
			for y in range(0,3):
				newlines2.append(lines[x].strip())
				x+=1
			
			for outline in newlines1:
				print >> outfile1, outline
			
			for outline in newlines2:
				print >> outfile2, outline
		
		outfile1.close()
		outfile2.close()
		#os.system("gunzip "+folder+"/unmapped.fastq")
		
		#newref=open(tmpname+".dna", "w")
		#print >> newref, ">Reference"
		#newref.close()
		#os.system('grep -v "^>" '+reffile+" >> "+tmpname+".dna")
		
		
		os.system("/nfs/pathogen/sh16_scripts/multiple_mappings_to_bam.py -r "+reffile+" -M -f -p bwa -L "+tmpname+"_[12].fastq >  /dev/null 2>&1")
		os.system("/nfs/pathogen/sh16_scripts/get_unmapped_fastq_from_bam.py "+tmpname+"_bwa/"+tmpname+".bam")
		#os.system("gunzip "+tmpname+"/"+tmpname+"_unmapped.fastq")
		
		
		os.system('mv '+tmpname+"_bwa/"+tmpname+"_unmapped.fastq "+tmpname+"_bwa/"+folder.replace("_bwa","")+"_unmapped.fastq")
		os.system('/nfs/pathogen/sh16_scripts/velvet_assembly.sh -n -e '+str(int(exp/2))+' -o "-min_contig_lgth 500" -p -f '+tmpname+"_bwa/"+folder.replace("_bwa","")+"_unmapped.fastq")
	#	if os.path.isfile(tmpname+"_bwa/"+tmpname+"_unmapped.fastq"):
	#		os.system("gzip "+tmpname+"_bwa/"+tmpname+"_unmapped.fastq")
		os.system("cat "+tmpname+"_bwa/"+folder.replace("_bwa","")+"_unmapped_velvet/contigs.fa >> "+reffile)
		os.system("rm -rf "+tmpname+"_bwa 1tmp*_sbs.sh")
	
	os.system("mv "+reffile+" accessory_genome.mfa")	
	os.system("rm -rf "+tmpname+"*")
	