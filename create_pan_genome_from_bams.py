#!/usr/bin/env python

#Compares mummer snp output files.  NOTE: must be from the same reference sequence

##################
# Import modules #
##################

import string, re
import os, sys, random, math, time
from optparse import OptionParser, OptionGroup
from random import *
from Bio import SeqIO
from Bio.Seq import Seq


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
	group.add_option("-e", "--embl", action="store", dest="embl", help="Reference embl fasta", default="", metavar="FILE")
	group.add_option("-o", "--output_prefix", action="store", dest="prefix", help="Output prefix", default="")
	group.add_option("-m", "--no_map", action="store_false", dest="mapping", help="Do not map data against reference using bwa first [default = do mapping]", default=True)
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
	if options.embl!='' and not os.path.isfile(options.embl):
		DoError('Cannot find file '+options.embl)
	if options.prefix=='':
		DoError('No output prefix specified')
	
	
	return
	
##############################
# Function to filter repeats #
##############################
def filter_repeats(contigfile, tmpname):
	ide = 97;
	print "Filtering contigs contained in larger contigs with at least "+str(ide)+"% match"
	
	seq_records=[]
	lengths={}
	
	for seq_record in SeqIO.parse(open(contigfile), "fasta"):
		seq_records.append(seq_record)
		lengths[seq_record.id]=len(seq_record.seq)

	if len(seq_records)==0:
		DoError("Amos has made an empty fasta file...bah")
		return
	
	os.system("nucmer --maxmatch --minmatch 200 --mincluster 200 --nosimplify --prefix="+tmpname+" "+contigfile+" "+contigfile)#+" >  /dev/null 2>&1")
	os.system("show-coords -r -T -o "+tmpname+".delta > "+tmpname+".coords")
	coords_lines=open(tmpname+".coords", "rU").readlines()[4:]

	matchlengths={}
	toremove=[]
	for line in coords_lines:
		words=line.strip().split()
	   
		if words[8] not in toremove and words[7]!=words[8] and lengths[words[8]]<lengths[words[7]]:

			if not words[7] in matchlengths:
				matchlengths[words[7]]={}
			if not words[8] in matchlengths[words[7]]:
				matchlengths[words[7]][words[8]]=[0.0,0.0]

			matchlengths[words[7]][words[8]][0]+=int(words[5])
			matchlengths[words[7]][words[8]][1]+=(float(words[6])/100)*int(words[5])
		  
			
			if matchlengths[words[7]][words[8]][1]>(matchlengths[words[7]][words[8]][0]/100)*ide and matchlengths[words[7]][words[8]][0]>(float(lengths[words[8]])/100)*ide:
				toremove.append(words[8])
	   
			 

	
	
	count=1
	my_new_records=[]
	for record in seq_records:
		if record.id not in toremove:
			#record.id="contig_"+str(count)
			my_new_records.append(record)
			count+=1
	SeqIO.write(my_new_records, open(contigfile,"w"), "fasta")
	

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
	
	reffile=tmpname+"_ref.mfa"
	
#	folder=sys.argv[1]
#	if os.path.isfile(folder+"_bwa/"+folder+"_unmapped.fastq.gz"):
#		os.system("gunzip -f "+folder+"_bwa/"+folder+"_unmapped.fastq.gz")

	if options.mapping:
		os.system("~sh16/scripts/multiple_mappings_to_bam_test.py -U long -f -p bwa -X -g -r "+options.ref+" "+' '.join(args)+' > '+tmpname+'jobstring')
		
		jobnum=open(tmpname+'jobstring', "rU").readlines()[-3].split(">")[0].split("<")[1]
		
		
		print "Mapping reads against reference"
		
		todo=1
		print "JOBID    ARRAY_SPEC  OWNER  NJOBS PEND DONE  RUN EXIT SSUSP USUSP PSUSP"
		#while os.path.isfile(tmpname+'waitfile'):
		while todo!=0:
			time.sleep(10)
			bjoblines=os.popen("bjobs -A "+jobnum).readlines()
			if len(bjoblines)<2:
				continue
			bjobsstring=bjoblines[1]
			print bjobsstring.strip()+"\r",
			sys.stdout.flush()
			pend=int(bjobsstring.split()[4])
			done=int(bjobsstring.split()[5])
			run=int(bjobsstring.split()[6])
			exited=int(bjobsstring.split()[7])
			todo=run+pend
	
	
	
	#sys.exit()
	
	#outfile=open(tmpname+".fastq","w")
	#if os.path.isfile(folder+"_bwa/"+folder+"_unmapped.fastq"):
	#	lines=open(folder+"_bwa/"+folder+"_unmapped.fastq", "rU").readlines()
	#else:
	#	print "Cannot find "+folder+"_bwa/"+folder+"_unmapped.fastq"
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
	donefolders=[]
	for arg in args:
		folder="_".join('.'.join(arg.split("/")[-1].split(".")[:-1]).split("_")[:-1])
		if folder in donefolders:
			continue
		donefolders.append(folder)
		#if not os.path.isfile(folder+"_bwa/"+folder+"_unmapped.fastq") and not os.path.isfile(folder+"_bwa/"+folder+"_unmapped.fastq.gz"):
		os.system("~sh16/scripts/get_unmapped_fastq_from_bam.py "+folder+"_bwa/"+folder+".bam")
		if os.path.isfile(folder+"_bwa/"+folder+"_Mapped.plot.gz"):
			os.system("gunzip -f "+folder+"_bwa/"+folder+"_Mapped.plot.gz")
			
		if os.path.isfile(folder+"_bwa/"+folder+"_Mapped.plot"):
			lines=open(folder+"_bwa/"+folder+"_Mapped.plot", "rU").readlines()
			
			exp=0.0
			for line in lines:
				exp+=int(line.strip())
			
			exp=exp/len(lines)
			depths.append([exp, folder])
			
			os.system("gzip "+folder+"_bwa/"+folder+"_Mapped.plot")
			
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
	if os.path.isfile(folder+"_bwa/"+folder+"_unmapped.fastq.gz"):
		os.system("gunzip -f "+folder+"_bwa/"+folder+"_unmapped.fastq.gz")
	if os.path.isfile(folder+"_bwa/"+folder+"_unmapped.fastq"):
		os.system('~sh16/scripts/velvet_assembly.sh -n -e '+str(int(float(exp)/2))+' -o "-min_contig_lgth 500" -p -f '+folder+"_bwa/"+folder+"_unmapped.fastq")
	else:
		print "Cannot find "+folder+"_bwa/"+folder+"_unmapped.fastq"
		sys.exit()
	os.system("gzip -f "+folder+"_bwa/"+folder+"_unmapped.fastq")
	os.system("mv "+folder+"_bwa/"+folder+"_unmapped_velvet "+tmpname+"_velvet")
	os.system("cp "+tmpname+"_velvet/contigs.fa "+reffile)
	
	for strains in depths[1:]:
	
		exp=strains[0]
		folder=strains[1]
		
		print folder, exp
		
		if exp >40:
			exp=40
		elif exp<30:
			exp=30
		
		
		
		outfile1=open(tmpname+"_1.fastq","w")
		outfile2=open(tmpname+"_2.fastq","w")
		
		if os.path.isfile(folder+"/"+folder+"_bwa/"+folder+"_unmapped.fastq.gz"):
			os.system("gunzip -f "+folder+"_bwa/"+folder+"_unmapped.fastq.gz")
		
		if os.path.isfile(folder+"_bwa/"+folder+"_unmapped.fastq"):
			lines=open(folder+"_bwa/"+folder+"_unmapped.fastq", "rU").readlines()
		else:
			print "Cannot find "+folder+"_bwa/"+folder+"_unmapped.fastq"
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
		
		
		os.system("~sh16/scripts/multiple_mappings_to_bam_test.py -r "+reffile+" -y -U long -f -X -p bwa -L "+tmpname+"_[12].fastq")# >  /dev/null 2>&1")
		os.system("~sh16/scripts/get_unmapped_fastq_from_bam.py "+tmpname+"_bwa/"+tmpname+".bam")
		sys.exit()
		#os.system("gunzip "+tmpname+"/"+tmpname+"_unmapped.fastq")
		
		
		os.system('mv '+tmpname+"_bwa/"+tmpname+"_unmapped.fastq "+tmpname+"_bwa/"+folder+"_unmapped.fastq")
		os.system('~sh16/scripts/velvet_assembly.sh -n -e '+str(int(float(exp)/2))+' -o "-min_contig_lgth 500" -p -f '+tmpname+"_bwa/"+folder+"_unmapped.fastq")
	#	if os.path.isfile(tmpname+"_bwa/"+tmpname+"_unmapped.fastq"):
	#		os.system("gzip "+tmpname+"_bwa/"+tmpname+"_unmapped.fastq")
		os.system("cat "+tmpname+"_bwa/"+folder+"_unmapped_velvet/contigs.fa >> "+reffile)
		os.system("rm -rf "+tmpname+"_bwa 1tmp*_sbs.sh")
	
	filter_repeats(reffile, tmpname)
	
	os.system("mv "+reffile+" "+options.prefix+"_accessory_genome.mfa")	
	os.system("rm -rf "+tmpname+"*")
	os.system("cat "+options.ref+" "+options.prefix+"_accessory_genome.mfa > "+options.prefix+"_pan_genome.mfa")
	
	#Make pseudofastq for reference
	if not os.path.isfile('.'.join(options.ref.split('.')[:-1])+"_1.fastq") and not os.path.isfile('.'.join(options.ref.split('.')[:-1])+"_2.fastq"):
		refname='.'.join(options.ref.split('.')[:-1])
	else:
		refname=tmpname
	os.system("~sh16/scripts/fasta2fastq_shredder.py "+options.ref+" "+refname+" 76 3 c 200")
		
	#final mapping of all isolates against ref+accessory
	if options.embl!="":
		os.system("~sh16/scripts/multiple_mappings_to_bam.py -M -p ssaha -f -a -x -P -t -g -r "+options.prefix+"_pan_genome.mfa -e "+options.embl+" "+' '.join(args)+' '+refname+'_[12].fastq')
	else:
		os.system("~sh16/scripts/multiple_mappings_to_bam.py -M -p ssaha -f -a -x -P -t -g -r "+options.prefix+"_pan_genome.mfa "+' '.join(args)+' '+refname+'_[12].fastq')
	
	
	
	
	
	