#!/usr/bin/env python
import string, re
import os, sys
from optparse import OptionParser, OptionGroup
import shlex, subprocess
#probably should add this to make sequences in real fasta format?
#from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord



####################
# Set some globals #
####################

SAMTOOLS_DIR="/nfs/users/nfs_s/sh16/samtools-0.1.17/"
BCFTOOLS_DIR="/nfs/users/nfs_s/sh16/samtools-0.1.17/bcftools/"
SSAHA_DIR="/nfs/users/nfs_s/sh16/ssaha2_v2.5.1_x86_64/"
BWA_DIR=""
MY_SCRIPTS_DIR="/nfs/users/nfs_s/sh16/scripts/"


#################################
# Simple Error Printing Funtion #
#################################

def DoError(ErrorString):
	print "!!!Error:", ErrorString,"!!!"
	sys.exit()


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	
	group = OptionGroup(parser, "IO Options")
	group.add_option("-b", "--bcf", action="store", dest="bcf", help="bcf file", default="", metavar="file")
	group.add_option("-v", "--vcf", action="store", dest="vcf", help="dindel vcf file [default is bcf]", default="", metavar="file")
	group.add_option("-o", "--output", action="store", dest="output", help="prefix for output files", default="", metavar="prefix")
#	group.add_option("-p", "--pseudosequence", action="store", dest="pseudosequence", help="output pseudosequence including indels [default is to write a psedoalignment to the reference and separate indel text file]", default=False)
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "Filtering Options")
	group.add_option("-d", "--depth", action="store", dest="depth", help="Minimum number of reads matching SNP [default= %default]", default=4, type="int", metavar="int")
	group.add_option("-D", "--stranddepth", action="store", dest="stranddepth", help="Minimum number of reads matching SNP per strand [default= %default]", default=2, type="int", metavar="int")
	group.add_option("-r", "--ratio", action="store", dest="ratio", help="Minimum ratio of first to second base call [default= %default]", default=0.8, type="float", metavar="float")
	#are the three above options still necessary???
	group.add_option("-q", "--QUAL", action="store", dest="QUAL", help="Minimum base quality [default= %default]", default=50.0, type="float", metavar="float")
	group.add_option("-m", "--MQUAL", action="store", dest="MQUAL", help="Minimum mapping quality [default= %default]", default=0.0, type="float", metavar="float")
	group.add_option("-a", "--AF1", action="store", dest="AF1", help="Minimum allele frequency (you would expect an AF of 1 for haploid SNPs). For non-SNP bases, the program will use 1- this number. [default= %default]", default=0.95, type="float")
	group.add_option("-c", "--CI95", action="store", dest="CI95", help="Maximum 95% confidence interval variation from AF. [default= %default]", default=0.0, type="float", metavar="float")
	group.add_option("-S", "--strandbias", action="store", dest="strand_bias", help="p-value cutoff for strand bias. [default= %default]", default=0.001, type="float", metavar="float")
	group.add_option("-Q", "--baseqbias", action="store", dest="baseq_bias", help="p-value cutoff for base quality bias. [default= %default]", default=0.0, type="float", metavar="float")
	group.add_option("-M", "--mappingbias", action="store", dest="mapping_bias", help="p-value cutoff for mapping bias. [default= %default]", default=0.001, type="float", metavar="float")
	group.add_option("-T", "--taildistancebias", action="store", dest="tail_bias", help="p-value cutoff for tail distance bias. [default= %default]", default=0.001, type="float", metavar="float")
	
	parser.add_option_group(group)
	
	
	
	
	return parser.parse_args()


################
# Main program #
################		

if __name__ == "__main__":


	#Get command line arguments

	(options, args) = main()


	#Do some checking of the input files
	
	if options.bcf=="":
		DoError("No input bcf file specified")
	elif not os.path.isfile(options.bcf):
		DoError("Cannot find input bcf file")
	if options.vcf=="":
		DoError("No input vcf file specified")
	elif not os.path.isfile(options.vcf):
		DoError("Cannot find input vcf file")
	if options.output=="":
		DoError("No output prefix specified")
		
	
	
	if options.stranddepth<0:
		print "Minimum number of reads matching SNP on each strand must be >=0. Resetting to 0"
		options.stranddepth=0
	if options.depth<(options.stranddepth*2):
		print "Minimum number of reads matching SNP must be at least double that for each strand. Resetting to", options.stranddepth*2
		options.stranddepth=options.stranddepth*2
	if options.ratio<0.5 or options.ratio>1:
		DoError("Ratio of first to second base (-r) must be greater than 0.5 and less than or equal to 1")
	if options.QUAL<0 or options.QUAL>99:
		DoError("Base quality (-q) must be between 0 and 99")
	if options.MQUAL<0 or options.MQUAL>99:
		DoError("Mapping quality (-m) must be between 0 and 99")
	if options.AF1<0 or options.AF1>1:
		DoError("Minimum allele frequency for SNPs (-a) must be between 0 and 1")
	if options.CI95<0 or options.CI95>1:
		DoError("Maximum 95% confidence interval of allele frequency (-c) must be between 0 and 1")
	if options.strand_bias<0 or options.strand_bias>1:
		DoError("p-value cutoff for strand bias (-S) must be between 0 and 1")
	if options.baseq_bias<0 or options.baseq_bias>1:
		DoError("p-value cutoff for base quality bias (-Q) must be between 0 and 1")
	if options.mapping_bias<0 or options.mapping_bias>1:
		DoError("p-value cutoff for mapping bias (-M) must be between 0 and 1")
	if options.tail_bias<0 or options.tail_bias>1:
		DoError("p-value cutoff for tail distance bias (-T) must be between 0 and 1")
	
	
	#create file of candidate indels
	
	os.system('grep -v "^#" '+options.vcf+' | awk \'{print $1 " " $2}\' > '+options.output+'tmp.indellist')
	
	
	bcftoolssarg = shlex.split(BCFTOOLS_DIR+"bcftools view -I -l "+options.output+"tmp.indellist "+options.bcf)

	returnval = subprocess.Popen(bcftoolssarg, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

	stdout, stderr  = returnval.communicate()
	
	if len(stderr)>0:
		bcftoolssarg = shlex.split(OLD_BCFTOOLS_DIR+"bcftools view -I -l "+options.output+"tmp.indellist "+options.bcf)
		returnval = subprocess.Popen(bcftoolssarg, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

		stdout, stderr  = returnval.communicate()

		if len(stderr)>0:
			print "Failed to open ", filename, "bcftools error:", stderr
			sys.exit()
	
	
	lines=stdout.split("\n")
	
	indelinfo={}
	
	infocolumn=-1
	poscolumn=-1
	for line in lines:
		words=line.strip().split()
		if len(words)==0 or len(words[0])==0:
			continue
		
		if words[0][0]=="#":
			if len(words[0])>1 and words[0][1]!="#":
				foundinfo=False
				foundpos=False
				for x, word in enumerate(words):
					if word=="INFO":
						foundinfo=True
						infocolumn=x
					elif word=="POS":
						foundpos=True
						poscolumn=x
				if not foundinfo:
					print "No info column found"
					sys.exit()
				if not foundpos:
					print "No pos column found"
					sys.exit()
				samplename=words[-1]
		elif infocolumn==-1:
			print "No info column found"
			sys.exit()
		elif poscolumn==-1:
			print "No pos column found"
			sys.exit()
		elif len(words)>infocolumn:
			info=words[infocolumn].split(";")
			
			indelinfo[words[poscolumn]]={}
			for i in info:
				if i.split("=")[0]=="DP4":
					indelinfo[words[poscolumn]]["DP4"]=i.split("=")[1].split(",")
				else:
					indelinfo[words[poscolumn]][i.split("=")[0]]=i.split("=")[1]
			
				

	infocolumn=-1
	poscolumn=-1
	indeloutputlines={}
	for line in open(options.vcf, "rU"):
		words=line.strip().split()
		if len(words)==0 or len(words[0])==0:
			continue
		elif words[0][0]=="#":
			if len(words[0])>1 and words[0][1]!="#":
				foundinfo=False
				foundpos=False
				for x, word in enumerate(words):
					if word=="INFO":
						foundinfo=True
						infocolumn=x
					elif word=="POS":
						foundpos=True
						poscolumn=x
				if not foundinfo:
					print "No info column found"
					sys.exit()
				if not foundpos:
					print "No pos column found"
					sys.exit()
					
		elif infocolumn==-1:
			print "No info column found"
			sys.exit()
		elif poscolumn==-1:
			print "No pos column found"
			sys.exit()
		elif len(words)>infocolumn:
			info=words[infocolumn].split(";")
			infodata=["INDEL"]
			
			if words[poscolumn] not in indelinfo:
				indelinfo[words[poscolumn]]={}
			if "DP4" not in indelinfo[words[poscolumn]]:
				indelinfo[words[poscolumn]]["DP4"]=["0","0","0","0"]
			if "DP" not in indelinfo[words[poscolumn]]:
				indelinfo[words[poscolumn]]["DP"]=["0"]
			for i in info:
				if i.split("=")[0]=="NF":
					indelinfo[words[poscolumn]]["DP4"][2]=i.split("=")[1]
				elif i.split("=")[0]=="NR":
					indelinfo[words[poscolumn]]["DP4"][3]=i.split("=")[1]
				elif i.split("=")[0]=="HP":	
					indelinfo[words[poscolumn]]["HP"]=i.split("=")[1]
			
			totalindelreads=int(indelinfo[words[poscolumn]]["DP4"][2])+int(indelinfo[words[poscolumn]]["DP4"][3])
			totalnonindelreads=int(indelinfo[words[poscolumn]]["DP"])-totalindelreads
			if totalnonindelreads<0:
				totalnonindelreads=0
				indelinfo[words[poscolumn]]["DP"]=str(totalindelreads)
			
			indelinfo[words[poscolumn]]["DP4"][0]=str(totalnonindelreads/2)
			indelinfo[words[poscolumn]]["DP4"][1]=str(totalnonindelreads-(totalnonindelreads/2))
			#indelinfo[words[poscolumn]]["AF1"]=str(1-int(indelinfo[words[poscolumn]]["AF1"]))
			for info in ["DP", "DP4","MQ","HP"]:#["DP", "AF1","AC1","DP4","MQ","FQ","PV4","HP"]:
				if info in indelinfo[words[poscolumn]]:
					if info=="DP4":
						infodata.append("DP4="+','.join(indelinfo[words[poscolumn]]["DP4"]))
					else:
						infodata.append(info+"="+indelinfo[words[poscolumn]][info])
					
				
			
			
			words[infocolumn]=";".join(infodata)
			if not words[0] in indeloutputlines:
				indeloutputlines[words[0]]={}
			words[6]="."
			
			if len(words[4])>len(words[3]):
				words[3]=words[3].lower()
				words[4]=words[4][:len(words[3])].lower()+words[4][len(words[3]):].upper()
			elif len(words[4])<len(words[3]):
				words[3]=words[3].lower()
				words[4]=words[4].lower()
			words[8]="."
			words[9]="."
			indeloutputlines[words[0]][words[poscolumn]]="\t".join(words)
	
	
	
	
	
	outfile=open(options.output+"tmp.vcf", "w")	
	
	print >> outfile, "##fileformat=VCFv4.1"
	print >> outfile, "##source=add_dindel_indels_to_bcf.py"
	print >> outfile, """##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Root-mean-square mapping quality of covering reads">
##INFO=<ID=FQ,Number=1,Type=Float,Description="Phred probability of all samples being the same">
##INFO=<ID=AF1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele frequency (assuming HWE)">
##INFO=<ID=AC1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele count (no HWE assumption)">
##INFO=<ID=G3,Number=3,Type=Float,Description="ML estimate of genotype frequencies">
##INFO=<ID=HWE,Number=1,Type=Float,Description="Chi^2 based HWE test P-value based on G3">
##INFO=<ID=CLR,Number=1,Type=Integer,Description="Log ratio of genotype likelihoods with and without the constraint">
##INFO=<ID=UGT,Number=1,Type=String,Description="The most probable unconstrained genotype configuration in the trio">
##INFO=<ID=CGT,Number=1,Type=String,Description="The most probable constrained genotype configuration in the trio">
##INFO=<ID=PV4,Number=4,Type=Float,Description="P-values for strand bias, baseQ bias, mapQ bias and tail distance bias">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=PC2,Number=2,Type=Integer,Description="Phred probability of the nonRef allele frequency in group1 samples being larger (,smaller) than in group2.">
##INFO=<ID=PCHI2,Number=1,Type=Float,Description="Posterior weighted chi^2 P-value for testing the association between group1 and group2 samples.">
##INFO=<ID=QCHI2,Number=1,Type=Integer,Description="Phred scaled PCHI2.">
##INFO=<ID=PR,Number=1,Type=Integer,Description="# permutations yielding a smaller PCHI2.">
##INFO=<ID=HP,Number=1,Type=Integer,Description="Reference homopolymer tract length">
##ALT=<ID=DEL,Description="Deletion">
##FILTER=<ID=q20,Description="Quality below 20">
##FILTER=<ID=hp10,Description="Reference homopolymer length was longer than 10">
##FILTER=<ID=fr0,Description="Non-ref allele is not covered by at least one read on both strands">
##FILTER=<ID=wv,Description="Other indel in window had higher likelihood">"""
	
	
	
	
	bcftoolssarg = shlex.split(BCFTOOLS_DIR+"bcftools view -I "+options.bcf)

	returnval = subprocess.Popen(bcftoolssarg, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

	stdout, stderr  = returnval.communicate()
	
	if len(stderr)>0:
		bcftoolssarg = shlex.split(OLD_BCFTOOLS_DIR+"bcftools view -I "+options.bcf)
		returnval = subprocess.Popen(bcftoolssarg, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

		stdout, stderr  = returnval.communicate()

		if len(stderr)>0:
			print "Failed to open ", filename, "bcftools error:", stderr
			sys.exit()
	
	
	lines=stdout.split("\n")
	
	
	print len(lines)
	for line in lines:
		if len(line)>1 and line[:2]=="##":
			continue
		elif len(line)>0 and line[0]=="#":
			print line.strip()
			print >> outfile, line.strip()
		elif len(line)>0:
			words=line.strip().split()
			words[8]="."
			words[9]="."
			#print >> outfile, '\t'.join(words)
			#print >> outfile, line.strip()
			#print line
			if words[0] in indeloutputlines and words[1] in indeloutputlines[words[0]]:
				print >> outfile, indeloutputlines[words[0]][words[1]]
				print indeloutputlines[words[0]][words[1]]
				del indeloutputlines[words[0]][words[1]]
			
	
	
	
	
	outfile.close()	
	os.system(BCFTOOLS_DIR+"bcftools view -S -b -D testlist "+options.output+"tmp.vcf > "+options.output+"tmp.bcf")
	os.system(BCFTOOLS_DIR+"bcftools index "+options.output+"tmp.bcf")
	
	
	
	
