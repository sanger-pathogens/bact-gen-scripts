#!/usr/bin/env python
import string, re
import os, sys
from optparse import OptionParser, OptionGroup
import pylab
#probably should add this to make sequences in real fasta format?
#from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord



####################
# Set some globals #
####################

SAMTOOLS_DIR=""
SSAHA_DIR="/nfs/users/nfs_s/sh16/ssaha2_v2.5.1_x86_64/"
BWA_DIR=""
MY_SCRIPTS_DIR="/nfs/pathogen/sh16_scripts/"


#################################
# Simple Error Printing Funtion #
#################################

def DoError(ErrorString):
	sys.stderr.write("\n!!!Error: "+ErrorString+"!!!\n\n")
	sys.exit()


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	
	group = OptionGroup(parser, "IO Options")
	group.add_option("-b", "--bcf", action="store", dest="bcf", help="bcf/vcf file", default="", metavar="file")
	group.add_option("-v", "--vcf", action="store_true", dest="vcf", help="variation input file is in vcf format [default is bcf]", default=False)
	group.add_option("-o", "--output", action="store", dest="output", help="output file name (file type will be guessed from the name. e.g. plot.pdf will be saved in pdf forat) if no output name is given, the plot will be shown on the screen", default="", metavar="file name ")
#	group.add_option("-p", "--pseudosequence", action="store", dest="pseudosequence", help="output pseudosequence including indels [default is to write a psedoalignment to the reference and separate indel text file]", default=False)
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "Filtering Options")
	group.add_option("-d", "--depth", action="store", dest="depth", help="Minimum number of reads matching SNP [default= %default]", default=4, type="int", metavar="int")
	group.add_option("-D", "--stranddepth", action="store", dest="stranddepth", help="Minimum number of reads matching SNP per strand [default= %default]", default=2, type="int", metavar="int")
	group.add_option("-r", "--ratio", action="store", dest="ratio", help="Minimum ratio of first to second base call [default= %default]", default=0.8, type="float", metavar="float")
	#are the three above options still necessary???
	group.add_option("-q", "--QUAL", action="store", dest="QUAL", help="Minimum base quality [default= %default]", default=50.0, type="float", metavar="float")
	group.add_option("-m", "--MQUAL", action="store", dest="MQUAL", help="Minimum mapping quality [default= %default]", default=50.0, type="float", metavar="float")
	group.add_option("-a", "--AF1", action="store", dest="AF1", help="Minimum allele frequency (you would expect an AF of 1 for haploid SNPs). For non-SNP bases, the program will use 1- this number. [default= %default]", default=0.95, type="float")
	group.add_option("-c", "--CI95", action="store", dest="CI95", help="Maximum 95% confidence interval variation from AF. [default= %default]", default=0.0, type="float", metavar="float")
	group.add_option("-S", "--strandbias", action="store", dest="strand_bias", help="p-value cutoff for strand bias. [default= %default]", default=0.05, type="float", metavar="float")
	group.add_option("-Q", "--baseqbias", action="store", dest="baseq_bias", help="p-value cutoff for base quality bias. [default= %default]", default=0.05, type="float", metavar="float")
	group.add_option("-M", "--mappingbias", action="store", dest="mapping_bias", help="p-value cutoff for mapping bias. [default= %default]", default=0.05, type="float", metavar="float")
	group.add_option("-T", "--taildistancebias", action="store", dest="tail_bias", help="p-value cutoff for tail distance bias. [default= %default]", default=0.05, type="float", metavar="float")
	
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
		DoError("No input file specified")
	elif not os.path.isfile(options.bcf):
		DoError("Cannot find input file")

	
	
	if options.stranddepth<0:
		print "Minimum number of reads matching SNP on each strand must be >=0. Resetting to 0"
		options.stranddepth=0
	if options.depth<(options.stranddepth*2):
		print "Minimum number of reads matching SNP must be at least double that for each strand. Resetting to", options.stranddepth*2
		options.stranddepth=options.stranddepth*2
	if options.ratio<=0.5 or options.ratio>1:
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
	
	
	
	
	if options.vcf:
		try:
			bcffile=open(options.bcf, "rU")
		except StandardError:
			DoError("Cannot open vcf file")
	else:
		try:
			bcffile=os.popen("bcftools view "+options.bcf)
		except StandardError:
			DoError("Cannot open bcf file")
	
	atleastoneread=0
	mapped=0
	snps=0
	deletions=0
	deletion_lengths=0
	insertions=0
	insertion_lengths=0
	indels=[]
	heterocount=0
	basequalfail=0
	mapqualfail=0
	depthfail=0
	fdepthfail=0
	rdepthfail=0
	ratiofail=0
	fratiofail=0
	rratiofail=0
	AFfail=0
	strandbiasfail=0
	baseqbiasfail=0
	mappingbiasfail=0
	tailbiasfail=0
	
	
	remcount=0
	fremcount=0
	ratios=[]
	fratios=[]
#	fratios=[]
#	rratios=[]

	
	for line in bcffile:

		
	
		words=line.split()
			
		if words[0][0]=="#" or float(words[5])<options.QUAL:
			continue
		
		DP=[0,0,0,0]
		MQ=0
		info=words[7].split(";")
		for i in info:
			if len(i.split("="))==2 and i.split("=")[0]=="DP4":
				DP=map(float,i.split("=")[1].split(","))
			elif len(i.split("="))==2 and i.split("=")[0]=="MQ":
				MQ=float(i.split("=")[1])
				
		
		if MQ<options.MQUAL:
			continue

		#Calculate the ref/alt ratios
		
		DPratio=0.0
		
#		DPratiof=0.0
#		DPratior=0.0
		if (DP[0]+DP[1]+DP[2]+DP[3])>=options.depth and DP[0]+DP[2]>=options.stranddepth and DP[1]+DP[3]>=options.stranddepth:
			try: DPratio=(float(DP[0])+float(DP[1]))/(DP[0]+DP[1]+DP[2]+DP[3])
			except ZeroDivisionError:
				DPratio=1.0
			
			ratios.append(DPratio*100)
			if DPratio<options.ratio and DPratio>(1.0-options.ratio):
				remcount+=1
			
			try: DPfratio=(float(DP[0]))/(DP[0]+DP[2])
			except ZeroDivisionError:
				DPfratio=1.0
				
			try: DPrratio=(float(DP[1]))/(DP[1]+DP[3])
			except ZeroDivisionError:
				DPrratio=1.0
			if DPfratio<0.5:
				DPfcomp=DPfratio
			else:
				DPfcomp=1.0-DPfratio
			if DPrratio<0.5:
				DPrcomp=DPrratio
			else:
				DPrcomp=1.0-DPrratio
			
			if DPrcomp<DPfcomp:
				extremeratio=DPrratio
			else:
				extremeratio=DPfratio
			
			
			
			fratios.append(extremeratio*100)
			if extremeratio<options.ratio and extremeratio>(1.0-options.ratio):
				fremcount+=1
			
#		if DP[0]+DP[2]>=options.stranddepth:
#			try: DPratiof=(float(DP[0]))/(DP[0]+DP[2])
#			except ZeroDivisionError:
#				DPratiof=1.0
#			fratios.append(DPratiof*100)
#		if DP[1]+DP[3]>=options.stranddepth:
#			try: DPratior=(float(DP[1]))/(DP[1]+DP[3])
#			except ZeroDivisionError:
#				DPratior=1.0
#			rratios.append(DPratior*100)
		
		

		
	data=[]
	label=[]
	if len(ratios)>0:
		data.append(ratios)
		label.append("All")
#	if len(fratios)>0:
#		data.append(fratios)
#		label.append("Forward")
#	if len(rratios)>0:
#		data.append(rratios)
#		label.append("Reverse")
	
	if len(data)==0:
		DoError("No bases pass quality cutoffs")
	
#	pylab.Figure()
#
#	pylab.axvspan((1.0-options.ratio)*100, (options.ratio)*100, facecolor='r', alpha=0.2, hatch='/')
#	n, bins, patches = pylab.hist(data, bins=25, log=True, label=label) 
#	#pylab.legend(loc="best")
#	pylab.axis(ymin=0.1)
#	pylab.axvline(x=(options.ratio)*100, color="r", linestyle="--")
#	pylab.axvline(x=(1.0-options.ratio)*100, color="r", linestyle="--")
#	if remcount!=1:
#		pylab.figtext(0.15, 0.85, str(remcount)+" bases filtered due to heterozygosity")
#	else:
#		pylab.figtext(0.15, 0.85, str(remcount)+" base filtered due to heterozygosity")
#	pylab.title("Base heterozygosity plot of "+options.bcf) 
#	pylab.xlabel("%age of reads called as reference base") 
#	pylab.ylabel("Frequency")
#	if options.output=="":
#		pylab.show()
#	else:
#		pylab.savefig(options.output)
		
	data=[]	
	if len(fratios)>0:
		data.append(fratios)
		label.append("All")
		
	pylab.Figure()

	#pylab.axvspan((1.0-options.ratio)*100, (options.ratio)*100, facecolor='r', alpha=0.2, hatch='/')
	n, bins, patches = pylab.hist(data, bins=25, log=True, label=label, alpha=0.7) 
	#pylab.legend(loc="best")
	pylab.axis(ymin=0.1)
	#pylab.axvline(x=(options.ratio)*100, color="r", linestyle="--")
	#pylab.axvline(x=(1.0-options.ratio)*100, color="r", linestyle="--")
	if remcount!=1:
		pylab.figtext(0.15, 0.8, str(fremcount)+" bases filtered due to heterozygosity")
	else:
		pylab.figtext(0.15, 0.8, str(fremcount)+" base filtered due to heterozygosity")
	pylab.title("Base heterozygosity plot of "+options.bcf) 
	pylab.xlabel("%age of reads called as reference base") 
	pylab.ylabel("Frequency")
	if options.output=="":
		pylab.show()
	else:
		pylab.savefig(options.output)
		
		
		
		
		
