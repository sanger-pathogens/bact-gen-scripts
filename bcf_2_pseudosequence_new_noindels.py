#!/usr/bin/env python
import string, re
import os, sys
from optparse import OptionParser, OptionGroup
#probably should add this to make sequences in real fasta format?
#from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord



####################
# Set some globals #
####################

SAMTOOLS_DIR=""
BCFTOOLS_DIR=""
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
	group.add_option("-b", "--bcf", action="store", dest="bcf", help="bcf/vcf file", default="", metavar="file")
	group.add_option("-v", "--vcf", action="store_true", dest="vcf", help="variation input file is in vcf format [default is bcf]", default=False)
	group.add_option("-B", "--bam", action="store", dest="bam", help="bam/sam file", default="", metavar="file")
	group.add_option("-s", "--sam", action="store_true", dest="sam", help="bam/sam input file is in sam format [default is bam]", default=False)
	group.add_option("-o", "--output", action="store", dest="output", help="prefix for output files", default="", metavar="prefix")
#	group.add_option("-p", "--pseudosequence", action="store", dest="pseudosequence", help="output pseudosequence including indels [default is to write a psedoalignment to the reference and separate indel text file]", default=False)
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "Filtering Options")
	group.add_option("-d", "--depth", action="store", dest="depth", help="Minimum number of reads matching SNP [default= %default]", default=4, type="int", metavar="int")
	group.add_option("-D", "--stranddepth", action="store", dest="stranddepth", help="Minimum number of reads matching SNP per strand [default= %default]", default=2, type="int", metavar="int")
	group.add_option("-r", "--ratio", action="store", dest="ratio", help="Minimum ratio of first to second base call [default= %default]", default=0.8, type="float", metavar="float")
	group.add_option("-R", "--qratio", action="store", dest="qratio", help="Minimum ratio of total depth that is good quality [default= %default]", default=0.5, type="float", metavar="float")
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
		DoError("No input file specified")
	elif not os.path.isfile(options.bcf):
		DoError("Cannot find input file")
	
	if options.output=="":
		DoError("No output prefix specified")
		
	if options.bam=="":
		DoError("sam or bam file from which bcf was made must be specified")
	elif options.sam:
		header=os.popen(SAMTOOLS_DIR+"samtools view -S -H "+options.bam).readlines()
	else:
		header=os.popen(SAMTOOLS_DIR+"samtools view -H "+options.bam).readlines()
	
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
	
	
	
	contigsizes={}
	contigstart={}
	contigorder=[]
	totallength=0
	
	for line in header:
		words=line.split()
		name=""
		length=0
		for word in words:
			if word.split(":")[0]=="SN":
				name=word.split(":")[1]
			elif word.split(":")[0]=="LN":
				length=int(word.split(":")[1])
		if name!="" and length!=0:
			contigsizes[name]=length
			contigstart[name]=totallength
			contigorder.append(name)
			totallength+=length
				
		
	if len(contigsizes)==0:
		DoError("No contigs found. Perhaps your sam/bam has no header?")
		
	
	contigs={}
	
	for contig in contigorder:
		contigs[contig]=["N"]*contigsizes[contig]
	
	if options.vcf:
		try:
			bcffile=open(options.bcf, "rU")
		except StandardError:
			DoError("Cannot open vcf file")
	else:
		try:
			bcffile=os.popen(BCFTOOLS_DIR+"bcftools view "+options.bcf)
		except StandardError:
			DoError("Cannot open bcf file")
	
	atleastoneread=0
	mapped=0
	snps=0
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
	qratiofail=0
	
	count=0
	total=0.0
	hundredth=float(totallength)/100
	
	skip=0
	
	
	for line in bcffile:
		
		
	
		words=line.split()
			
		if words[0][0]=="#":
			if words[0][1]!="#":
			
				headings=words
				headings[0]=headings[0][1:]
			continue

		
		if len(words)!=len(headings):
			print "words not equal to headings"
			print headings
			print words
			sys.exit()
		
		BASEINFO={}
		
		for x, heading in enumerate(headings):
		
			if heading=="INFO":
				
				BASEINFO[heading]={}
				
				try: info=words[x].split(";")
				except StandardError:
					print "Cannot split info string", words[x]
					sys.exit()
				for i in info:
					
					infotype=i.split("=")[0]
					
					if len(i.split("="))<2:
						if infotype=="INDEL":
							BASEINFO[heading][infotype]=True
					else:
						infodata=i.split("=")[1]
						try: BASEINFO[heading][infotype]=float(infodata)
						except StandardError:
							try: BASEINFO[heading][infotype]=map(float,infodata.split(","))
							except StandardError:
								BASEINFO[heading][infotype]=infodata
				
				
					
			else:
				try: BASEINFO[heading]=float(words[x])
				except StandardError:
					BASEINFO[heading]=words[x]
		
		if "INDEL" in BASEINFO["INFO"]:
			continue
		

		count+=1
		if count>=hundredth:
			total=float(BASEINFO["POS"])
			count=0
			print "%.0f%% complete\r" % (100*(total/totallength)),
			sys.stdout.flush()

		
		
		
		
		
		#Calculate the ref/alt ratios
		
		BASEINFO["INFO"]["DP4ratios"]={}
		if not "DP4" in BASEINFO["INFO"]:
			BASEINFO["INFO"]["DP4"]=[0,0,0,0]
			BASEINFO["INFO"]["DP4ratios"]["fref"]=0.0
			BASEINFO["INFO"]["DP4ratios"]["rref"]=0.0
			BASEINFO["INFO"]["DP4ratios"]["falt"]=0.0
			BASEINFO["INFO"]["DP4ratios"]["ralt"]=0.0
			BASEINFO["INFO"]["DP4rratio"]=1.0
			BASEINFO["INFO"]["AF1"]=0
			BASEINFO["INFO"]["MQ"]=0
		elif "DP4" in BASEINFO["INFO"]:
			try: BASEINFO["INFO"]["DP4ratios"]["fref"]=float(BASEINFO["INFO"]["DP4"][0])/(BASEINFO["INFO"]["DP4"][0]+BASEINFO["INFO"]["DP4"][2])
			except ZeroDivisionError:
				BASEINFO["INFO"]["DP4ratios"]["fref"]=0.0
			try: BASEINFO["INFO"]["DP4ratios"]["rref"]=float(BASEINFO["INFO"]["DP4"][1])/(BASEINFO["INFO"]["DP4"][1]+BASEINFO["INFO"]["DP4"][3])
			except ZeroDivisionError:
				BASEINFO["INFO"]["DP4ratios"]["rref"]=0.0
			try: BASEINFO["INFO"]["DP4ratios"]["falt"]=float(BASEINFO["INFO"]["DP4"][2])/(BASEINFO["INFO"]["DP4"][0]+BASEINFO["INFO"]["DP4"][2])
			except ZeroDivisionError:
				BASEINFO["INFO"]["DP4ratios"]["falt"]=0.0
			try: BASEINFO["INFO"]["DP4ratios"]["ralt"]=float(BASEINFO["INFO"]["DP4"][3])/(BASEINFO["INFO"]["DP4"][1]+BASEINFO["INFO"]["DP4"][3])
			except ZeroDivisionError:
				BASEINFO["INFO"]["DP4ratios"]["ralt"]=0.0
			try: BASEINFO["INFO"]["DP4rratio"]=(float(BASEINFO["INFO"]["DP4"][0])+float(BASEINFO["INFO"]["DP4"][1]))/(BASEINFO["INFO"]["DP4"][0]+BASEINFO["INFO"]["DP4"][1]+BASEINFO["INFO"]["DP4"][2]+BASEINFO["INFO"]["DP4"][3])
			except ZeroDivisionError:
				BASEINFO["INFO"]["DP4rratio"]=1.0
		
		
		atleastoneread+=1
		
		#filter the call
		
		keep=True
		SNP=True
		INDEL=False
		failedfilters=[]
		
		if BASEINFO["ALT"]=="." or BASEINFO["INFO"]["DP4rratio"]>=0.5:
			SNP=False
		
		if BASEINFO["QUAL"]<options.QUAL:
			keep=False
			basequalfail+=1
			failedfilters.append("Q"+str(options.QUAL))
		if  BASEINFO["INFO"]["MQ"]<options.MQUAL:
			keep=False
			mapqualfail+=1
			failedfilters.append("MQ"+str(options.MQUAL))
		if not SNP and BASEINFO["INFO"]["DP4"][0]+BASEINFO["INFO"]["DP4"][1]<options.depth:
			keep=False
			depthfail+=1
			failedfilters.append("D"+str(options.depth))
		if not SNP and BASEINFO["INFO"]["DP4"][0]+BASEINFO["INFO"]["DP4"][1]<(options.qratio*BASEINFO["INFO"]["DP"]):
			keep=False
			qratiofail+=1
			failedfilters.append("QRR"+str(options.qratio))
		if not SNP and BASEINFO["INFO"]["DP4"][0]<options.stranddepth:
			keep=False
			fdepthfail+=1
			failedfilters.append("FD"+str(options.stranddepth))
		if not SNP and BASEINFO["INFO"]["DP4"][1]<options.stranddepth:
			keep=False
			rdepthfail+=1
			failedfilters.append("RD"+str(options.stranddepth))
		if not SNP and BASEINFO["INFO"]["DP4ratios"]["fref"]<options.ratio:
			keep=False
			fratiofail+=1
			failedfilters.append("FRR"+str(options.ratio))
		if not SNP and BASEINFO["INFO"]["DP4ratios"]["rref"]<options.ratio:
			keep=False
			rratiofail+=1
			failedfilters.append("RRR"+str(options.ratio))
		if SNP and BASEINFO["INFO"]["DP4"][2]+BASEINFO["INFO"]["DP4"][3]<options.depth:
			#print BASEINFO["POS"], "6"
			keep=False
			depthfail+=1
			failedfilters.append("D"+str(options.depth))
		if SNP and BASEINFO["INFO"]["DP4"][2]+BASEINFO["INFO"]["DP4"][3]<(options.qratio*BASEINFO["INFO"]["DP"]):
			#print BASEINFO["POS"], "6"
			keep=False
			qratiofail+=1
			failedfilters.append("QRA"+str(options.qratio))
		if SNP and BASEINFO["INFO"]["DP4"][2]<options.stranddepth:
			#print BASEINFO["POS"], "7"
			keep=False
			fdepthfail+=1
			failedfilters.append("FD"+str(options.stranddepth))
		if SNP and BASEINFO["INFO"]["DP4"][3]<options.stranddepth:
			#print BASEINFO["POS"], "8"
			keep=False
			rdepthfail+=1
			failedfilters.append("RD"+str(options.stranddepth))
		if SNP and BASEINFO["INFO"]["DP4ratios"]["falt"]<options.ratio:
			#print BASEINFO["POS"], "9"
			keep=False
			fratiofail+=1
			failedfilters.append("FRA"+str(options.ratio))
		if SNP and BASEINFO["INFO"]["DP4ratios"]["ralt"]<options.ratio:
			#print BASEINFO["POS"], "10"
			keep=False
			rratiofail+=1
			failedfilters.append("RRA"+str(options.ratio))
		if not SNP and "AF1" in BASEINFO["INFO"] and BASEINFO["INFO"]["AF1"]>(1-options.AF1):
			keep=False
			AFfail+=1
			failedfilters.append("AF"+str(1-options.AF1))
		if SNP and BASEINFO["INFO"]["AF1"]<options.AF1:
			#print BASEINFO["POS"], "13"
			keep=False
			AFfail+=1
			failedfilters.append("AF"+str(options.AF1))
		if SNP and "PV4" in BASEINFO["INFO"]:
			if BASEINFO["INFO"]["PV4"][0]<=options.strand_bias:
				#print BASEINFO["POS"], "14"
				keep=False
				strandbiasfail+=1
				failedfilters.append("SB"+str(options.strand_bias))
			if BASEINFO["INFO"]["PV4"][1]<=options.baseq_bias:
				#print BASEINFO["POS"], "15"
				keep=False
				baseqbiasfail+=1
				failedfilters.append("QB"+str(options.baseq_bias))
			if BASEINFO["INFO"]["PV4"][2]<=options.mapping_bias:
				#print BASEINFO["POS"], "16"
				keep=False
				mappingbiasfail+=1
				failedfilters.append("MQB"+str(options.mapping_bias))
			if BASEINFO["INFO"]["PV4"][3]<=options.tail_bias:
				#print BASEINFO["POS"], "17"
				keep=False
				tailbiasfail+=1
				failedfilters.append("TB"+str(options.tail_bias))
			
		
		
		#find hetrozygous SNP calls
		if len(BASEINFO["ALT"].split(","))>1:
			HETERO=True
			heterocount+=1
			BASEINFO["ALT"]=BASEINFO["ALT"].split(",")[0]
			#keep=False
		
		
		#make the pseudosequence
		if keep:
			
			if SNP:
				snpline=int(BASEINFO["INFO"]["DP"])
				snps+=1
				contigs[BASEINFO["CHROM"]][int(BASEINFO["POS"])-1]=BASEINFO["ALT"][0]
				mapped+=1
				
			else:
				contigs[BASEINFO["CHROM"]][int(BASEINFO["POS"])-1]=BASEINFO["REF"][0]
				mapped+=1
				snpline=0
			

	os.system("gzip -f "+options.output+"_filtered_mapping.plot")
	
	out=open(options.output+".mfa","w")
	
	name=options.output.split("/")[-1].split('.')[0]
	
	for x, contig in enumerate(contigorder):
		print >> out, ">"+name+"_"+''.join(contig)
		print >> out, ''.join(contigs[contig])
	out.close()
	
	out=open(options.output+"_SNP_filter.log","w")
	print >> out,'Total bases in reference:', totallength
	print >> out,'Mapped:%d (%.2f%%)' % (mapped, (float(mapped)/totallength)*100)
	print >> out,'SNPs:', snps
	print >> out,'Reasons for base call rejections:'
	print >> out,'  Base quality:', basequalfail
	print >> out,'  Mapping quality:', mapqualfail
	print >> out,'  Depth:', depthfail
	print >> out,'  Ratio of high quality reads:', qratiofail
	print >> out,'  Forward depth:', fdepthfail
	print >> out,'  Reverse depth:', rdepthfail
	print >> out,'  Forward SNP ratio:', fratiofail
	print >> out,'  Reverse SNP ratio:', rratiofail
	print >> out,'  Allele frequency:', AFfail
	print >> out,'Reasons for SNP specific base call rejections:'
	print >> out,'  Heterozygous SNP calls:', heterocount
	print >> out,'  Strand bias:', strandbiasfail
	print >> out,'  Base quality bias:', baseqbiasfail
	print >> out,'  Mapping bias:', mappingbiasfail
	print >> out,'  Tail distance bias:', tailbiasfail
	out.close()

