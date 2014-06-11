#!/usr/bin/env python

SMALT_DIR=""
SAMTOOLS_DIR=""
BCFTOOLS_DIR=""

##################
# Import modules #
##################

import os, sys, string
from random import randint, choice
from optparse import OptionParser
import pysam
from numpy import min, max, median, mean, std
from math import sqrt, pow

##############################################
## Function to reverse complement a sequence #
##############################################

def revcomp(sequence):
	rev=sequence[::-1]
	revcomp=''
	d={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'g':'c', 'c':'g', "n":"n", "N":"N"}
	for i in rev:
		if d.has_key(i):
			revcomp=revcomp+d[i]
		else:
			revcomp=revcomp+i
	
	return revcomp



##############################
# Get command line arguments #
##############################

def get_user_options():
	usage = "usage: %prog [options]"
	version="%prog 1.0. Written by Simon Harris, Wellcome Trust Sanger Institute, 2012"
	parser = OptionParser(usage=usage, version=version)
	
	#do not allow arguments to be interspersed. e.g. -a -b arg1 agr2. MUST be -a arg1 -b arg2.
	parser.disable_interspersed_args()
	

	parser.add_option("-s", "--serotypes", action="store", dest="serotypes", help="multifasta containing example serotypes", default="", metavar="FILE")
	parser.add_option("-o", "--output", action="store", dest="output", help="prefix for output files", default="", metavar="FILE")
	parser.add_option("-f", "--forward", action="store", dest="forward", help="forward fastq file (may be zipped, but must end .gz)", default="", metavar="FILE")
	parser.add_option("-r", "--reverse", action="store", dest="reverse", help="reverse fastq file (may be zipped, but must end .gz)", default="", metavar="FILE")
	parser.add_option("-m", "--map", action="store_false", dest="map", help="Do not map (bam file exists)", default=True)
	parser.add_option("-b", "--bam", action="store", dest="bam", help="Name of bam file", default="", metavar="FILE")
	parser.add_option("-c", "--cutoff", action="store", dest="cutoff", help="Cutoff for including secondary hits. Any hit within this %age of the best hit will be included. [default=%default]", default=10, type=float)
	

	return parser.parse_args()

(options, args)=get_user_options()

if not os.path.isfile(options.forward):
	print "Could not find forward fastq file"
	sys.exit()
if not os.path.isfile(options.reverse):
	print "Could not find reverse fastq file"
	sys.exit()
if options.bam!="" and not os.path.isfile(options.bam):
	print "Could not find bam file"
	sys.exit()
if options.cutoff>100 or options.cutoff<0:
	print "Cutoff (-c) must be a percentage between 0 and 100"
	sys.exit()
if options.output=="":
	print "No prefix for output files (-o) provided"
	sys.exit()
if not os.path.isfile(options.serotypes):
	print "Could not find serotypes file"
	sys.exit()
try:
	contigsfile=open(options.serotypes, "rU").read()
except StandardError:
	print "Could not open serotypes file"
	sys.exit()




chars = string.ascii_letters + string.digits
tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))

if options.map:

	#funzip=False
	#runzip=False
	#if options.forward.split(".")[-1]=="gz":
	#	print "Unzipping forward fastq file"
	#	os.system("zcat "+options.forward+" > "+tmpname+"_1.fastq")
	#	forward=tmpname+"_1.fastq"
	#	funzip=True
	#else:
	#	forward=options.forward
	#if options.reverse.split(".")[-1]=="gz":
	#	print "Unzipping reverse fastq file"
	#	os.system("zcat "+options.reverse+" > "+tmpname+"_2.fastq")
	#	reverse=tmpname+"_2.fastq"
	#	runzip=True
	#else:
	#	reverse=options.reverse
	
	forward=options.forward
	reverse=options.reverse
	
	if not os.path.isfile(options.serotypes+".fai"):
		os.system(SAMTOOLS_DIR+"samtools faidx "+options.serotypes)
	if not os.path.isfile(options.serotypes+".index.sma") or not os.path.isfile(options.serotypes+".index.smi"):
		os.system(SMALT_DIR+"smalt index -k 13 -s 1 "+options.serotypes+".index "+options.serotypes)
	 
	os.system(SMALT_DIR+"smalt map  -y 0.5 -d 0 -f samsoft -o "+tmpname+".sam "+options.serotypes+".index "+forward+" "+reverse)
	#os.system(SMALT_DIR+"smalt map  -y 0.5 -d -1 -f samsoft -o "+tmpname+".sam "+options.serotypes+".index "+forward+" "+reverse)
	os.system(SAMTOOLS_DIR+"samtools view -F 4 -b -S "+tmpname+".sam -t "+options.serotypes+".fai > "+tmpname+".1.bam")
	os.system("rm -f "+tmpname+".sam")
	os.system(SAMTOOLS_DIR+"samtools sort "+tmpname+".1.bam "+options.output+"_mapping")
	os.system("rm -f "+tmpname+".1.bam")
	os.system(SAMTOOLS_DIR+"samtools index "+options.output+"_mapping.bam")
	
	
	
	#if funzip:
	#	os.system("rm -f "+forward)
	#if runzip:
	#	os.system("rm -f "+reverse)
else:
	print "Skipping mapping"
	if options.bam!="" and os.path.isfile(options.bam):
		os.system("cp "+options.bam+" "+options.output+"_mapping.bam")
		os.system(SAMTOOLS_DIR+"samtools index "+options.output+"_mapping.bam")
	if not os.path.isfile(options.output+"_mapping.bam"):
		print "Could not find "+options.output+"_mapping.bam"
		sys.exit()
	

os.system("~sh16/scripts/reportlabtest.py -d area -f 3 -l 1 "+options.output+"_mapping.bam "+options.serotypes+" -o "+options.output+"_serotype_mapping.pdf")

#
#try:
#	genesfile=open(options.serotypes, "rU").read()
#except StandardError:
#	print "Could not open contigs file"
#	sys.exit()
#	
#genes={}
#geneorder=[]
#for line in genesfile.split(">")[1:]:
#	bits=line.split("\n")
#	geneorder.append(bits[0].split()[0])
#	genes[bits[0].split()[0]]=''.join(bits[1:]).upper()
#
try:
	samfile = pysam.Samfile( options.output+"_mapping.bam", "rb" )
except StandardError:
	print "Failed to read bam file"
	sys.exit()

samrefs=samfile.references
samlengths=samfile.lengths
samfile.close()
#genedepths={}
#geneerrors={}
#geneinsertions={}
#genedeletions={}
#
#for x, ref in enumerate(samrefs):
#	genedepths[ref]=[0]*samlengths[x]
#	geneerrors[ref]=[0]*samlengths[x]
#	geneinsertions[ref]=[0]*samlengths[x]
#	genedeletions[ref]=[0]*samlengths[x]
#count=0
#for read in samfile:
#
#	if not read.is_unmapped:# and not read.is_reverse:
#		start=read.pos
#		readpos=0
#		refpos=start
#		refname=samfile.getrname(read.rname)
#		readseq=read.seq.upper()
#		
#		for cig in read.cigar:
#			if cig[0]==0:
#				for x in range(0,cig[1]):
#					if readseq[readpos]!=genes[refname][refpos] and genes[refname][refpos] in ["A", "C", "G", "T"]:
##						print readseq[readpos], genes[refname][refpos]
#						geneerrors[refname][refpos]+=1
#					readpos+=1
#					
#					genedepths[refname][refpos]+=1
#					refpos+=1
#			elif cig[0]==1:
#				geneinsertions[refname][refpos]+=1
##				refstats[refname]["insertions"]+=1
#				readpos+=cig[1]
#			elif cig[0]==2:
#				genedeletions[refname][refpos]+=1
##				refstats[refname]["deletions"]+=1
#				refpos+=cig[1]
#			elif cig[0]==4:
#				readpos+=cig[1]
#			else:
#				print cig
#
#reforder=[]
#for ref in samrefs:
#	reforder.append([100*(float(len(genedepths[ref])-genedepths[ref].count(0))/len(genedepths[ref])),ref])
#reforder.sort()
#reforder.reverse()
#
#
#toprint=[]
#
#for x in reforder:
#	ref=x[1]
#	if x==reforder[0]:
#		toprint=[options.forward.split("/")[-1].replace(".gz","").replace("_1.fastq",""), ref, 100*(float(len(genedepths[ref])-genedepths[ref].count(0))/len(genedepths[ref])), min(genedepths[ref]), max(genedepths[ref]), mean(genedepths[ref]), std(genedepths[ref])/mean(genedepths[ref]), []]
#		refmatch=100*(float(len(genedepths[ref])-genedepths[ref].count(0))/len(genedepths[ref]))
#	elif x[0]>refmatch-options.cutoff:
#		toprint[-1].append(','.join(map(str,[ref, 100*(float(len(genedepths[ref])-genedepths[ref].count(0))/len(genedepths[ref])), mean(genedepths[ref])])))
#		
#	
#output=open(options.output+"_serotyping.txt", "w")
#print >> output, '\t'.join(["File", "Best Hit", "% Covered", "Min Coverage", "Max Coverage", "Mean Coverage", "Coverage CV", "Secondary Hits"])
#toprint[-1]='; '.join(toprint[-1])
#print >> output, '\t'.join(map(str,toprint))
#
#output.close()
#samfile.close()	

os.system("samtools mpileup -d 1000 -DSugBf  "+options.serotypes+" "+options.output+"_mapping.bam > "+tmpname+".mpileup")
os.system("bcftools view -bcg "+tmpname+".mpileup > "+options.output+".bcf")
os.system("rm -rf "+tmpname+".mpileup")


try:
	bcffile=os.popen(BCFTOOLS_DIR+"bcftools view "+options.output+".bcf")
except StandardError:
	print "Cannot open bcf file"
	sys.exit()

serotypeinfo={}
for x, ref in enumerate(samrefs):
	serotypeinfo[ref]={"length": float(samlengths[x]), "match": 0.0, "SNP": 0, "q_base_coverage": 0.0, "bases_covered": 0.0, "base_coverage": 0.0, "depth": 0.0, "heterozygous": 0.0, "insertion": 0.0, "deletion": 0.0}
	
for line in bcffile:
	line=line.strip()
	if len(line)>0 and line[0]=="#":
		continue
	words=line.split()
	serotype=words[0]
	position=words[1]
	refbase=words[3]
	altbase=words[4]
	quality=float(words[5])
	info=words[7].split(";")
	baseinfo={}
	for i in info:
		infotype=i.split("=")[0]
					
		if len(i.split("="))<2:
			if infotype=="INDEL":
				baseinfo[infotype]=True
		else:
			infodata=i.split("=")[1]
			try: baseinfo[infotype]=float(infodata)
			except StandardError:
				try: baseinfo[infotype]=map(float,infodata.split(","))
				except StandardError:
					baseinfo[infotype]=infodata
	if "INDEL" in baseinfo:
		continue
	if "AF1" in baseinfo and baseinfo["AF1"]!=0 and baseinfo["AF1"]!=1:
		serotypeinfo[serotype]["heterozygous"]+=1
	elif "DP4" in baseinfo and baseinfo["DP4"][0]+baseinfo["DP4"][1]<baseinfo["DP4"][2]+baseinfo["DP4"][3]:
		serotypeinfo[serotype]["SNP"]+=1
	elif "DP4" in baseinfo and baseinfo["DP4"][0]+baseinfo["DP4"][1]>baseinfo["DP4"][2]+baseinfo["DP4"][3]:
		serotypeinfo[serotype]["match"]+=1
	if "DP4" in baseinfo:
		serotypeinfo[serotype]["q_base_coverage"]+=(baseinfo["DP4"][0]+baseinfo["DP4"][1]+baseinfo["DP4"][2]+baseinfo["DP4"][3])
	
	if "DP" in baseinfo:
		serotypeinfo[serotype]["base_coverage"]+=baseinfo["DP"]
		if baseinfo["DP"]>0:
			serotypeinfo[serotype]["bases_covered"]+=1

serotypepercentid=[]
for serotype in serotypeinfo:
	serotypepercentid.append([serotypeinfo[serotype]["match"]/serotypeinfo[serotype]["length"], serotypeinfo[serotype]["base_coverage"], serotype])
	
serotypepercentid.sort()
serotypepercentid.reverse()
toprint=[]
for serotype in serotypepercentid:
	ref=serotype[2]
	if ref==serotypepercentid[0][2]:
		try:
			part4=serotypeinfo[ref]["base_coverage"]/serotypeinfo[ref]["bases_covered"]
		except ZeroDivisionError:
			part4=0
		toprint=[options.forward.split("/")[-1].replace(".gz","").replace("_1.fastq",""), ref, 100*serotype[0], 100*(serotypeinfo[ref]["bases_covered"]/serotypeinfo[ref]["length"]),  part4,  serotypeinfo[ref]["base_coverage"]/serotypeinfo[ref]["length"], serotypeinfo[ref]["SNP"], []]
		refmatch=100*serotype[0]
	elif 100*serotype[0]>refmatch-options.cutoff:
		
		try:
			part4=serotypeinfo[ref]["base_coverage"]/serotypeinfo[ref]["bases_covered"]
		except ZeroDivisionError:
			part4=0
		
		toprint[-1].append(','.join(map(str,[ref, 100*serotype[0], 100*(serotypeinfo[ref]["bases_covered"]/serotypeinfo[ref]["length"]),  part4,  serotypeinfo[ref]["base_coverage"]/serotypeinfo[ref]["length"], serotypeinfo[ref]["SNP"]])))
	
	
	

#
#for x in reforder:
#	ref=x[1]
#	if x==reforder[0]:
#		toprint=[options.forward.split("/")[-1].replace(".gz","").replace("_1.fastq",""), ref, 100*(float(len(genedepths[ref])-genedepths[ref].count(0))/len(genedepths[ref])), min(genedepths[ref]), max(genedepths[ref]), mean(genedepths[ref]), std(genedepths[ref])/mean(genedepths[ref]), []]
#		refmatch=100*(float(len(genedepths[ref])-genedepths[ref].count(0))/len(genedepths[ref]))
#	elif x[0]>refmatch-options.cutoff:
#		toprint[-1].append(','.join(map(str,[ref, 100*(float(len(genedepths[ref])-genedepths[ref].count(0))/len(genedepths[ref])), mean(genedepths[ref])])))

output=open(options.output+"_serotyping.txt", "w")
print >> output, '\t'.join(["File", "Best Hit", "%ID", "% Covered", "Mean Coverage of covered bases", "Mean Coverage of all bases", "SNPs", "Secondary Hits"])
toprint[-1]='; '.join(toprint[-1])
print >> output, '\t'.join(map(str,toprint))
output.close()
sys.exit()




