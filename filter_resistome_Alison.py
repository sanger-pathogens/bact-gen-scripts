#!/usr/bin/env python

SMALT_DIR=""
SAMTOOLS_DIR=""

##################
# Import modules #
##################

import os, sys, string
from random import randint, choice
from optparse import OptionParser
import pysam
from numpy import min, max, median, mean, std
from scipy.stats import mannwhitneyu, ttest_ind
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
	

	parser.add_option("-c", "--contigs", action="store", dest="contigs", help="multifasta containing contigs to search in", default="", metavar="FILE")
	parser.add_option("-g", "--genes", action="store", dest="genes", help="multifasta containing genes to search for", default="", metavar="FILE")
	parser.add_option("-b", "--bamfile", action="store", dest="bamfile", help="bamfile of mapped genes", default="", metavar="FILE")
	parser.add_option("-o", "--output", action="store", dest="output", help="prefix for output files", default="", metavar="FILE")
	parser.add_option("-f", "--forward", action="store", dest="forward", help="forward fastq file (may be zipped, but must end .gz)", default="", metavar="FILE")
	parser.add_option("-r", "--reverse", action="store", dest="reverse", help="reverse fastq file (may be zipped, but must end .gz)", default="", metavar="FILE")
	parser.add_option("-i", "--id", action="store", dest="id", help="minimum id to report match (excluding clipping due to contig breaks) [default = %default]", default=0.9, type="float", metavar="float")
	

	return parser.parse_args()

(options, args)=get_user_options()

if not os.path.isfile(options.contigs):
	print "Could not find contigs file"
	sys.exit()	
if not os.path.isfile(options.bamfile):
	print "Could not find bam file"
	sys.exit()	
if not os.path.isfile(options.genes):
	print "Could not find genes file"
	sys.exit()
try:
	contigsfile=open(options.contigs, "rU").read()
except StandardError:
	print "Could not open contigs file"
	sys.exit()

if options.id<0 or options.id>1:
	print "percent id (-i) must be between 0 and 1"
	sys.exit()

contigs={}
genes_present=[]
for line in contigsfile.split(">")[1:]:
	bits=line.split("\n")
	contigs[bits[0].split()[0]]=''.join(bits[1:])

filename=options.bamfile


try:
	genesfile=open(options.genes, "rU").read()
except StandardError:
	print "Could not open contigs file"
	sys.exit()
	
genes={}
geneorder=[]
for line in genesfile.split(">")[1:]:
	bits=line.split("\n")
	geneorder.append(bits[0].split()[0])
	genes[bits[0].split()[0]]=''.join(bits[1:]).upper()
refstarts={}
refends={}
count=0
for x, ref in enumerate(geneorder):
	refstarts[ref]=count
	count+=len(genes[ref])
	refends[ref]=count


if filename.split(".")[-1]=="bam":
	samfile = pysam.Samfile( filename, "rb" )
elif filename.split(".")[-1]=="sam":
	samfile = pysam.Samfile( filename, "r" )
else:
	print filename+" not a bam file"
	sys.exit() 

refs=samfile.references
lengths=samfile.lengths
lengthdict={}
for x in range(0, len(refs)):
	lengthdict[refs[x]]=lengths[x]

for read in samfile:
	if read.is_unmapped:
		continue
	if read.is_reverse:
		strand="-"
	else:
		strand="+"
	start=read.pos
	readpos=0
	refpos=start
	insertions=0
	inslength=0
	deletions=0
	dellength=0
	SNPs=0
	clipped=0
	cliplength=0
	lcliplen=0
	rcliplen=0
	for cignum, cig in enumerate(read.cigar):
		  
		if cig[0]==0:
			for x in range(0,cig[1]):
				if read.seq[readpos].upper()!=contigs[refs[read.tid]][refpos].upper():
					SNPs+=1
				readpos+=1
				refpos+=1
		elif cig[0]==1:
			insertions+=1
			inslength+=cig[1]
			readpos+=cig[1]
		elif cig[0]==2:
			deletions+=1
			dellength+=cig[1]
			refpos+=cig[1]
		elif cig[0]==4:
			clipped+=1
			cliplength+=cig[1]
			if cignum==0:
				lcliplen+=cig[1]
			elif cignum==(len(read.cigar)-1):
				rcliplen+=cig[1]
			else:
				print "Internal clipping?!"
				print read.cigar
			readpos+=cig[1]
		else:
			print cig
	end=refpos
	   
	   
	adjustedcliplen=cliplength
	at_contig_break=0
	   
	if lcliplen>start:
		adjustedcliplen=adjustedcliplen-(lcliplen-start)
		at_contig_break+=1
	elif (rcliplen+end)>lengths[read.tid]:
		adjustedcliplen=adjustedcliplen-((rcliplen+end)-lengths[read.tid])
		at_contig_break+=1
		   
	matchlength=len(read.seq)-(cliplength)
	matchpercent=(float(matchlength)/len(read.seq))*100
   
	percentid=((float(len(read.seq))-(SNPs+adjustedcliplen))/len(read.seq))
#	   if read.is_reverse:
#	   	readseq=revcomp(read.seq)
#	   else:
#	   	readseq=read.seq
	   
	if percentid>=options.id:
			genes_present.append([read.qname, samfile.getrname(read.rname), lengthdict[samfile.getrname(read.rname)], start, end, strand, len(read.seq), SNPs, insertions, deletions, clipped, inslength, dellength, cliplength, at_contig_break, 100*percentid, matchlength, matchpercent, genes[read.qname]])
	   
	   #print filename, read.qname, len(read.seq), SNPs, insertions, deletions, clipped, inslength, dellength, cliplength

#for gene in genes_present:
#	print gene
if len(genes_present)>0:

	bestgene=genes_present[0]
	secondary_genes=[]
	secondary_fragments=[]
	outputlines=[]
	bestgene_sequences={}
	
	for gene in genes_present[1:]:
	    
	    #print gene
	
		if gene[1]==bestgene[1]:
			if (bestgene[4]-gene[3])>(bestgene[6]*0.5) or (bestgene[4]-gene[3])>(gene[6]*0.5):
				#print "overlap"
				if float(gene[7]+gene[8]+gene[9])/gene[6] < float(bestgene[7]+bestgene[8]+bestgene[9])/bestgene[6]:
					secondary_genes.append(';'.join(map(str, bestgene[:-1])))
					bestgene=gene
				else:
					secondary_genes.append(';'.join(map(str, gene[:-1])))
			else:
				outputlines.append(bestgene[:-1]+[', '.join(secondary_genes)])
				bestgene_sequences[bestgene[0]]=bestgene[-1]
		    	#print outputlines
				bestgene=gene
				secondary_genes=[]
				secondary_fragments=[]
		else:
			outputlines.append(bestgene[:-1]+[', '.join(secondary_genes)])
			bestgene_sequences[bestgene[0]]=bestgene[-1]
			bestgene=gene
			secondary_genes=[]
	
	outputlines.append(bestgene[:-1]+[', '.join(secondary_genes)])
	bestgene_sequences[bestgene[0]]=bestgene[-1]
	
	
	if options.output!="":
		output=open(options.output+"_hits.txt","w")
		plotlines=[]
		print >> output, "\t".join(["gene", "contig", "contig length", "start", "end", "strand", "length", "SNPs", "No. insertions", "No. deletions", "No. clipped regions", "total insertion length", "total deletion length", "clipped length", "No. contig breaks", "Percent id", "Match length", "Match length percent", "Overlapping secondary gene hits"])
		for line in outputlines:
			print >> output, '\t'.join(map(str,line))
	
			plotlines.append([refstarts[line[0]],refends[line[0]],((float(line[6])-float(line[7]+line[8]+line[9]))/line[6])*100])
	
		output.close()
#		sys.exit()
		plotout=open(options.output+"_hits.plot","w")
		print >> plotout, "#BASE MATCH"
		plotlines.sort()
		for line in plotlines:
			for x in xrange(line[0],line[1]):
				print >> plotout, x, line[2]
		plotout.close()
		
		if options.forward != "" and options.reverse!="":
			
			bestgeneorder=[]
			for gene in geneorder:
				if gene in bestgene_sequences:
					bestgeneorder.append(gene)
			
			accout=open(options.output+"_hits.mfa", "w")
			for gene in bestgeneorder:
				print >> accout, ">"+gene
				print >> accout, bestgene_sequences[gene]
			accout.close()
			
			#make random name for files
			chars = string.ascii_letters + string.digits
			tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))
			funzip=False
			runzip=False
			if options.forward.split(".")[-1]=="gz":
				print "Unzipping forward fastq file"
				os.system("zcat "+options.forward+" > "+tmpname+"_1.fastq")
				options.forward=tmpname+"_1.fastq"
				funzip=True
			if options.reverse.split(".")[-1]=="gz":
				print "Unzipping reverse fastq file"
				os.system("zcat "+options.reverse+" > "+tmpname+"_2.fastq")
				options.reverse=tmpname+"_2.fastq"
				runzip=True
			
			os.system(SAMTOOLS_DIR+"samtools faidx "+options.output+"_hits.mfa")
			os.system(SMALT_DIR+"smalt index -k 13 -s 1 "+options.output+"_hits.mfa.index "+options.output+"_hits.mfa")
			#os.system(SMALT_DIR+"smalt map -y "+str(options.id)+" -r 12345 -f samsoft -o "+options.output+".sam "+options.output+"_hits.mfa.index "+options.forward+" "+options.reverse)
			os.system(SMALT_DIR+"smalt map -y 0.5 -r 12345 -f samsoft -o "+options.output+".sam "+options.output+"_hits.mfa.index "+options.forward+" "+options.reverse)
			os.system(SAMTOOLS_DIR+"samtools view -F 4 -b -S "+options.output+".sam -t "+options.output+"_hits.mfa.fai > "+options.output+".1.bam")
			os.system(SAMTOOLS_DIR+"samtools sort "+options.output+".1.bam "+options.output+"_mapping")
			os.system(SAMTOOLS_DIR+"samtools index "+options.output+"_mapping.bam")
			
			os.system("rm -f "+options.output+".sam "+options.output+".1.bam "+options.output+"_hits.mfa.index.smi "+options.output+"_hits.mfa.index.sma "+options.output+"_hits.mfa.fai")
			
			if funzip:
				os.system("rm -f "+options.forward)
			if runzip:
				os.system("rm -f "+options.reverse)
	
			
			
			try:
				samfile = pysam.Samfile( options.output+"_mapping.bam", "rb" )
			except StandardError:
				print "Failed to read bam file"
				sys.exit()
				
			samrefs=samfile.references
			samlengths=samfile.lengths
			genedepths={}
			geneerrors={}
			
			for x, ref in enumerate(samrefs):
				genedepths[ref]=[0]*samlengths[x]
				geneerrors[ref]=[0]*samlengths[x]
			count=0
			for read in samfile:
			
		
				if not read.is_unmapped:# and not read.is_reverse:
					start=read.pos
					readpos=0
					refpos=start
					refname=samfile.getrname(read.rname)
					readseq=read.seq.upper()
	#				
	#				if refname=="fusC":
	#					print read.is_reverse
	#					print " "*start+readseq
	#					print bestgene_sequences[refname]
	#					print genes[refname]
	#					count+=1
	#				if count==10:	
	#					sys.exit()
					
					for cig in read.cigar:
						if cig[0]==0:
							for x in range(0,cig[1]):
								if readseq[readpos]!=genes[refname][refpos] and genes[refname][refpos] in ["A", "C", "G", "T"]:
	#								print readseq[readpos], genes[refname][refpos]
									geneerrors[refname][refpos]+=1
								readpos+=1
								
								genedepths[refname][refpos]+=1
								refpos+=1
						elif cig[0]==1:
	#						insertions+=1
	#						refstats[refname]["insertions"]+=1
							readpos+=cig[1]
						elif cig[0]==2:
	#						deletions+=1
	#						refstats[refname]["deletions"]+=1
							refpos+=cig[1]
						elif cig[0]==4:
							readpos+=cig[1]
						else:
							print cig
				
		
			samfile.close()	
			count=0
			filecount=0
			maxseqlen=0
			tmpfilelist=[]
			output=open(options.output+"_coverage.plot", "w")
			print >> output, "#BASE COVERAGE"
			for gene in bestgeneorder:
				filecount+=1
				
				tmpseqout=open(tmpname+".fasta", "w")
	#			tmpfilelist.append(tmpname+".fasta")
				print >> tmpseqout, ">"+str(gene)
				print >> tmpseqout, genes[gene]
				tmpseqout.close()
				tmpout=open(tmpname+".plot", "w")
				print >> tmpout, "#BASE COVERAGE SNPs"
	#			tmpfilelist.append(tmpname+str(filecount)+".plot")
				if len(genedepths[gene])>maxseqlen:
					maxseqlen=len(genedepths[gene])
				for x in xrange(len(genedepths[gene])):
					count+=1
					if genedepths[gene][x]!=0 or geneerrors[gene][x]!=0 or x==len(genedepths[gene])-1:
						print >> output, count, genedepths[gene][x], geneerrors[gene][x]
						print >> tmpout, x+1, genedepths[gene][x], geneerrors[gene][x]
				tmpout.close()
				
				os.system(SAMTOOLS_DIR+"~sh16/scripts/reportlabtest.py -H 6 -w -d area -4 "+str(len(genedepths[gene]))+"  -Y 0 -l 1 -o "+tmpname+str(filecount)+".pdf "+tmpname+".plot "+tmpname+".fasta")
				tmpfilelist.append(tmpname+str(filecount)+".pdf")
				
				os.system("rm -f "+tmpname+".plot "+tmpname+".fasta")
				
			output.close()
			
			
			os.system("gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile="+options.output+"_per_gene_SNP_plot.pdf "+' '.join(tmpfilelist))
			#os.system(SAMTOOLS_DIR+"~sh16/scripts/reportlabtest.py -H 4 -d area -4 "+str(maxseqlen)+" -O portrait -Y 0 -l 4 -o "+options.output+"_accessory_mapping_test.pdf "+' '.join(tmpfilelist[::-1]))
			
			os.system("rm -f "+' '.join(tmpfilelist))
			
			os.system(SAMTOOLS_DIR+"~sh16/scripts/reportlabtest.py -H 8 -w -d area -Y 0 -l 1 -o "+options.output+"_coverage_plot.pdf "+options.output+"_coverage.plot "+options.output+"_hits.mfa")
			
			count=0
			output=open(options.output+"_presence.plot", "w")
			print >> output, "#BASE PRESENT"
			for gene in geneorder:
				if gene in genedepths:
					for x in xrange(len(genedepths[gene])):
						count+=1
						if genedepths[gene][x]!=0 or geneerrors[gene][x]!=0:
							print >> output, count, 1
						elif x==len(genedepths[gene])-1:
							print >> output, count, 0
				else:
					count+=len(genes[gene])
					print >> output, count, 0
			
			output.close()
			
			output=open(options.output+"_coverage.txt", "w")
			print >> output, '\t'.join(['gene', 'min', 'max', 'median', 'mean', 'std'])
			print '\t'.join(['gene', 'min', 'max', 'median', 'mean', 'std'])
			for gene in bestgeneorder:
				print >> output, '\t'.join(map(str,[gene, min(genedepths[gene]), max(genedepths[gene]), median(genedepths[gene]), mean(genedepths[gene]), std(genedepths[gene])]))
				
	#			MLSTdata=genedepths[bestgeneorder[-7]]+genedepths[bestgeneorder[-6]]+genedepths[bestgeneorder[-5]]+genedepths[bestgeneorder[-4]]+genedepths[bestgeneorder[-3]]+genedepths[bestgeneorder[-2]]+genedepths[bestgeneorder[-1]]
	#			
	#			pooledsd=sqrt((((len(genedepths[gene])-1)*pow(std(genedepths[gene]),2))+((len(MLSTdata)-1)*pow(std(MLSTdata),2)))/(len(genedepths[gene])+len(MLSTdata)))
	#			
	#			print pooledsd
	#			print (mean(genedepths[gene])-mean(genedepths[bestgeneorder[-1]]))/pooledsd
	#			
				print '\t'.join(map(str,[gene, min(genedepths[gene]), max(genedepths[gene]), median(genedepths[gene]), mean(genedepths[gene]), std(genedepths[gene])]))
	#			print mean(genedepths[gene]), std(genedepths[gene]), mean(genedepths[bestgeneorder[-4]]), std(genedepths[bestgeneorder[-4]]), mannwhitneyu(genedepths[gene], genedepths[bestgeneorder[-7]]+genedepths[bestgeneorder[-6]]+genedepths[bestgeneorder[-5]]+genedepths[bestgeneorder[-4]]+genedepths[bestgeneorder[-3]]+genedepths[bestgeneorder[-2]]+genedepths[bestgeneorder[-1]])
	#			print mean(genedepths[gene]), std(genedepths[gene]), mean(genedepths[bestgeneorder[-4]]), std(genedepths[bestgeneorder[-4]]), mannwhitneyu(numpy.array(genedepths[gene]), numpy.array(genedepths[bestgeneorder[-7]]+genedepths[bestgeneorder[-6]]+genedepths[bestgeneorder[-5]]+genedepths[bestgeneorder[-4]]+genedepths[bestgeneorder[-3]]+genedepths[bestgeneorder[-2]]+genedepths[bestgeneorder[-1]]))
	#			print mean(genedepths[gene]), std(genedepths[gene]), mean(genedepths[bestgeneorder[-4]]), std(genedepths[bestgeneorder[-4]]), ttest_ind(genedepths[gene], genedepths[bestgeneorder[-7]]+genedepths[bestgeneorder[-6]]+genedepths[bestgeneorder[-5]]+genedepths[bestgeneorder[-4]]+genedepths[bestgeneorder[-3]]+genedepths[bestgeneorder[-2]]+genedepths[bestgeneorder[-1]])
				
	    	output.close()
	else:
		print "\t".join(["gene", "contig", "contig length", "start", "end", "strand", "length", "SNPs", "No. insertions", "No. deletions", "No. clipped regions", "total insertion length", "total deletion length", "clipped length", "No. contig breaks", "Percent id", "Match length", "Match length percent", "Overlapping secondary gene hits"])
		for line in outputlines:
			print '\t'.join(map(str,line))

else:
	print "Not hits found"
	count=0
	output=open(options.output+"_presence.plot", "w")
	print >> output, "#BASE PRESENT"
	for gene in geneorder:
		count+=len(genes[gene])
		print >> output, count, 0
	
	output.close()
