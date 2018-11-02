#!/usr/bin/env python
import string, re
import os, sys
from optparse import OptionParser
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
	print "!!!Error:", ErrorString,"!!!"
	sys.exit()


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-p", "--pileup", action="store", dest="pileup", help="pileup file", default="")
	parser.add_option("-b", "--bam", action="store", dest="bam", help="bam file", default="")
	parser.add_option("-s", "--sam", action="store", dest="sam", help="sam file", default="")
	parser.add_option("-d", "--depth", action="store", dest="depth", help="Minimum depth [default= %default]", default=1, type="int")
	parser.add_option("-q", "--quality", action="store", dest="quality", help="Minimum mapping Quality [default= %default]", default=120, type="int")
	parser.add_option("-Q", "--strandquality", action="store", dest="strandquality", help="Minimum per strand mapping Quality [default= %default]", default=60, type="int")
	parser.add_option("-r", "--ratio", action="store", dest="ratio", help="SNP/Mapping quality ratio cutoff [default= %default]", default=0.75, type="float")
	
	parser.add_option("-o", "--output", action="store", dest="output", help="output file prefix", default="")
	
	return parser.parse_args()


################
# Main program #
################		

if __name__ == "__main__":


	#Get command line arguments

	(options, args) = main()


	#Do some checking of the input files
	
	if options.pileup=="":
		DoError("No pileup file specified")
	
	if options.output=="":
		DoError("No output prefix specified")
		
	if options.sam=="" and options.bam=="":
		DoError("sam or bam file from which pileup was made must be specified")
	elif options.sam!="":
		header=os.popen(SAMTOOLS_DIR+"samtools view -S -H "+options.sam).readlines()
	else:
		header=os.popen(SAMTOOLS_DIR+"samtools view -H "+options.bam).readlines()
	
	contigsizes={}
	contigorder=[]
	
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
			contigorder.append(name)
				
		
	if len(contigsizes)==0:
		DoError("No contigs found. Perhaps your sam/bam has no header?")
		
	
	contigs={}
	
	for contig in contigorder:
		contigs[contig]=["N"]*contigsizes[contig]
	

	pileupfile=open(options.pileup, "rU")
#	printnext=False
	for line in pileupfile:
#		words=line.split()
#		if words[3].upper()in ["A", "C", "G", "T"] and int(words[4])>options.quality and (words[2].upper()==words[3].upper() or (float(words[5])/float(words[4]))>=options.ratio) and words[2]!="*":
#			contigs[words[0]][int(words[1])-1]=words[3]
	
	
		words=line.split()
		
#		if len(words)<10:
#			print "ERROR, a line in your pileup file is corrupted"
#			if words[2]==words[3]:
#				contigs[contig][snplocation-1]=refbase
#			else:
#				contigs[contig][snplocation-1]="N"
#			continue
			
		
		
		contig=words[0]
		snplocation=int(words[1])
		refbase=words[2]
		snpbase=words[3]
		readdepth=int(words[7])
		reads=words[8]
		qualities=words[9]
		
		if readdepth==0:
			contigs[contig][snplocation-1]='N'
			continue
		
		elif readdepth<options.depth:
			contigs[contig][snplocation-1]='N'
			continue
		
		if refbase=="*":#need to do more with these!!
#			if printnext==True:
			#print line
			continue
#		elif refbase not in ["A","C","G","T","a","c","t","g"]:
#			contigs[contig][snplocation-1]=='N'
#			continue
			
		
		
#		baselocs={"A":5, "C":6, "G":7, "T":8}
		
		#remove from here down to next commented line (which needs to be added back in) to get old version
#		locqual=0.0
#		if readdepth>=mindepth:# and tempseq[snplocation-1]!='?':
#			locqual=((int(words[baselocs[refbase]])-int(words[baselocs[refbase]+6]))*5)+((float(words[baselocs[refbase]+6]))*2.5)
#	
#
#		if locqual>=quality:

		qualx=0
		snpbasecount=[{"A":0, "G":0, "C":0, "T":0},{"A":0, "G":0, "C":0, "T":0},{"A":0, "G":0, "C":0, "T":0}]
		snpqualcount=[{"A":0, "G":0, "C":0, "T":0},{"A":0, "G":0, "C":0, "T":0},{"A":0, "G":0, "C":0, "T":0}]
		insertcount={}
		deletioncount={}
		refbasecount=[0,0,0]
		refqualcount=[0,0,0]
		totalqual=[0,0,0]
		
		#print reads, qualities
		
		x=0
		while x<len(reads):
			#print reads[x], x, qualx
			read=reads[x].upper()
			
			if reads[x] in ['a','c','g','t',',']:
				strand=2
			else:
				strand=1
			
			if read in [",","."]:
				refbasecount[strand]+=1
				refqualcount[strand]+=ord(qualities[qualx])-33
				totalqual[strand]+=ord(qualities[qualx])-33
				
				refbasecount[0]+=1
				refqualcount[0]+=ord(qualities[qualx])-33
				totalqual[0]+=ord(qualities[qualx])-33
				qualx+=1
			elif snpbasecount[0].has_key(read):
				snpbasecount[strand][read]+=1
				snpqualcount[strand][read]+=ord(qualities[qualx])-33
				totalqual[strand]+=ord(qualities[qualx])-33
				
				snpbasecount[0][read]+=1
				snpqualcount[0][read]+=ord(qualities[qualx])-33
				totalqual[0]+=ord(qualities[qualx])-33
				qualx+=1
			elif read=="^":
				x+=1
			elif read in ["-","+"]:
				x+=1
				try:
					int(reads[x+1])
				except:
					skip=int(reads[x])
				else:
					#print reads[x:x+2]
					skip=int(reads[x:x+2])+1
				
				x+=skip
#				if read=="-":
#					if not deletioncount.has_key(skip):
#						deletioncount[skip]=0
#					deletioncount[skip]+=1
#				elif read=="+":
#					if not insertcount.has_key(skip):
#						insertcount[skip]=0
#					insertcount[skip]+=1
				
			x+=1
#		printnext=False
#		for insert in insertcount.keys():
#			if float(insertcount[insert])/readdepth>options.ratio:
#				print "i", insert, insertcount[insert], readdepth
#				print line
#				printnext=True
#		for deletion in deletioncount.keys():
#			if float(deletioncount[deletion])/readdepth>options.ratio:
#				print "d", deletion, deletioncount[deletion], readdepth
#				print line
#				printnext=True
				
		
		
		#added=False
		if refqualcount[0]>options.quality and refqualcount[1]>options.strandquality and refqualcount[2]>options.strandquality and (float(refqualcount[1])/totalqual[1])>=options.ratio and (float(refqualcount[2])/totalqual[2])>=options.ratio:
			#print snplocation, refqualcount, float(refqualcount)/totalqual
			#added=True
			contigs[contig][snplocation-1]=refbase
		else:
			for snpqual in snpqualcount[0]:
				if snpqualcount[0][snpqual]>options.quality and snpqualcount[1][snpqual]>options.strandquality and snpqualcount[2][snpqual]>options.strandquality and (float(snpqualcount[1][snpqual])/totalqual[1])>=options.ratio and (float(snpqualcount[2][snpqual])/totalqual[2])>=options.ratio:
					#print snplocation, snpqualcount[snpqual], float(snpqualcount[snpqual])/totalqual, snpqual, refbase, reads, qualities
					contigs[contig][snplocation-1]=snpqual
					#added=True
					break
#		if not added:
#			print snplocation, snpqualcount, refqualcount, totalqual, refbase, snpbase


#		if readdepth>=mindepth and snpbase in ["A","C","G","T","a","c","t","g"]:
#	
#			basefoundcount=0
#			for base in words[5:9]:
#				if int(base)!=0:
#					basefoundcount=basefoundcount+1
#			if basefoundcount>1:
#				hetcount=hetcount+1
#			if graphs=='y':
#				print >> poolcovynout, 1
#		else:
#			contigs[contig][snplocation-1]=='N'


	
	
	
	
	
			
	oneseq=""
	if len(contigs)>1:
		out=open(options.output+".mfa","w")
		for contig in contigorder:
			print >> out, ">"+''.join(contig)
			print >> out, ''.join(contigs[contig])
			oneseq=oneseq+''.join(contigs[contig])
		out.close()
	else:
		for contig in contigorder:
			oneseq=oneseq+''.join(contigs[contig])
	
	
	
	out=open(options.output+".dna","w")	
	print >> out, ">"+options.output.split("/")[-1].split(".")[0]
	print >> out, oneseq



