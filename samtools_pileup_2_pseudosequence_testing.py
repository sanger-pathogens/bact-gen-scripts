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

	parser.add_option("-p", "--pileup", action="store", dest="pileup", help="pileup file (may be gzipped as long as the suffix is .gz)", default="")
	parser.add_option("-b", "--bam", action="store", dest="bam", help="bam file", default="")
	parser.add_option("-s", "--sam", action="store", dest="sam", help="sam file", default="")
	parser.add_option("-d", "--depth", action="store", dest="depth", help="Minimum number of reads matching SNP [default= %default]", default=4, type="int")
	parser.add_option("-D", "--stranddepth", action="store", dest="stranddepth", help="Minimum number of reads matching SNP per strand [default= %default]", default=2, type="int")
	#parser.add_option("-q", "--quality", action="store", dest="quality", help="Minimum base quality [default= %default]", default=120, type="int")
	parser.add_option("-R", "--RMS", action="store", dest="RMS", help="Minimum root mean squared mapping quality [default= %default]", default=25, type="int")
	#parser.add_option("-Q", "--strandquality", action="store", dest="strandquality", help="Minimum per strand base quality [default= %default]", default=60, type="int")
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
	zipped=False
	if options.pileup=="":
		DoError("No pileup file specified")
	elif not os.path.isfile(options.pileup):
		DoError("Cannot find pileup file")	
	elif options.pileup.split('.')[-1]=="gz":
		os.system("gunzip "+options.pileup)
		options.pileup='.'.join(options.pileup.split('.')[:-1])
		zipped=True
	
	if options.output=="":
		DoError("No output prefix specified")
		
	if options.sam=="" and options.bam=="":
		DoError("sam or bam file from which pileup was made must be specified")
	elif options.sam!="":
		header=os.popen(SAMTOOLS_DIR+"samtools view -S -H "+options.sam).readlines()
	else:
		header=os.popen(SAMTOOLS_DIR+"samtools view -H "+options.bam).readlines()
	
	if options.stranddepth<0:
		print "Minimum number of reads matching SNP on each strand must be >=0. Resetting to 0"
		options.stranddepth=0
	if options.depth<(options.stranddepth*2):
		print "Minimum number of reads matching SNP must be at least double that for each strand. Resetting to", options.stranddepth*2
		options.stranddepth=options.stranddepth*2
	if options.RMS<0 or options.RMS>90:
		DoError("RMS must be between 0 and 90???")
	
	contigsizes={}
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
			contigorder.append(name)
			totallength+=length
				
		
	if len(contigsizes)==0:
		DoError("No contigs found. Perhaps your sam/bam has no header?")
		
	
	contigs={}
	
	for contig in contigorder:
		contigs[contig]=["N"]*contigsizes[contig]
	
	try:
		pileupfile=open(options.pileup, "rU")
	except StandardError:
		DoError("Cannot open pileup file")
	
	mapped=0
	snps=0
	deletions=0
	insertions=0
	indels=[]
	
	count=0
	total=0.0
	hundredth=float(totallength)/100
		
			
	
	
#	printnext=False
	for line in pileupfile:
#		words=line.split()
#		if words[3].upper()in ["A", "C", "G", "T"] and int(words[4])>options.quality and (words[2].upper()==words[3].upper() or (float(words[5])/float(words[4]))>=options.ratio) and words[2]!="*":
#			contigs[words[0]][int(words[1])-1]=words[3]
		
		count=count+1
		if count>=hundredth:
			total=total+count
			count=0
			print "%.0f%% complete\r" % (100*(total/totallength)),
			sys.stdout.flush()
		
	
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
		RMS=int(words[6])
		readdepth=int(words[7])
		reads=words[8]
		qualities=words[9]
		
		if readdepth==0:
			contigs[contig][snplocation-1]='N'
			continue
		
		elif readdepth<options.depth:
			contigs[contig][snplocation-1]='N'
			continue
		
		if refbase=="*":#need to do more with these??
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

		#qualx=0
		snpbasecount=[{"A":0, "G":0, "C":0, "T":0},{"A":0, "G":0, "C":0, "T":0},{"A":0, "G":0, "C":0, "T":0}]
		#snpqualcount=[{"A":0, "G":0, "C":0, "T":0},{"A":0, "G":0, "C":0, "T":0},{"A":0, "G":0, "C":0, "T":0}]
		#insertcount={}
		#deletioncount={}
		refbasecount=[0,0,0]
		#refqualcount=[0,0,0]
		#totalqual=[0,0,0]
		totalcount=[0,0.0,0.0]
		
		#print reads, qualities
		
		x=0
		while x<len(reads):
			#print reads[x], x, qualx
			read=reads[x]
			
			if reads[x] in ['a','c','g','t',',']:
				strand=2
			else:
				strand=1
			
			read=reads[x].upper()
			
			if read in [",","."]:
				refbasecount[strand]+=1
				#refqualcount[strand]+=ord(qualities[qualx])-33
				#totalqual[strand]+=ord(qualities[qualx])-33
				totalcount[strand]+=1
				
				refbasecount[0]+=1
				#refqualcount[0]+=ord(qualities[qualx])-33
				#totalqual[0]+=ord(qualities[qualx])-33
				totalcount[0]+=1
				#qualx+=1
			elif snpbasecount[0].has_key(read):
				snpbasecount[strand][read]+=1
				#snpqualcount[strand][read]+=ord(qualities[qualx])-33
				#totalqual[strand]+=ord(qualities[qualx])-33
				totalcount[strand]+=1
				
				snpbasecount[0][read]+=1
				#snpqualcount[0][read]+=ord(qualities[qualx])-33
				#totalqual[0]+=ord(qualities[qualx])-33
				totalcount[0]+=1
				#qualx+=1
			elif read=="*":
				totalcount[0]+=1
				totalcount[0]+=0.5
				totalcount[0]+=0.5
			elif read=="^":
				x+=1
#			elif read=='N':
#				qualx+=1
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
				if read=="+" and str(snplocation)=="1418117":
					print skip, reads[x-skip:x+1]
				
				if reads[x] in ['a','c','g','t']:
					strand=2
				else:
					strand=1
				indelstring=read+reads[x-skip+1:x+1].upper()
				
				if not snpbasecount[0].has_key(indelstring):
					snpbasecount[0][indelstring]=0
					#snpqualcount[0][indelstring]=0
					snpbasecount[1][indelstring]=0
					#snpqualcount[1][indelstring]=0
					snpbasecount[2][indelstring]=0
					#snpqualcount[2][indelstring]=0
				snpbasecount[0][indelstring]+=1
				#snpqualcount[0][indelstring]+=ord(qualities[qualx-1])-33
				snpbasecount[strand][indelstring]+=1
				#snpqualcount[strand][indelstring]+=ord(qualities[qualx-1])-33
				#print snpbasecount
				
				
			x+=1

		#if Root mean squared mapping quality>threshold set by user
		if RMS>options.RMS:


			#Get indel info if we've added any 'bases' to the basecount
			#Need to do this before identifying SNPs so that deletion lines are skipped
			
			if len(snpbasecount[0])>4:
			
				for snpqual in snpbasecount[0].keys():
					if len(snpqual)>1:
						ratios=[]
						for y in range(0,len(snpbasecount)):
							if totalcount[y]==0:
								ratios.append(1.0)
							else:
								ratios.append(float(snpbasecount[y][snpqual])/totalcount[y])
						#get the next line in the file
						
						#check if the indel bases meet the thresholds for quality set in the options
						if snpbasecount[0][snpqual]>=options.depth and snpbasecount[1][snpqual]>=options.stranddepth and snpbasecount[2][snpqual]>=options.stranddepth and ratios[1]>=options.ratio and ratios[2]>=options.ratio:
											
							indelline=pileupfile.next()
							#check the next line is an indel line
							if indelline.split()[2]!="*":
								print line, indelline, pileupfile.next()
								DoError("Invalid indel line found")
							
							
							#check the indel I found is the same as that in the indel line
							if indelline.split()[8].upper()==snpqual:
								if snpqual[0]=="-":
									
									#print "Deletion of", snpqual, 'at', str(snplocation)
									for y in range(1,len(snpqual)):
										count=count+1
										delline=pileupfile.next()
										delbases=delline.split()[8]
										delnum=0
										nondelnum=0
										f=0
										while f<len(delbases):
											if delbases[f]=="*":#if a read has a deletion add 1 to delnum
												delnum+=1
											elif delbases[f]=="^":#skip any bases that are the first base in a read
												f+=2
											elif delbases[f]=="$":#skip any bases that are the last base in a read
												f+=1
											elif delbases[f].upper() in ['A', 'C', 'G', 'T', ',', '.']:#if a read has no deletion, add 1 to nondelnum
												nondelnum+=1
											f+=1
										
										if (float(delnum)/(delnum+nondelnum))>options.ratio and int(delline.split()[6])>options.RMS:
											#contigs[contig][int(delline.split()[1])-1]='-'
											mapped+=1
									deletions+=1
									indels.append([contig, str(snplocation-1), "-", snpqual[1:]])
									
								elif snpqual[0]=="+":
									#print "Insertion of", snpqual, 'after', str(snplocation-1)
									insertions+=1
									indels.append([contig, str(snplocation-1), "+", snpqual[1:]])
								else:
									print "something's gone wrong", snpqual
									sys.exit()
						



			
			ratios=[]
			for y in range(0,len(refbasecount)):
				if totalcount[y]==0:
					ratios.append(1.0)
				else:
					#ratios.append(float(refqualcount[y])/totalqual[y])
					ratios.append(float(refbasecount[y])/totalcount[y])
		
			#if bases that match the reference conform to parameters set, change the base to that of the reference
			
			if refbasecount[0]>=options.depth and refbasecount[1]>=options.stranddepth and refbasecount[2]>=options.stranddepth and ratios[1]>=options.ratio and ratios[2]>=options.ratio:
				contigs[contig][snplocation-1]=refbase
				mapped+=1
			#if not, check each possible snp to see if that meets the criterai. If it does, change the base to that of the SNP
			else:
				for snpqual in snpbasecount[0].keys():
					if len(snpqual)>1:
						continue
					ratios=[]
					for y in range(0,len(snpbasecount)):
						
						if totalcount[y]==0:
							ratios.append(1.0)
						else:
							ratios.append(float(snpbasecount[y][snpqual])/totalcount[y])
				
					if snpbasecount[0][snpqual]>=options.depth and snpbasecount[1][snpqual]>=options.stranddepth and snpbasecount[2][snpqual]>=options.stranddepth and ratios[1]>=options.ratio and ratios[2]>=options.ratio:
						contigs[contig][snplocation-1]=snpqual
						mapped+=1
						snps+=1
						
						break
			
			

	
	
	if zipped:
		os.system("gzip "+options.pileup)
	
	
			
	out=open(options.output+".mfa","w")
	for contig in contigorder:
		print >> out, ">"+''.join(contig)
		print >> out, ''.join(contigs[contig])
	out.close()
	
	
	print 'Total bases in reference:', totallength
	print 'Mapped:%d (%.2f%%)' % (mapped, (float(mapped)/totallength)*100)
	print 'SNPs:', snps
	print 'Insertions:', insertions
	print 'Deletions:', deletions
	out=open(options.output+"_indels.txt","w")
	for indel in indels:
		print >> out, "\t".join(indel)
	out.close()

