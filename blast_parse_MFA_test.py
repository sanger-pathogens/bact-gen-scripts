#!/usr/bin/env python
# !/usr/bin/env python
import string
import os, sys
from Bio import Seq
from Bio import AlignIO
from optparse import OptionParser, OptionGroup
import tre

sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_general import *
from Si_SeqIO import *
from Si_SNPs_temp import *
#from  multiprocessing import cpu_count



import time



		#os.system("blastall -e 100 -p blastn -a "+str(options.processors)+" -d /lustre/scratch103/sanger/sh16/Chlamydia/primerdesigns/softmasked_dusted.fa -i "+tmpname+"primerseqs.fasta -o "+tmpname+"blast.out -m 3 -W 7 -F F")
print "\nParsing blast results"
sys.stdout.flush()

currseq=""
primers={}
readfile=open(sys.argv[1],"rU")
for line in readfile:
	line.strip()
	if len(line.split())==2 and line.split()[0]=="Query=":
		currseq=line.split()[1]
		try:
			currseqlen=int(readfile.next().strip().split()[0].replace("(",""))
		except StopIteration:
			break
		primers[currseq]=1000
		
		
		while not ( len(line.split())==4 and len(line.split()[0])>2 and line.split()[0][-2:]=="_0") and not (len(line.split())>0 and line.split()[0]=="BLASTN"):
			try:
				stop=False
				line=readfile.next().strip()
			except StopIteration:
				stop=True
				break
				
		if stop or (len(line.split())>0 and line.split()[0]=="BLASTN"):
			continue
		
		start=string.find(line, line.split()[2])
		end=start+len(line.split()[2])
		bestmatch=4
		missingcbit=currseqlen-int(line.strip().split()[3])
		line=readfile.next().strip()
		
		
		while len(line.split())==4 and line.split()[0]!="BLASTN":
			matches=0
			cmatches=0
			score=0.0
			matchlen=len(line[end:start:-1].strip())
			matchpos=0
			started=False
			#print matches, cmatches, score, matchpos, started, matchlen
			for x, hit in enumerate(line[end:start:-1]):
				if started:
					matchpos+=1
				elif hit!=" ":
					started=True
				if hit==".":
					matches+=1
					if x+missingcbit<5:
						cmatches+=1
					score+=1.0
				elif started and hit==" " and matchpos<=matchlen:
					score-=2
				elif started and hit=="N" and matchpos<=matchlen:
					score-=0.25
				elif started and matchpos<=matchlen:
					score-=1.0
			#print matches, cmatches, score, matchpos, started, matchlen
#			if ((5-cmatches)<2 and (currseqlen-score)<6) or (currseqlen-score)<5:
#				primers[currseq]["HUMAN_MIN_MISPRIMING_PASS"]=False

			if (currseqlen-score)<primers[currseq]:	
				primers[currseq]=(currseqlen-score)
				
			line=readfile.next().strip()
		
count=0
badcount=0
for primername in primers.keys():
	if primers[primername]==1000:
		print primername
		badcount+=1
	count+=1
	
	
print count, badcount

#os.system("formatdb -i "+tmpname+"primerseqs.fasta -p F")


