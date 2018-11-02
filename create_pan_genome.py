#!/usr/bin/env python

#Compares mummer snp output files.  NOTE: must be from the same reference sequence



##################
# Import modules #
##################

import string, re
import os, sys, getopt, random, math, time
from random import *


#if len(sys.argv)!=3:
#	print "Usage: create_pan_genome.py <ssaha_folders>"
#	sys.exit()
	

chars = string.ascii_letters + string.digits

tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))

reffile=tmpname+"_ref.mfa"

folder=sys.argv[1]
os.system("gunzip "+folder+"/unmap.fastq.gz")

outfile=open(tmpname+".fastq","w")
lines=open(folder+"/unmap.fastq", "rU").readlines()
x=0
while x<len(lines):
	newlines=[lines[x].strip()]
	x+=1
	for y in range(0,3):
		newlines.append(lines[x].strip())
		x+=1
	
	if x==len(lines) or lines[x].strip()[:-2]!=newlines[0].strip()[:-2]:
		continue
		
	for y in range(0,4):
		newlines.append(lines[x].strip())
		x+=1
	
	for outline in newlines:
		print >> outfile, outline

outfile.close()

os.system("gzip "+folder+"/unmap.fastq")


os.system('/nfs/pathogen/sh16_scripts/velvet_assembly.sh -n -e 15 -o "-min_contig_lgth 500" -p -f '+tmpname+".fastq")

os.system("cp "+tmpname+"_velvet/contigs.fa "+reffile)

for folder in sys.argv[2:]:
	
	os.system("gunzip "+folder+"/unmap.fastq.gz")
	
	outfile1=open(tmpname+"_1.fastq","w")
	outfile2=open(tmpname+"_2.fastq","w")
	lines=open(folder+"/unmap.fastq", "rU").readlines()
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
	#os.system("gunzip "+folder+"/unmap.fastq")
	
	newref=open(tmpname+".dna", "w")
	print >> newref, ">Reference"
	newref.close()
	os.system('grep -v "^>" '+reffile+" >> "+tmpname+".dna")
	
	
	os.system("/nfs/pathogen/sh16_scripts/run_multiple_mappings_auto.py -r "+tmpname+".dna -M -I -P ssaha -L -p "+tmpname+"_[12].fastq")
	
	os.system("gunzip "+tmpname+"_ssaha/unmap.fastq.gz")
	os.system('/nfs/pathogen/sh16_scripts/velvet_assembly.sh -n -e 15 -o "-min_contig_lgth 500" -p -f '+tmpname+"_ssaha/unmap.fastq")
	os.system("gzip "+tmpname+"_ssaha/unmap.fastq")
	os.system("cat "+tmpname+"_ssaha/unmap_velvet/contigs.fa >> "+reffile)
	os.system("rm -rf "+tmpname+"_ssaha 1tmp*_sbs.sh")

os.system("mv "+reffile+" accessory_genome.mfa")	
os.system("rm -rf "+tmpname+"*")
	