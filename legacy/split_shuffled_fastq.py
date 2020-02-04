#!/usr/bin/env python
import os, sys

if len(sys.argv)!=3 or "-h" in sys.argv:
	print "split_shuffled_fastq.py <input shuffled fastq> <output prefix>"
	sys.exit()

lines=open(sys.argv[1],"rU").readlines()

outfile1=open(sys.argv[2]+"_1.fastq","w")
outfile2=open(sys.argv[2]+"_2.fastq","w")
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
