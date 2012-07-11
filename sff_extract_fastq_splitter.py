#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math

def revcomp(sequence):
	rev=sequence[::-1]
	revcomp=''
	d={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'g':'c', 'c':'g'}
	for i in rev:
		if d.has_key(i):
			revcomp=revcomp+d[i]
		else:
			revcomp=revcomp+i
	
	return revcomp
	
def rev(sequence):
	rev=sequence[::-1]
	
	return rev

if (len(sys.argv)!=2):
	print "sff_extract_fastq_splitter.py <fastq file from sff_extract>"
	sys.exit()

lines=open(sys.argv[1], 'rU').readlines()

single=open(sys.argv[1].split('.')[0]+'.single.tmp', 'w')
forward=open(sys.argv[1].split('.')[0]+'.forward.tmp', 'w')
reverse=open(sys.argv[1].split('.')[0]+'.reverse.tmp', 'w')
paired=open(sys.argv[1].split('.')[0]+'.paired.tmp', 'w')


for x in range(0,len(lines),4):
	if lines[x].strip().split('.')[-1]=='fn':
		print >> single, lines[x].strip().replace('.fn','')
		print >> single, revcomp(lines[x+1].strip())
		print >> single, lines[x+2].strip().replace('.fn','')
		print >> single, rev(lines[x+3].strip())
		
		print >> forward, lines[x].strip().replace('.fn','')
		print >> forward, revcomp(lines[x+1].strip())
		print >> forward, lines[x+2].strip().replace('.fn','')
		print >> forward, rev(lines[x+3].strip())
		
		print >> reverse, lines[x].strip().replace('.fn','')
		print >> reverse, "N"*len(lines[x+1].strip())
		print >> reverse, lines[x+2].strip().replace('.fn','')
		print >> reverse, "&"*len(lines[x+3].strip())
		
	if lines[x].strip().split('.')[-1]=='part1':
		if lines[x+4].strip().split('.')[-1]=='part2' and lines[x+4].strip().split('.')[0]==lines[x].strip().split('.')[0]:
			print >> paired, lines[x].strip().replace('.part1','/1')
			print >> paired, revcomp(lines[x+1].strip())
			print >> paired, lines[x+2].strip().replace('.part1','/1')
			print >> paired, rev(lines[x+3].strip())
			print >> paired, lines[x+4].strip().replace('.part2','/2')
			print >> paired, lines[x+5].strip()
			print >> paired, lines[x+6].strip().replace('.part2','/2')
			print >> paired, lines[x+7].strip()
			print >> forward, lines[x].strip().replace('.part1','/1')
			
			if len(lines[x+1].strip())<len(lines[x+5].strip()):
				print >> forward, revcomp(lines[x+1].strip())+"N"*(len(lines[x+5].strip())-len(lines[x+1].strip()))
				print >> forward, lines[x+2].strip().replace('.part1','/1')
				print >> forward, revcomp(lines[x+3].strip())+"&"*(len(lines[x+5].strip())-len(lines[x+1].strip()))
			else:
				print >> forward, revcomp(lines[x+1].strip())
				print >> forward, lines[x+2].strip().replace('.part1','/1')
				print >> forward, rev(lines[x+3].strip())
			
			print >> reverse, lines[x+4].strip().replace('.part2','/2')
			if len(lines[x+1].strip())>len(lines[x+5].strip()):
				print >> reverse, lines[x+5].strip()+"N"*(len(lines[x+1].strip())-len(lines[x+5].strip()))
				print >> reverse, lines[x+6].strip().replace('.part2','/2')
				print >> reverse, lines[x+7].strip()+"&"*(len(lines[x+1].strip())-len(lines[x+5].strip()))
			else:
				print >> reverse, lines[x+5].strip()
				print >> reverse, lines[x+6].strip().replace('.part2','/2')
				print >> reverse, lines[x+7].strip()
		else:
			print >> single, lines[x].strip().replace('.part1','')
			print >> single, revcomp(lines[x+1].strip())
			print >> single, lines[x+2].strip().replace('.part1','')
			print >> single, rev(lines[x+3].strip())
		
			print >> forward, lines[x].strip().replace('.part1','')
			print >> forward, revcomp(lines[x+1].strip())
			print >> forward, lines[x+2].strip().replace('.part1','')
			print >> forward, rev(lines[x+3].strip())
			
			print >> reverse, lines[x].strip().replace('.part1','')
			print >> reverse, "N"*len(lines[x+1].strip())
			print >> reverse, lines[x+2].strip().replace('.part1','')
			print >> reverse, "&"*len(lines[x+3].strip())

forward.close()
reverse.close()
#paired.close()
single.close()