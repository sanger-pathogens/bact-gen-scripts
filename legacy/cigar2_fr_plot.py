#!/usr/bin/env python
import string, re, gzip
import os, sys


if len(sys.argv)!=4 or '-h' in sys.argv[1:]:
	print "cigar2plot.py <reference fasta file> <cigar file (can be zipped)> <output file>"
	sys.exit()

if sys.argv[1].split('.')[-1]=='gz':
	reflen=len(''.join(gzip.open(sys.argv[1], 'r').read().split('>')[1].split('\n')[1:]))	
else:
	reflen=len(''.join(open(sys.argv[1], 'rU').read().split('>')[1].split('\n')[1:]))


forward=[0]*reflen
reverse=[0]*reflen



if sys.argv[2].split('.')[-1]=='gz':
	lines=gzip.open(sys.argv[2], 'r').readlines()
else:
	lines=open(sys.argv[2], 'r').readlines()
	
for loc, line in enumerate(lines):
	words=line.split()
	if len(words)<5:
		continue
	if words[4]=='+':
		for x in range(int(words[6])-1,int(words[7])):
			forward[x]=forward[x]+1
	elif words[4]=='-':
		for x in range(int(words[6])-1,int(words[7])):
			reverse[x]=reverse[x]+1

output=open(sys.argv[3],'w')

for x in range(0,len(forward)):
	print >> output, forward[x], reverse[x]

output.close()


#The following code plots the reads from the forward and reverse fastq files
#forward=[0]*reflen
#reverse=[0]*reflen
#for line in lines:
#	words=line.split()
#	if len(words)<5:
#		continue
#	if words[1].split('.')[-1]=='F':
#		for x in range(int(words[6])-1,int(words[7])):
#			forward[x]=forward[x]+1
#	elif words[1].split('.')[-1]=='R':
#		for x in range(int(words[6])-1,int(words[7])):
#			reverse[x]=reverse[x]+1
#
#output=open("test.plot",'w')
#
#for x in range(0,len(forward)):
#	print >> output, forward[x], reverse[x]
#
#output.close()

		