#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math


if (len (sys.argv)!=3 and len (sys.argv)!=4) or '-h' in sys.argv[1:]:
	print "MUMmer_tiling_to_tab.py <mcl cluster file> <contig fasta file> <outfile prefix>"
	sys.exit()

prefix='tribemcl'
tribefile=sys.argv[1]	
contigfile=sys.argv[2]
if len(sys.argv)==4:
	prefix=sys.argv[3]

lines=open(contigfile, "rU").read().split('>')[1:]

contigs={}
straincount={}

for line in lines:
	words=line.strip().split('\n')
	contigs[words[0].split()[0]]=''.join(words[1:])
	if not straincount.has_key('_'.join(words[0].split('_')[:2])):
		straincount['_'.join(words[0].split('_')[:2])]=[0]

lines=open(tribefile, "rU").readlines()
outputb=open("lengths_in_clusters.txt", 'w')

curblock='0'

print 'cluster',
print >>outputb, 'cluster',
strainsort=straincount.keys()
strainsort.sort()
for key in strainsort:
	print key,
	print >> outputb, key,
print
print >> outputb, '\n',

for curblock, line in enumerate(lines):
	words=line.strip().split()
	os.system("mkdir cluster_"+str(curblock))
	print "cluster", curblock
	curcluster={}
	
	lengths={}
	
	for word in words:
		output=open("cluster_"+str(curblock)+"/"+"_".join(word.split('_')[:2])+".fasta", 'a')
		print >> output, '>'+word
		print >> output, contigs[word]
		if not lengths.has_key('_'.join(word.split('_')[:2])):
			lengths['_'.join(word.split('_')[:2])]=len(contigs[word])
		else:
			lengths['_'.join(word.split('_')[:2])]=lengths['_'.join(word.split('_')[:2])]+len(contigs[word])
		output.close()
	
	
	strainsort.sort()
	
	print >>outputb, curblock,
	for name in strainsort:
		if lengths.has_key(name):
			print >>outputb, lengths[name],
		else:
			print >> outputb, '0',
			
	print >> outputb, '\n',
		
outputb.close()		
		
		
		
		
		
		
		
		
		
		