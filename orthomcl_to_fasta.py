#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math


if (len (sys.argv)!=3 and len (sys.argv)!=4) or '-h' in sys.argv[1:]:
	print "Ortho_mcl_to_fasta.py <orthomcl cluster file> <contig fasta file> <outfile prefix>"
	sys.exit()

prefix='orthomcl'
orthofile=sys.argv[1]	
contigfile=sys.argv[2]
if len(sys.argv)==4:
	prefix=sys.argv[3]

lines=open(contigfile, "rU").read().split('>')[1:]

contigs={}
straincount={}

for line in lines:
	words=line.strip().split('\n')
	contigs[words[0].split()[0]]=''.join(words[1:])
	if not straincount.has_key(words[0].split('_')[1]):
		straincount[words[0].split('_')[1]]=[0]

lines=open(orthofile, "rU").readlines()
#
#curblock='0'
#curcluster={}
#
#print 'cluster',
#strainsort=straincount.keys()
#strainsort.sort()
#for key in strainsort:
#	print key,
#print


for line in lines:
	cluster_bit=line.strip().split(":")[0]
	taxabit=line.strip().split(":")[1]
	block=cluster_bit.split("(")[0]
	genes=[]
	for gene in taxabit.split():
		genename=gene.split("(")[0]
		genes.append(genename)
	output=open(prefix+'_'+block+'.fasta','w')
	for gene in genes:
		print >> output, ">"+gene
		print >> output, contigs[gene]
	output.close()
	print block
#	
#	
#	if block!=curblock:
#		output=open(prefix+'_'+curblock+'.fasta','w')
#		keys = [ (value, key) for key, value in curcluster.iteritems() ]
#		keys.sort(reverse=True)
#		#keys.reverse()
#
#		for key in keys:
#			straincount[key[1].split('_')[1]][int(curblock)]=straincount[key[1].split('_')[1]][int(curblock)]+1
#			print >> output, '>'+key[1]
#			print >> output, contigs[key[1]]
#		output.close()
#		print curblock,
#		for key in strainsort:
#			print straincount[key][int(curblock)],
#			
#			straincount[key].append(0)
#		print
#		
#	curcluster={}
#	curblock=block
#	if len(words[1].split('_'))>3:
#		print words[1].split('_')
#		curcluster[words[1]]=int(words[1].split('_')[4])
#	else:
#		curcluster[words[1]]=int(words[1].replace('misc_feature',''))
#	
#
#output=open(prefix+'_'+curblock+'.fasta','w')
#keys = [ (value, key) for key, value in curcluster.iteritems() ]
#keys.sort()
#keys.reverse()
#
#for key in keys:
#	straincount[key[1].split('_')[1]][int(curblock)]=straincount[key[1].split('_')[1]][int(curblock)]+1
#	print >> output, '>'+key[1]
#	print >> output, contigs[key[1]]
#output.close()
#print curblock,
#for key in strainsort:
#	print straincount[key][int(curblock)],











