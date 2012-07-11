#!/usr/bin/env python
import string, re, gzip
import os, sys, getopt, random, math


if len (sys.argv)!=2 or '-h' in sys.argv[1:]:
	print "Cluster_by_alignment.py <multifasta file of sequences>"
	sys.exit()

lines=open(sys.argv[1],'rU').read().split('>')[1:]

sequences={}
names=[]
for line in lines:
	sequences[line.split('\n')[0].split()[0]]=''.join(line.split('\n')[1:])
	names.append(line.split('\n')[0].split()[0])

comparisons={}

for name in names:
	comparisons[name]={}
	for nameb in names:
		if name!=nameb:
			comparisons[name][nameb]=0

for x, name in enumerate(names):	
	output=open('seq1.fasta','w')
	print >> output, ">"+name
	print >> output, sequences[name]
	output.close()
	for y in range(x+1,len(names)):
		output=open('seq2.fasta','w')
		print >> output, ">"+names[y]
		print >> output, sequences[names[y]]
		output.close()
		os.system("progressiveMauve --output=test --backbone-output=test.bb seq1.fasta seq2.fasta")
		if os.path.isfile('test.bb'):
			bbs=open('test.bb','rU').readlines()[1:]
			
			seq1=0
			seq2=0
			for bb in bbs:
				words=bb.split()
				if int(words[1].replace('-',''))-int(words[0].replace('-',''))!=0 and int(words[3].replace('-',''))-int(words[2].replace('-',''))!=0:
					seq1=seq1+(int(words[1].replace('-',''))-int(words[0].replace('-','')))
					seq2=seq2+(int(words[3].replace('-',''))-int(words[2].replace('-','')))
			
			os.system('rm -f test.bb *.sml')
			
			comparisons[name][names[y]]=(float(seq1)/len(sequences[name]))*100
			comparisons[names[y]][name]=(float(seq2)/len(sequences[names[y]]))*100
		else:
			comparisons[name][names[y]]=0.0
			comparisons[names[y]][name]=0.0

output = open("testing.csv","w")
print >> output, '\t'+'\t'.join(names),

for name in names:
	print >> output, '\n'+name,
	for nameb in names:
		if name!=nameb:
			print >> output, '\t'+str(comparisons[name][nameb]),
		else:
			print >> output, '\t-',
		


output.close()