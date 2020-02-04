#!/usr/bin/env python
import string, re, gzip
import os, sys
import pylab
import numpy

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

if (len(sys.argv)!=4 and len(sys.argv)!=5 and len(sys.argv)!=6 and len(sys.argv)!=7) or '-h' in sys.argv[1:]:
	print "kmer_counter.py <reference fasta file (can be zipped)> <read length> <is genome linear/circular (l/c)> <step>"
	sys.exit()

if sys.argv[1].split('.')[-1]=='gz':
	ref=''.join(gzip.open(sys.argv[1], 'r').read().split('>')[1].split('\n')[1:]).upper()
else:
	ref=''.join(open(sys.argv[1], 'rU').read().split('>')[1].split('\n')[1:]).upper()

try:
	readlen=int(sys.argv[2])
	if readlen>61 or readlen<1:
		print 'kmer must be an integer between 1 and 61'
		sys.exit()
except ValueError:
	print 'kmer must be an integer'
	sys.exit()

step=1
if len(sys.argv)==5 or len(sys.argv)==7:
	try:
		step=int(sys.argv[4])
	except ValueError:
		print 'step must be an integer'
		sys.exit()
lc=sys.argv[3]

if lc not in ['l','c']:
	print 'linear/circular option must be either "l" or "c"'
	sys.exit()



repeatcounts=[]
kmerlist=[]
for kmerlength in range(15,readlen+1,2):
	count=0
	kmers={}
	
	for x in range(0,len(ref),step):
		if x+kmerlength>len(ref) and lc=="l":
			continue
		count=count+1
			
		if x+kmerlength<len(ref):
			
			dna=ref[x:x+kmerlength]
		else:
			dna=ref[x:]+ref[:(x+kmerlength)-len(ref)]
			
		if kmers.has_key(dna):
			kmers[dna]=kmers[dna]+1
		else:
			kmers[dna]=1	
	
	kmerkeys=kmers.keys()
	kmerkeys.sort()
	
	repeatcount=0
	for key in kmerkeys:
		if kmers[key]>1:
			#print key, kmers[key]
			repeatcount=repeatcount+1
	
	print kmerlength, repeatcount, len(kmers.keys()), count
	repeatcounts.append((float(repeatcount)/len(kmers.keys())*100))
	kmerlist.append(kmerlength)

kmers={}
kmerkeys=[]

print kmerlist, repeatcounts
pylab.Figure()
pylab.plot(kmerlist,repeatcounts)
pylab.title("Kmer Repeat Counts") 
pylab.xlabel("Kmer length") 
pylab.ylabel("Percentage of kmers that are repeated") 
pylab.show()
	
	
		
		