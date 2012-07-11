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



if (len (sys.argv)!=2) or '-h' in sys.argv[1:]:
	print "clusters_from_blast.py <fasta file>"
	sys.exit()
	
lines=open(sys.argv[1], "rU").read().split('>')[1:]

contigseqs={}
seed=lines[0].split('\n')[0].split()[0]

for line in lines:
	words=line.strip().split('\n')
	contigseqs[words[0].split()[0]]=''.join(words[1:]).replace('-','')
	


seqs=contigseqs.keys()
seqs.sort()

for seq in seqs:
	output = open(seq+'.fasta','w')
	print >> output, '>'+seq
	print >> output, contigseqs[seq]
	
	output.close()
	#print "muscle -in temp.fasta -out "+seq+".aln"
	#sys.stdout.flush()
	#os.system("muscle -in temp.fasta -out "+seq+".aln")