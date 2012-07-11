#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math


if (len (sys.argv)!=4) or '-h' in sys.argv[1:]:
	print "strain_snps.py <fasta file 1> <fasta file 2> <list of pairs>"
	sys.exit()

data=open(sys.argv[1], "rU").read()[1:].split('\n>')

seqs1={}

for datum in data:
	if len(datum.strip())==0:
		continue
	name=datum.strip().split('\n')[0].split()[0]
	seqs1[name]=''.join(datum.strip().split('\n')[1:])

data=open(sys.argv[2], "rU").read().split('>')[1:]

seqs2={}

for datum in data:
	if len(datum.strip())==0:
		continue
	name=datum.strip().split('\n')[0].split()[0]
	seqs2[name]=''.join(datum.strip().split('\n')[1:])
	
lines=open(sys.argv[3], "rU").readlines()

kaksin=open("kaks.input", "w")
missing=open("missing_seqs.txt", "w")
for line in lines:
	if len(datum.strip())==0:
		continue
	words=line.strip().split()
	
	somemissing='n'
	if not words[0] in seqs1.keys():
		print >> missing, words[0]
		somemissing='y'
	if not words[1] in seqs2.keys():
		print >> missing, words[1]
		somemissing='y'
		
	if somemissing=='y':
		continue
	
	musclein=open("temp.fna", "w")
	print >> musclein, ">"+words[0]
	print >> musclein, seqs1[words[0]]
	print >> musclein, ">"+words[1]
	print >> musclein, seqs2[words[1]]
	musclein.close()
	
	os.system("~sh16/scripts/Nuc_to_aa.py -i temp.fna -o temp.faa")
	
	
	os.system("muscle -in temp.faa -out temp.muscleout")
	
	os.system("~sh16/scripts/Protein_align_to_nuc.py -n temp.fna -p temp.muscleout -o temp.nucaln")
	
	
	data=open("temp.nucaln", "rU").read().split('>')[1:]
	
	seq1=''.join(data[0].strip().split('\n')[1:])
	seq2=''.join(data[1].strip().split('\n')[1:])
	
	print >> kaksin, words[0]+","+words[1]
	print >> kaksin, seq1
	print >> kaksin, seq2+"\n"
	

missing.close()
kaksin.close()

os.system("~/KaKs_Calculator1.2/src/KaKs_Calculator -m NG -i kaks.input -o kaks.out")
os.system("rm -f temp.*")

