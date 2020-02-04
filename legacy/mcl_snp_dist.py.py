#!/usr/bin/env python
import string, re, gzip
import os, sys, getopt, random, math

if len(sys.argv)!=3 or '-h' in sys.argv[1:]:
	print "plot2tab.py <all.snps file> <reference file>"
	sys.exit()


snpfile=sys.argv[1]
if not os.path.isfile(snpfile):
	print "Cannot find file:", snpfile

if snpfile.split('.')[-1]=='gz':
	sfile=gzip.open(snpfile, 'r')
	lines=sfile.readlines()
else:
	lines=open(snpfile, 'rU').readlines()

reffile=sys.argv[2]
if not os.path.isfile(reffile):
	print "Cannot find file:", reffile

if reffile.split('.')[-1]=='gz':
	rfile=gzip.open(reffile, 'r')
	reflen=len(''.join(rfile.read().split("\n")[1:]))
else:
	reflen=len(''.join(open(reffile, 'rU').read().split("\n")[1:]))


maxdist=0

output=open("snp_dist.txt","w")

for x in range(len(lines)-1):
	words=lines[x].split()
	for y in range(len(lines)):
		wordsb=lines[y].split()
		if words[2]<30 or wordsb[2]<30 or x==y:
			continue
		
		if int(wordsb[3])>int(words[3]):
			if (int(wordsb[3])-int(words[3]))>maxdist:
				maxdist=(int(wordsb[3])-int(words[3]))
		else:
			if (int(words[3])-int(wordsb[3]))>maxdist:
				maxdist=(int(words[3])-int(wordsb[3]))
		
print maxdist
for x in range(len(lines)-1):
	words=lines[x].split()
	#for y in range(x,len(lines)):
	#	wordsb=lines[y].split()
	#	if words[2]<30 or wordsb[2]<30 or x==y:
	#		continue
	#	if int(wordsb[3])>int(words[3]):
			#print >> output, words[3], wordsb[3], (float(wordsb[3])-float(words[3]))
	#		print >> output, int(wordsb[3])-int(words[3]),int(wordsb[3]),int(words[3])
		#else:
			#print >> output, words[3], wordsb[3], (float(words[3])-float(wordsb[3]))
	wordsb=lines[x+1].split()
	if words[2]<30 or wordsb[2]<30:
		continue
	print >> output, int(wordsb[3])-int(words[3])

output.close()		



#os.system("~/mcl-09-261/bin/mcl snp_dist.txt  --abc -I 5 -o mclout.txt")
