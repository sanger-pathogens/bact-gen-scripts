#!/usr/bin/env python
import os, sys

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
	
if len(sys.argv)!=3 or "-h" in sys.argv:
	print "split_pacbio_fastq.py <input pacbio fastq> <output prefix>"
	sys.exit()



outfile1=open(sys.argv[2]+"_1.fastq","w")
outfile2=open(sys.argv[2]+"_2.fastq","w")

readfile=open(sys.argv[1],"rU")

for line in readfile:
	newlines=[line.strip()]
	for y in range(0,3):
		newlines.append(readfile.next().strip())
	
	
	#print newlines[1]
	if len(newlines[1])>500:
	
		print >> outfile1, newlines[0]+":"+str(len(newlines[1]))+"/1"
		print >> outfile2, newlines[0]+":"+str(len(newlines[1]))+"/2"
		
		print >> outfile1, newlines[1][:201]
		print >> outfile2, revcomp(newlines[1])[:201]
		
		print >> outfile1, newlines[2]
		print >> outfile2, newlines[2]
		
		print >> outfile1, newlines[3][:201]
		print >> outfile2, newlines[3][::-1][:201]

outfile1.close()
outfile2.close()
