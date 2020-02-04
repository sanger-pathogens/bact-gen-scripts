#!/usr/bin/env python
import string, re, gzip
import os, sys


if '-h' in sys.argv[1:]:
	print "hom_hez_ratio.py <list of all.snp files (may be zipped)>"
	sys.exit()


translate={"A":7, "G":8, "C":9, "T":10}

print "Name\tHomozygous\tHeterozygous\tHet/Hom\tError_rate\to-e/e"

for filename in sys.argv[1:]:

	if sys.argv[1].split('.')[-1]=='gz':
		lines=gzip.open(filename, 'r').readlines()
	else:
		lines=open(filename, 'rU').readlines()

	homcount=0
	hezcount=0
	SNPcount=0
	othercount=0
	
	for line in lines:
		if line.split()[0]=="SNP_hom:":
			homcount=homcount+1
			SNPcount=SNPcount+int(line.split()[translate[line.split()[6]]])
		elif line.split()[0]=="SNP_hez:":
			hezcount=hezcount+1
			maximum=0
			total=0
			for base in line.split()[6].split("/"):
				if int(line.split()[translate[base]])>maximum:
					maximum=int(line.split()[translate[base]])
				total=total+int(line.split()[translate[base]])
			SNPcount=SNPcount+maximum
			othercount=othercount+(total-maximum)
	
	if homcount>0:
		hethez=str(float(hezcount)/homcount)
	else:
		hethez="inf"
	
	if SNPcount>0:
		ratio=str(float(othercount)/SNPcount)
		oe=str(((float(othercount)/SNPcount)-0.01)/0.01)
	else:
		ratio="inf"
		oe="-"
	
	print filename+"\t"+str(homcount)+"\t"+str(hezcount)+"\t"+hethez+"\t"+ratio+"\t"+oe
		