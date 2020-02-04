#!/usr/bin/env python

import os, sys
import numpy
import gzip

print '\t'.join(["file", "mean mapping depth", "SDev mapping depth", "% of reference mapped"])

for f in sys.argv[1:]:

	if f.split(".")[-1]=="gz":
		os.system("gunzip "+f)
		g=".".join(f.split(".")[:-1])
	else:
		g=f

	values=[]
	linecount=0
	nonzerocount=0
	for line in open(g,"r"):
		values.append(int(line.strip()))
		if int(line.strip())>0:
			nonzerocount+=1
		linecount+=1

	print '\t'.join([str(f), str(numpy.mean(values)), str(numpy.std(values)), str((float(nonzerocount)/linecount)*100)])
	
	if f.split(".")[-1]=="gz":
		os.system("gzip "+g)