#!/usr/bin/env python

import os, sys, numpy

if len(sys.argv)<2 or sys.argv[1]=="-h":
	print "Usage:"
	print "plotstats.py <plot files>"
	print "Outputs mean and standard deviation of plots"

print "File Unmapped_bases %mapped Mean Std Median Max Min"

for filename in sys.argv[1:]:

	zipped=False
	if filename.split(".")[-1]=="gz":
		os.system("gunzip "+filename)
		zipped=True
		filename='.'.join(filename.split(".")[:-1])
		
	try:
		plotfile=open(filename, "rU")
	except StandardError:
		print "Cannot open file", filename
		continue
	
	data=[]
	for line in plotfile:
		if len(line)>0:
			data.append(float(line.strip()))
	
	print filename, data.count(0.0), 100-((float(data.count(0.0))/len(data))*100), numpy.mean(data), numpy.std(data), numpy.median(data), numpy.max(data), numpy.min(data)
	
	if zipped:
		os.system("gzip "+filename)
