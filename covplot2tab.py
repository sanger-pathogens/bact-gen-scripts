#!/usr/bin/env python
import string, re, gzip
import os, sys, getopt, random, math

for x in sys.argv[1:]:
	print x
	if x.split('.')[-1]=='gz':
		lines=gzip.open(x, 'r').readlines()
	else:
		lines=open(x, 'rU').readlines()
	output=open(x.split('.')[0]+'.tab', 'w')
	print >> output, "ID   coverage"
	count=1
	inblock='n'
	featureline=''
	for line in lines:
		cov=int(line.strip().split()[0])
		if cov==0 and inblock=='n':
			featureline="FT   misc_feature    "+str(count)+".."
			inblock='y'
		elif cov!=0 and inblock=='y':
			featureline=featureline+str(count-1)
			print >> output, featureline
			inblock='n'
		
		count=count+1
	if inblock=='y':
		featureline=featureline+str(count-1)
		print >> output, featureline

	output.close()
	