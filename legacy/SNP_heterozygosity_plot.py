#!/usr/bin/env python
import string, re, gzip
import os, sys
import pylab
import numpy



if '-h' in sys.argv[1:]:
	print "hom_hez_ratio.py <list of all.snp files (may be zipped)>"
	sys.exit()


colours=["b", "g", "r", "c"]
legend_colours=[]

translate={"A":[7,13], "C":[8,14], "G":[9,15], "T":[10,16]}


if sys.argv[1].split('.')[-1]=='gz':
	lines=gzip.open(sys.argv[1], 'r').readlines()
else:
	lines=open(sys.argv[1], 'rU').readlines()


max_percent_array=[[],[]]

for line in lines:
	if line.split()[0]=="SNP_hez:":
		maximum=0.0
		total=0
		
		for base in line.split()[6].split("/"):
			if int(line.split()[translate[base][0]])>maximum:
				maximum=float(line.split()[translate[base][0]])
			total+=int(line.split()[translate[base][0]])
		if total>0:
			max_percent_array[0].append((maximum/total)*100)

		
		
		maximum=0.0
		total=0
		for base in line.split()[6].split("/"):
			if int(line.split()[translate[base][1]])>maximum:
				maximum=float(line.split()[translate[base][1]])
				
			total+=int(line.split()[translate[base][1]])
		if total>0:
			max_percent_array[1].append((maximum/total)*100)
		
		
		
		
	else:
		max_percent_array[0].append(100)
		max_percent_array[1].append(100)
		
		
		
legend_colours.append(pylab.Rectangle((0, 0), 1, 1, fc=colours[0]))
legend_colours.append(pylab.Rectangle((0, 0), 1, 1, fc=colours[1]))


pylab.Figure()


n, bins, patches = pylab.hist(max_percent_array, bins=20) 



pylab.legend(legend_colours, ["High quality reads","Low quality reads"])
#pylab.legend(legend_colours, ["PT2019 Multiplexed","PT2019_2 whole lane with PCR","PT2019_4 whole lane no PCR"], loc="upper left")
pylab.title("Plot of %age of reads with most common base") 
pylab.xlabel("%age of reads with most common base") 
pylab.ylabel("Frequency") 
pylab.show()
		