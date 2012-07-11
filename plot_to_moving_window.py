#!/usr/bin/env python

#Compares mummer snp output files.  NOTE: must be from the same reference sequence



#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------

import string, re
import os, sys, getopt, random, math
from scipy.stats import chi2



def Usage():
	print 'plot_to_moving_window.py Usage:'
	print 'plot_to_moving_window.py [options] <input alignment file>'
	print 'Options:'
	print '-w\t\twindow size'
	print '-s\t\tstep size'
	print '-r\t\treference tree'
	print '-f\n\nforce (overwrites previous moving_tree files)'


#------------------------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------------------------


def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "hw:s:fr:", ["window=", "step=", "force=", "ref="])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	window=1000
	step=500
	force='n'
	ref=''
	
	for opt, arg in opts:
	
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-w", "--window"):
			window=int(arg)
		elif opt in ("-s", "--step"):
			step=int(arg)
	
	inputfile=args[0]
	
	
	if inputfile==[]:
		print 'Error: No input file selected!'
		Usage()
		sys.exit()

	return window, step, inputfile
	


	
if __name__ == "__main__":
	argv=sys.argv[1:]
	window, step, inputfile=getOptions(argv)
	
	
	lines=open(inputfile,'rU').readlines()
	
	values=[]
	for x, v in enumerate(lines):
		if v>0:
			values.append([x,float(v)])
	
	
	output=open('.'.join(inputfile.split('.')[:-1])+"_moving_window.plot", "w")
	
	lastbase=-1
	for base in range(0,len(values),int(step)):
		start=base-(window/2)
		end=base+(window/2)
		total=0.0
		
		if start>=0 and end<len(values):	
			for f in values[start:end]:
				total=total+f[1]
			total=total/((window/2)*2)
		elif start<0 and end<len(values):
			for f in values[len(values)+start:]:
				total=total+f[1]
			for f in values[:end]:
				total=total+f[1]
			total=total/((window/2)*2)
		elif start>=0 and end>=len(values):
			for f in values[start:]:
				total=total+f[1]
			for f in values[:end-len(values)]:
				total=total+f[1]
			total=total/((window/2)*2)
		else:
			print "Error! Your window size is larger than your number of snp sites!"
			sys.exit()
	
	
		while lastbase<values[base][0]:
			print >> output, total
			lastbase=lastbase+1
	
	
	
	
	
	
	
	
	