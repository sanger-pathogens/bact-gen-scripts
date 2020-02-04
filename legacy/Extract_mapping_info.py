#!/usr/bin/env python

#Compares mummer snp output files.  NOTE: must be from the same reference sequence



##################
# Import modules #
##################

import string, re, gzip
import os, sys, getopt, random, math
from scipy.stats import chi2


#########
# Usage #
#########

def Usage():
	os.system('clear')
	print '\mapping_vs_element_summary.py Usage:\n'

	

##############################
# Get command line arguments #
##############################


def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "hr:o:", ["ref=", "out="])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	outfile=''
	inputdirs=[]

	for opt, arg in opts:
	
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-o", "--out"):
			outfile=arg


	inputdirs=args
	
	if inputdirs==[]:
		print 'Error: No input files selected!'
		Usage()
		sys.exit()
	if outfile=='':
		outfile=ref.split('.')[0]+"_mapping_stats"
	
	
	
	return inputdirs, outfile

def mean(numberList):
	if len(numberList)==0:
		return 0
	else:
		return float(sum(numberList)) / len(numberList)
	
########
# Main #
########

	
if __name__ == "__main__":
	argv=sys.argv[1:]
	inputdirs, outfile=getOptions(argv)
	
	snpsnmapping={}
	
	for inputdir in inputdirs:
		if not os.path.isdir(inputdir):
			print "No such directory as "+inputdir+". Skipping..."
			continue
		
		elif not os.path.isfile(inputdir+"/"+inputdir+"_q60_d5_summary.out"):
			print "No such directory as "+inputdir+"/"+inputdir+"_q60_d5_summary.out. Skipping..."
			continue
		
		elif not os.path.isfile(inputdir+"/"+inputdir+"_q60_d5_summary.out"):
			print "No such directory as "+inputdir+"/"+inputdir+"_q60_d5_summary.out. Skipping..."
			continue
		
		snpsnmapping[inputdir]={}
		
		lines=open(inputdir+"/"+inputdir+"_q60_d5_summary.out", "rU").readlines()[1:-1]
		
		for line in lines:
			words=line.strip().split()
			snpsnmapping[inputdir][words[0]]=[words[1],words[15]]
		
#		for key in snpsnmapping[inputdir].keys():
#			lines=open(inputdir+"/"+key+"_ssaha/", "rU").read().split('\n')
		
	
	elements=snpsnmapping.keys()
	elements.sort()
	output=open(outfile, "w")
	print >> output, '\t'+'\t'.join(elements)
	
	strains=snpsnmapping[elements[0]].keys()
	strains.sort()
	
	
	
	for strain in strains:
	
		strainmapping=[]
		
		for element in elements:
			strainmapping.append(snpsnmapping[element][strain][0])
		
		print >> output, strain+"\t"+'\t'.join(strainmapping)
		
	
	print >> output, "\n\n"
	
	
	for strain in strains:
	
		strainmapping=[]
		
		for element in elements:
			strainmapping.append(snpsnmapping[element][strain][1])
		
		print >> output, strain+"\t"+'\t'.join(strainmapping)
		
	
	
	
	
	
	output.close()
	
	
	
	
	
	
	
	
	
	