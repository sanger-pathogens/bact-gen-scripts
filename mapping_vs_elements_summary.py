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
	
	ref=''
	outfile=''
	inputdirs=[]

	for opt, arg in opts:
	
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-r", "--ref"):
			ref=arg
		elif opt in ("-o", "--out"):
			outfile=arg


	inputdirs=args
	
	if inputdirs==[]:
		print 'Error: No input files selected!'
		Usage()
		sys.exit()
	if outfile=='':
		outfile=ref.split('.')[0]+"_mapping_stats"
	


	if ref=='':
		print 'Error: No reference dna file selected!'
		Usage()
		sys.exit()
	
	
	return ref, inputdirs, outfile

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
	ref, inputdirs, outfile=getOptions(argv)
	
	
	sequs=open(ref, 'rU').read().split(">")[1:]
	elements=[]
	lengths=[]
	for seq in sequs:
		elements.append(seq.split('\n')[0].split()[0])
		lengths.append(len(''.join(seq.replace(' ','').strip().split('\n')[1:])))
	
	output=open(outfile, 'w')
	
	print >> output, '\t'+'\t'.join(elements)
	
	for inputdir in inputdirs:
		name=inputdir.split('_')[1]
		data=gzip.open(inputdir+"/allcoverage.plot.gz", "r").read().strip().split('\n')
		data=map(int, data)
		
		tabout=open(name+'_mapping.tab', 'w')
		print >> tabout, 'ID   mapping'
		
		start=0
		
		percents=[]
		print name
		for elnum, element in enumerate(elements):
			end=start+lengths[elnum]
			zeros=data[start:end].count(0)
			if end-start<lengths[elnum]:
				zeros=zeros+(lengths[elnum]-(end-start))
			nonzeros=lengths[elnum]-zeros
			average=mean(data[start:end])
			percentmapped=(float(nonzeros)/lengths[elnum])*100
			percents.append(percentmapped)
			#print element, average, zeros, nonzeros, percentmapped
			
			inblock='n'
			unmapstart=0
			for x, base in enumerate(data[start:end]):
				y=x+start+1
				
				if base==0 and inblock=='n':
					unmapstart=y
					inblock='y'
				elif inblock=='y':
					print >> tabout, 'FT   misc_feature    '+str(unmapstart)+'..'+str(y-1)
					inblock='n'
			if inblock=='y':
				print >> tabout, 'FT   misc_feature    '+str(unmapstart)+'..'+str(y-1)
					
	
			start=end
		print >> output, name+'\t'+'\t'.join(map(str, percents))
		tabout.close()
	
	output.close()
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
        
