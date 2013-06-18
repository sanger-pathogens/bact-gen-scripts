#!/usr/bin/env python

#Compares mummer snp output files.  NOTE: must be from the same reference sequence



#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------

import string, re
import os, sys, getopt, random, math


def Usage():
	print 'Pairwise_SNPs.py Usage:'
	print 'Pairwise_SNPs.py -i [input SNP file] -o [output file]'
	print 'Copyright Simon R Harris, Wellcome Trust Sanger Institute, UK. 2008'

	

#------------------------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------------------------


def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "hi:o:", ["help", "in=", "out="])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	outfile=''
	inputfile=''

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-i", "--in"):
			inputfile=arg
		elif opt in ("-o", "--out"):
			outfile=arg
	


	if inputfile=='':
		print 'Error: No input file selected!'
		Usage()
		sys.exit()
	if outfile=='':
		outfile=inputfile+'.snps'

	return inputfile, outfile



if __name__ == "__main__":
	argv=sys.argv[1:]
	inputfile, outfile=getOptions(argv)

	lines=open(inputfile, 'rU').readlines()
	
	refpos=-1
	synpos=-1
	nostrains=0
	words=lines[0].split()
	for x, word in enumerate(words):
		if word=="Synonymous/Non-synonymous":
			synpos=x
		elif word=="Ref_base":
			refpos=x
			
	names=[]
	incref='n'
	if words[0]!='Position_in_alignment':
		names.append(words[0].replace('Position_in_',''))
		incref='y'
		
	for name in words[refpos+3:]:
			names.append(name)
			
	
	nostrains=len(names)
	
	if refpos==-1:
		print "Error: No reference column found"
		sys.exit()
	
	pairsnps=[[0]*nostrains for n in xrange(nostrains) ]
	synpairsnps=[[0]*nostrains for n in xrange(nostrains) ]
	intergenicsnps=[[0]*nostrains for n in xrange(nostrains) ]
	gapcount=[[0]*nostrains for n in xrange(nostrains) ]
	
	for line in lines[1:]:
		words=line.split()
		
		site=[]
		if incref=='y':
			site.append(words[refpos])
		
		for x in words[refpos+3:]:
			site.append(x)
		
			
		for x in range(len(site)):
			for y in range(x,len(site)):
				if x!=y:
					#print x,y, site[x],site[y], pairsnps[x][y]
					if site[x]!=site[y] and site[x]!='-' and site[y]!='-':
						pairsnps[x][y]=pairsnps[x][y]+1
						if words[synpos].replace('/1','').replace('1/','')=='N':
							synpairsnps[x][y]=synpairsnps[x][y]+1
						elif words[synpos]=='-':
							intergenicsnps[x][y]=intergenicsnps[x][y]+1
							
					if site[x]!='-' and site[y]!='-':
						gapcount[x][y]=gapcount[x][y]+1
					#print pairsnps[x][y], pairsnps[x]
	
	
#	
#	output=open(outfile,'w')
#	print >> output, ' ',
#	for name in names:
#		print >> output, name,
#	print >> output
#	for x in range(nostrains):
#		print >> output, names[x],
#		for y in range(nostrains):
#			if x<y:
#				print >> output, pairsnps[x][y],
#			elif y<x:
#				if synpairsnps[y][x]>0:
#					print >> output, float(pairsnps[y][x])/gapcount[y][x]*100,
#				else:
#					print >> output, 0,
#			else:
#				print >> output, '-',
#			
#		print >> output
#		
#		
#	print >> output, '\n ',
#	for name in names:
#		print >> output, name,
#	print >> output
#	for x in range(nostrains):
#		print >> output, names[x],
#		for y in range(nostrains):
#			if x<y:
#				print >> output, synpairsnps[x][y],
#			elif y<x:
#				if synpairsnps[y][x]>0:
#					print >> output, float(synpairsnps[y][x])/gapcount[y][x]*100,
#				else:
#					print >> output, 0,
#			else:
#				print >> output, '-',
#			
#		print >> output
#	output.close()
	
	
	
	output=open(outfile,'w')
	for x in range(nostrains):
		for y in range(x, nostrains):
			if x!=y:
				print >> output, names[x],
				print >> output, names[y],
				if synpos==-1:
					print >> output, pairsnps[x][y]
				else:
					print >> output, pairsnps[x][y],synpairsnps[x][y],intergenicsnps[x][y],
					
					if synpairsnps[x][y]>0:
						print >> output, float(pairsnps[x][y])/gapcount[x][y]*100
					else:
						print >> output, 0

	output.close()
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
