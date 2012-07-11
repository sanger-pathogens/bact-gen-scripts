#!/usr/bin/env python
import string, re, gzip
import os, sys, getopt, random, math

def Usage():
	print 'Protein_align_to_nuc.py Usage:'
	print 'Protein_align_to_nuc.py -n [input nucleotide fasta file] -p [input protein alignment] -o [output nucleotide alignment file name]'
	print 'Written by: Simon R Harris, Wellcome Trust Sanger Institute, UK. 2009'

	

#------------------------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------------------------


def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "hf:", ["help", "fasta="])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	fasta=''
	crunchies=[]

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-f", "--fasta"):
			fasta=arg

	crunchies=args
	


	if fasta=='':
		print 'Error: Missing input file!'
		Usage()
		sys.exit()

	return fasta, crunchies



if __name__ == "__main__":
        argv=sys.argv[1:]
        fasta, crunchies=getOptions(argv)

lines=open(fasta,'rU').read().split('>')[1:]

elements={}
names=[]

output=open("plasmidmappingpercents.csv", "w")

for line in lines:
	elements[line.split('\n')[0].split()[0]]=''.join(line.split('\n')[1:])
	names.append(line.split('\n')[0].split()[0])


print >> output, 'isolates,'+','.join(names)

for crunch in crunchies:
	lines=open(crunch,'rU').readlines()
	mapping={}
	
	print crunch
	
	for name in names:
		mapping[name]=[0]*len(elements[name])
		#print name, len(elements[name])
	
	for line in lines:
		words=line.strip().split()
		if int(words[0])<100:
			continue
		#print words[4]
		for x in range(int(words[2])-1,int(words[3])-1):
			#print x, mapping[words[4]][x], int(words[1])
			if mapping[words[4]][x]<int(words[1]):
				mapping[words[4]][x]=int(words[1])
			#print x, mapping[words[4]][x], int(words[1])
		
	
	outstring=crunch
	for name in names:
		total=0
		number=0
		for x in mapping[name]:
			#print x, name
			if x>0:
				number=number+1
				total=total+x
		
		percent=float(total)/len(elements[name])
		
		outstring=outstring+','+str(percent)
	print >> output, outstring

output.close()
		