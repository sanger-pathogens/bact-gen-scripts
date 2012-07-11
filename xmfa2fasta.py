#!/usr/bin/env python

#Compares mummer snp output files.  NOTE: must be from the same reference sequence



#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------

import string, re, math
import os, sys, getopt


def Usage():
	print "-r <reference name> -x <xmfa file> -s <number of sequences>"
	print '\nWritten by: Simon R Harris, Wellcome Trust Sanger Institute, UK. 2009\n'



#------------------------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------------------------

def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "hr:x:s:", ["help", "ref=", "xmfa="])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	ref=''
	align=''

	
	for opt, arg in opts:
		
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-r", "--ref"):
			ref=arg
		elif opt in ("-x", "--xmfa"):
			align=arg
		elif opt in ("-x", "--xmfa"):
			noseq=int(arg)
	
	contigs=args
	if ref=='':
		print 'no reference file selected'
		Usage()
		sys.exit()
	elif align=='':
		print 'no xmfa file selected'
		Usage()
		sys.exit()
	

	return ref, align

def rev(sequence):
	rev=sequence[::-1]
	
	return rev


if __name__ == "__main__":
	argv=sys.argv[1:]
	alnfiles=''
	ref, align=getOptions(argv)
	refdata=open(ref,'rU').read().split('\n')[1:]
	maxlen=len(''.join(refdata))
	
	contigs=[]

	#os.system('progressiveMauve --collinear --output='+ref+'.mauvealn '+ref+' '+alnfiles)
	alines=open(align, 'rU').read().split('=')
	sequfragments={}
	print "Reading Mauve output..."
	for aline in alines[:-1]:
		posn=-1
		
		lines=aline.split('>')
		
		for line in lines[0]:
			if line[0]=="#" and len(line.split())>1 and line[1]!='F':
				contigs.append(line.split()[1])
		
		
		if int(lines[1].strip().split(':')[0])!=1:
			continue
		
		dirn='+'
		
		for line in lines[1:]:
		
			if len(line)==0:
				continue
			number=int(line.strip().split(':')[0])-1
			
			if number==0:
				posn=int(line.strip().split(':')[1].split('-')[0])
				sequfragments[posn]=[]
				dirn=line.strip().split()[1]
				for i in range(len(contigs)+1):
						sequfragments[posn].append('')
			
			if dirn=='+':
				sequfragments[posn][number]=sequfragments[posn][number]+''.join(line.strip().split('\n')[1:])
			elif dirn=='-':
				sequfragments[posn][number]=sequfragments[posn][number]+rev(''.join(line.strip().split('\n')[1:]))
			
	
	
	keys=sequfragments.keys()
	keys.sort()
	
	sequences={}
	seqnames=[]
	for name in contigs:
		sequences[name]=''
		seqnames.append(name)
	
	print "Creating fasta alignment..."
	
	for key in keys:
		fraglen=0
		for x, fragment in enumerate(sequfragments[key]):
			if x==0:
				sequences[ref.split('.')[0]]=sequences[ref.split('.')[0]]+fragment
				fraglen=len(fragment)
			elif len(fragment)==fraglen:
				sequences[contigs[x-1]]=sequences[contigs[x-1]]+fragment
			else:
				sequences[contigs[x-1]]=sequences[contigs[x-1]]+'-'*fraglen
	
	
	
	
	
	
	alnout=open(ref+'_mauve.aln','w')
	
	for key in seqnames:
		print >> alnout, '>'+key
		print >> alnout, sequences[key]
		
	
	alnout.close()

	
	print "Done\n"		
				
