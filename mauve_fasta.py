#!/usr/bin/env python

#Compares mummer snp output files.  NOTE: must be from the same reference sequence



#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------

import string, re, math
import os, sys, getopt


def Usage():
	print '\nget_MUMmer_out.py Usage:\n'
	print '\t/nfs/team81/sh16/scripts/get_MUMmer_out.py [options] input files\nOptions:\n\t-a\tAlignment output file\n\t-c\tReference sequence is circular\n\t-d\tDirty mode: Does not remove MUMmer or Mauve files created during the analysis\n\t-e\tEMBL file(s) for reference (comma seperated list)\n\t-i\tInclude indels\n\t-m\tUse previous MUMmer SNP analysis output files\n\t-o\tOutput file name\n\t-t\tArtemis tab output file name\n\t-r\tReference sequence fasta file\n\t-h\tPrint this help'
	print '\nWritten by: Simon R Harris, Wellcome Trust Sanger Institute, UK. 2008\n'



#------------------------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------------------------

def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "hr:", ["help", "ref="])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	ref=''
	contigs=[]

	
	for opt, arg in opts:
		
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-r", "--ref"):
			ref=arg
	
	contigs=args

	if contigs==[]:
		print 'no contig files selected'
		Usage()
		sys.exit()
	elif ref=='':
		print 'no reference file selected'
		Usage()
		sys.exit()

	return ref, contigs

def rev(sequence):
	rev=sequence[::-1]
	
	return rev


if __name__ == "__main__":
	argv=sys.argv[1:]
	alnfiles=''
	ref, contigs=getOptions(argv)
	refdata=open(ref,'rU').read().split('\n')[1:]
	maxlen=len(''.join(refdata))
	
	for contig in contigs:
		alnfiles=alnfiles+' '+contig

	print "Creating multiple alignment with Mauve..."
	os.system('progressiveMauve --output='+ref+'.mauvealn --backbone-output='+ref+'.backbone '+ref+' '+alnfiles)
	#os.system('progressiveMauve --collinear --output='+ref+'.mauvealn '+ref+' '+alnfiles)
	alines=open(ref+'.mauvealn', 'rU').read().split('=')
	sequfragments={}
	print "Reading Mauve output..."
	for aline in alines[:-1]:
		posn=-1
		
		lines=aline.split('>')
		
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
	
	sequences={ref.split('.')[0]:''}
	seqnames=[ref.split('.')[0]]
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
	
	
	
	
	
	
	alnout=open(ref+'_mauve_iteration.aln','w')
	
	for key in seqnames:
		print >> alnout, '>'+key
		print >> alnout, sequences[key]
		
	
	alnout.close()

	
	print "Done\n"		
				
