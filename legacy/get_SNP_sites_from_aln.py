#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math

def Usage():
	print 'compare_mummer_out.py Usage:'
	print 'compare_mummer_out.py -r=[reference sequence] [input alignment(s)] > [output file pool] {-h}'
	print 'Copyright Simon R Harris, Wellcome Trust Sanger Institute, UK. 2008'

	

#------------------------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------------------------


def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "hs:a:o:r:", ["help", "snps=", "ref=", "out="])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	snpfile=''
	outfile=''
	alnfile=''
	ref=''

	for opt, arg in opts:
	
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-s", "--snps"):
			snpfile=arg
		elif opt in ("-a", "--align"):
			alnfile=arg
		elif opt in ("-o", "--out"):
			outfile=arg
		elif opt in ("-r", "--ref"):
			ref=arg


	inputdirs=args
	


	if snpfile=='':
		print 'Error: No SNP file selected!'
		Usage()
		sys.exit()
	elif alnfile=='':
		print 'Error: No info alignment file selected!'
		Usage()
		sys.exit()
	elif ref=='':
		print 'Error: No reference selected!'
		Usage()
		sys.exit()
	elif outfile=='':
		outfile=alnfile+'.new'
	

	return snpfile, alnfile, outfile, ref

if __name__ == "__main__":
	argv=sys.argv[1:]
	snpfile, alnfile, outfile, ref=getOptions(argv)
	newaln={}
	alignment={}
	snpaln={}
	alines=open(alnfile, 'rU').readlines()
	seqname=''
	for line in alines:
		line=line.strip()
		if len(line)>0 and line[0]=='>':
			seqname=line.split()[0][1:]
			alignment[seqname]=''
			newaln[seqname]=''
			snpaln[seqname]=''
		else:
			alignment[seqname]=alignment[seqname]+line
	
	for key in alignment.keys():
		print key, len(alignment[key])
		alignment[key]=list(alignment[key])
	
	for x, y in enumerate(alignment[ref]):
		if y=='-':
			for key in alignment.keys():
				alignment[key][x]='?'
	
	for key in alignment.keys():
		alignment[key]=''.join(alignment[key])
	
	for key in alignment.keys():
			newaln[key]=alignment[key].replace('?','')
	
	snplines=open(snpfile,'rU').readlines()
	
	for line in snplines[1:]:
		words=line.split()
		posn=int(words[0])-1
		for key in newaln.keys():
			snpaln[key]=snpaln[key]+newaln[key][posn]
	
	output=open(outfile, 'w')
	for key in snpaln.keys():
		if key!=ref:
			print >> output, key+' '*(10-len(key))+snpaln[key]
	output.close()
        
        
        
        	