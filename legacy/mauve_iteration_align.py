#!/usr/bin/env python

#Compares mummer snp output files.  NOTE: must be from the same reference sequence



#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------

import string, re, math
import os, sys, getopt


def Usage():
	print '\nget_MUMmer_out.py Usage:\n'
	print '\t/nfs/pathogen/sh16_scripts/get_MUMmer_out.py [options] input files\nOptions:\n\t-a\tAlignment output file\n\t-c\tReference sequence is circular\n\t-d\tDirty mode: Does not remove MUMmer or Mauve files created during the analysis\n\t-e\tEMBL file(s) for reference (comma seperated list)\n\t-i\tInclude indels\n\t-m\tUse previous MUMmer SNP analysis output files\n\t-o\tOutput file name\n\t-t\tArtemis tab output file name\n\t-r\tReference sequence fasta file\n\t-h\tPrint this help'
	print '\nWritten by: Simon R Harris, Wellcome Trust Sanger Institute, UK. 2008\n'



#------------------------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------------------------

def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "hr:c:", ["help", "ref=", "contigs"])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	ref=''
	outfile=''
	infile=''

	
	for opt, arg in opts:
		
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-r", "--ref"):
			ref=arg
		elif opt in ("-c", "--contigs"):
			infile=arg
	

	if infile=='':
		print 'no contig file selected'
		Usage()
		sys.exit()

	return ref, infile

def rev(sequence):
	rev=sequence[::-1]
	
	return rev


if __name__ == "__main__":
	argv=sys.argv[1:]
	ref, infile=getOptions(argv)
	
	newcontigs=open(infile.split('.')[0]+'_new_contigs.fasta', 'w')
	contigs=open(infile,'rU').read().split('>')[1:]
	for contig in contigs:
		curcontigfile=open('temp_contig.fasta', 'w')
		print >> curcontigfile, '>'+contig
		contigname=contig.split('\n')[0].split()[0]
		
		print "Creating multiple alignment with Mauve for "+contigname+"..."
		os.system('progressiveMauve --output='+ref+'.mauvealn '+ref+' temp_contig.fasta')
		#os.system('progressiveMauve --collinear --output='+ref+'.mauvealn '+ref+' '+alnfiles)
		sys.exit()
		os.system('rm temp_contig.fasta *.sml')
		
		
		print "Reading Mauve output..."
		alines=open(ref+'.mauvealn', 'rU').read().split('=')
		sequfragments={}
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
				
				if number!=2:
					continue
				
				posn=int(line.strip().split(':')[1].split('-')[0])
				dirn=line.strip().split()[1]
			
				if dirn=='+':
					sequfragments[posn]=sequfragments[posn][number]+''.join(line.strip().split('\n')[1:])
				elif dirn=='-':
					sequfragments[posn]=sequfragments[posn][number]+rev(''.join(line.strip().split('\n')[1:]))
				
		#sys.exit()
		
		keys=sequfragments.keys()
		keys.sort()
		
		print "Adding new contigs to fasta file..."
		
		count=65
		for key in keys:
			if len(keys)>1:
				print >> newcontigs, ">"+contigname+chr(count)
			else:
				print >> newcontigs, ">"+contigname
			print >> newcontigs, sequfragments[key]
			count=count+1
		
		
	newcontigs.close()		
	print "Done\n"		
				
