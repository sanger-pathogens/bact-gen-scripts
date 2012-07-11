#!/usr/bin/env python

#Compares mummer snp output files.  NOTE: must be from the same reference sequence



##################
# Import modules #
##################

import string, re
import os, sys, getopt


#########
# Usage #
#########

def Usage():
	print '\nremove_nc3_recomb_from_aln.py Usage:'
	print '\nremove_nc3_recomb_from_aln.py -a <alignment file> -o <output file name> <list of tab files containing features to be removed>'
	print '\nWritten by Simon R. Harris, Wellcome Trust Sanger Institute, UK. 2009\n'

##############################
# Get command line arguments #
##############################


def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "ha:o:", ["alignment=", "output="])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	alnfile=""
	tabfiles=[]
	outfile=""

	for opt, arg in opts:
	
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-a", "--alignment"):
			alnfile=arg
		elif opt in ("-o", "--output"):
			outfile=arg
	
	tabfiles=args
	
	if alnfile=='':
		print 'Error: Missing alignment file'
		Usage()
		sys.exit()
	if outfile=='':
		print 'Error: Missing output file'
		Usage()
		sys.exit()
	elif tabfiles==[]:
		print 'Error: No tab files given'
		Usage()
		sys.exit()
	
	
	return alnfile, tabfiles, outfile


if __name__ == "__main__":
	argv=sys.argv[1:]
	alnfile, tabfiles, outfile=getOptions(argv)
	
	
	print '\nReading input alignment...'
	sys.stdout.flush()
	
	sequences={}

	currseq=''

	if os.path.getsize(alnfile)<2000000000:
		lines=open(alnfile, "rU").read().split('>')[1:]
	else:
		lines=[]
		count=-1
		for linea in open(alnfile, "rU"):
			if linea[0]==">":
				count=count+1
				lines.append(linea.split()[0][1:]+'\n')
			else:	
				lines[count]=lines[count]+linea
		linesa=[]
	
	for line in lines:
		words=line.strip().split('\n')
		sequences[words[0]]=''.join(words[1:])
	
	
	print '\nReading tab files...'
	sys.stdout.flush()
	
	for g in tabfiles:
		f=g.split("/")[-1].replace("_rec.tab","")
		if not sequences.has_key(f):
			print "Cannot find",f,"in alignment"
			continue
			
		print "Found",f,"in alignment"

		for line in open(g,'rU'):
			if len(line.strip().split())>2 and line.strip().split()[1]!="source" and len(line.strip().split()[2].split(".."))>1:
				location=line.split()[2]
				start=int(location.split('..')[0])
				end=int(location.split('..')[1])
				sequences[f]=sequences[f][:start-1]+"-"*(end-start)+sequences[f][end-1:]
				
	print '\nPrinting output alignment...'
	sys.stdout.flush()
	
	output=open(outfile,"w")
	
	for name in sequences.keys():
		print >> output, ">"+name
		print >> output, sequences[name]
	
	output.close()
		



				
				
			
				

