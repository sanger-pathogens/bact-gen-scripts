#!/usr/bin/env python

#Compares mummer snp output files.  NOTE: must be from the same reference sequence



#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------

import string, re, math
import os, sys, getopt


def Usage():
	print '\nmauve2clonalframe.py\n'
	print '\nConverts Mauve output to ClonalFrmae input in cases where not all blocks in the Mauve output are present in all sequences\n'
	print '\nmauve2clonalframe.py Usage:\n'
	print '\t/nfs/pathogen/sh16_scripts/mauve2clonalframe.py [options]\nOptions:\n\t-i\tInput file name [Mauve xmfa file]\n\t-o\tOutput file name\n\t-h\tPrint this help'
	print '\nWritten by: Simon R Harris, Wellcome Trust Sanger Institute, UK. 2009\n'



#------------------------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------------------------

def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "hi:o:", ["help", "out=", "infile="])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)

	outfile=''
	infile=''
	mauve='n'
	
	for opt, arg in opts:
		
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-o", "--out"):
			outfile=arg
		elif opt in ("-i", "--infile"):
			infile=arg
		
	
	
	
	if infile=='':
		print 'no input file selected'
		Usage()
		sys.exit()
	if outfile=='':
		outfile=infile+".out"
	
	return outfile, mauve, infile
	
	


if __name__ == "__main__":
	argv=sys.argv[1:]
	outfile, mauve, infile=getOptions(argv)

	#if mauve=='y':
	#	print "Creating alignment with Mauve..."
	#	os.system('progressiveMauve --output='+ref+'.mummeraln '+ref+' '+alnfiles)
	
	lines=open(infile, 'rU').read().split('=')
	Blockseqs={}
	sequencenames=[]
	

	print "Reading Mauve output..."
	for line in lines[:-1]:
		
		sequences=line.strip().split('>')[1:]
		
		block=int(sequences[0].split('\n')[0].split(':')[1].split('-')[0])

		Blockseqs[block]={}		
		
		for sequence in sequences:
		
			if len(sequence.strip())>0:
				if not sequence.split()[0].split(':')[0] in sequencenames:
					sequencenames.append(sequence.split()[0].split(':')[0])
				Blockseqs[block][sequence.split()[0].split(':')[0]]=['>'+sequence.split('\n')[0],'\n'.join(sequence.strip().split('\n')[1:])]
	
	
	output=open(outfile, 'w')
	
	blocksort=Blockseqs.keys()
	noseqs=len(sequencenames)
	blocksort.sort()

	for block in blocksort:
		blocklen=len(Blockseqs[block][Blockseqs[block].keys()[0]][1].replace('\n',''))
		blocksequs=len(Blockseqs[block])
		if blocksequs!=noseqs:
			continue
		linelen=len(Blockseqs[block][Blockseqs[block].keys()[0]][1].split('\n')[0])
		for name in sequencenames:
			if Blockseqs[block].has_key(name):
				print >> output, Blockseqs[block][name][0]
				print >> output, Blockseqs[block][name][1]
			else:
				print "Error, "+name+" missing from block! This shouldn't happen! Quitting."
				sys.exit()
					
		print >> output, "="
	
	output.close()
	
					

#for x in Blockseqs.keys():
#	print x
#	for y in Blockseqs[x].keys():
#		print y, Blockseqs[x][y][0]
#sys.exit()
#	
#	keys=sequfragments.keys()
#	keys.sort()
#	
#	sequences={ref.split('.')[0]:''}
#	for name in subjectnames:
#		sequences[name]=''
#	
#	print "Creating fasta alignment..."
#	
#	for key in keys:
#		fraglen=0
#		for x, fragment in enumerate(sequfragments[key]):
#			if x==0:
#				sequences[ref.split('.')[0]]=sequences[ref.split('.')[0]]+fragment
#				fraglen=len(fragment)
#			elif len(fragment)==fraglen:
#				sequences[subjectnames[x-1]]=sequences[subjectnames[x-1]]+fragment
#			else:
#				sequences[subjectnames[x-1]]=sequences[subjectnames[x-1]]+'-'*fraglen
#	
#	alnout=open(align,'w')
#	
#	for key in sequences.keys():
#		print >> alnout, '>'+key
#		print >> alnout, sequences[key]
#		
#	
#	alnout.close()
#	
#	
