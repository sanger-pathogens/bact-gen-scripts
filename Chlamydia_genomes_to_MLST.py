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
	print '\nChlamydia_genomes_to_MLST.py Usage:'
	print '\nChlamydia_genomes_to_MLST.py -p <tab delimited primer list file> <Genome fasta files/multifasta>'
	print '\nWritten by Simon R. Harris, Wellcome Trust Sanger Institute, UK. 2009\n'

##############################
# Get command line arguments #
##############################


def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "hp:", ["primers="])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	primerfile=""
	fastafiles=[]

	for opt, arg in opts:
	
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-p", "--primers"):
			primerfile=arg
	
	fastafiles=args
	
	if primerfile=='':
		print 'Error: Missing primer file'
		Usage()
		sys.exit()
	elif fastafiles==[]:
		print 'Error: No fasta files given'
		Usage()
		sys.exit()
	
	
	return primerfile, fastafiles


############################################
# Function to reverse complement sequences #
############################################

def revcomp(sequence):
	rev=sequence[::-1]
	revcomp=''
	d={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'g':'c', 'c':'g'}
	for i in rev:
		if d.has_key(i):
			revcomp=revcomp+d[i]
		else:
			revcomp=revcomp+i
	
	return revcomp



if __name__ == "__main__":
	argv=sys.argv[1:]
	primerfile, fastafiles=getOptions(argv)
	
	
	seqs={}
	for f in fastafiles:

		lines=open(f,'rU').read().split(">")[1:]
		for line in lines:
			if seqs.has_key(line.split('\n')[0].split()[0]):
				print "Error: Two of your sequences have the same name:", name
				sys.exit()
			seqs[line.split('\n')[0].split()[0]]=''.join(line.split('\n')[1:]).upper().replace("-","")

		

	lines=open(primerfile, 'rU').readlines()
	
	header=lines[0]
	
	#genes=header.strip().split()[1:]
	
	lines=lines[1:]
	
	for seq in seqs.keys():
		output=open(seq+"_genome_regions.fasta",'w')
		print seq
		for line in lines:
			if len(line.strip().split())<3:
				continue
			words=line.strip().split()
			currprimer=words[0]
			
			print >> output, ">"+seq+"g_"+currprimer
			#print seqs[seq]
			#print words[0], words[1].upper(), revcomp(words[2].upper())
			
			start=string.find(seqs[seq], words[1].upper())
			end=string.find(seqs[seq], revcomp(words[2].upper()))
			
			if start==-1:
				start=string.find(seqs[seq], revcomp(words[1].upper()))
			if end==-1:
				end=string.find(seqs[seq], words[2].upper())
			
			if start==-1 or end==-1:
				print "Error: Cannot find primers for", currprimer, "in sequence", seq
				continue
			if end>start:
				print >> output, seqs[seq][start+len(words[1]):end]
				#print seqs[seq][start:start+20], seqs[seq][end-20:end]
			else:
				print >> output, revcomp(seqs[seq][end+len(words[2]):start])
				#print seqs[seq][end:end+20], seqs[seq][start-20:start]
			
			#print seq, start, end
		
		
	
	
		output.close()
	
	
	
	
	
	


