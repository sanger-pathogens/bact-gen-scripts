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
	print '\nsSGA.py Usage:\n'
	print '\nSGA.py -s <fasta file 1> -e <embl file 1> -S <fasta file 2> -E <embl file 2> -o <output alignment file name>'
	print '\nCopyright Simon R Harris, Wellcome Trust Sanger Institute, UK. 2009\n'

	

##############################
# Get command line arguments #
##############################


def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "hs:e:S:E:", ["seq1=", "seq2=", "embl1=", "embl2="])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	seq1=''
	embl1=''
	seq2=''
	embl2=''
	outfile=''

	for opt, arg in opts:
	
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-s", "--ref"):
			seq1=arg
		elif opt in ("-e", "--out"):
			embl1=arg
		elif opt in ("-S", "--ref"):
			seq2=arg
		elif opt in ("-E", "--out"):
			embl2=arg

	
	if seq1=='' orn embl1='' or seq2='' or embl2='':
		print 'Error: Missing input file'
		Usage()
		sys.exit()
	if outfile=='':
		outfile=seq1.split('/')[-1].split('.')[0]+"_"+seq2.split('/')[-1].split('.')[0]
	if align=='n' and raxml=='y':
		align='y'
	
	
	return seq1, embl1, seq2, embl2, outfile



def getCDSseq(embldata, CDSstring, direction, CDScount, feature, CDSposns, CorP):
	
	if '(' in CDSstring:
		if CDSstring.split('(')[0]=='complement':
			if direction=='f':
				direction='r'
			else:
				direction='f'
		parts='('.join(CDSstring.split('(')[1:])
		parts=')'.join(parts.split(')')[:-1])
		getCDSseq(embldata, parts, direction, CDScount, feature, CDSposns, CorP)
		
	else:
		joinlist=CDSstring.split(',')
		for join in joinlist:
			a=int(join.split('..')[0])
			b=int(join.split('..')[1])
			if direction=='f':
				start=a
				end=b+1
				y=1
			else:
				start=b
				end=a-1
				y=-1
			for x in range(start,end,y):
				
				embldata[donelen+x]=CorP
				if words[1]=='CDS':
					CDSposns[CDScount].append((x, direction))


	return embldata, CDSposns



if __name__ == "__main__":
        argv=sys.argv[1:]
        seq1, embl1, seq2, embl2=getOptions(argv)



	# Reading embl file
	
	if embl!='':
		print "\nReading EMBL file...",
		sys.stdout.flush()
		CDSposns={}
		CDSnames={}
		embldata=[]
		donelen=0
		CDScount=0
		pseudoposns={}
		pseudonames={}
		pseudocount=0
		for i in range(reflennogaps+1):
			embldata.append('I')
	
		lines=open(embl,'rU').readlines()
		
		for x, line in enumerate(lines):
			words=line.strip().split()
			if len(words)>2 and words[1] in ('CDS', 'rRNA', 'tRNA') and '..' in words[2]:
				#print words
				location=words[2].replace('<','').replace('>','')
				extralocationlines=1
				foundpool='n'
				y=x+1
				locustag=''
				pseudo='n'
				while foundpool=='n' and lines[y].split()[1] not in ('CDS', 'rRNA', 'tRNA', 'gene', 'misc_feature', 'sig_peptide', 'repeat_unit', 'repeat_region', 'misc_RNA'):
					if lines[y].split()[1][0]!='/' and y==x+extralocationlines:
						location=location+lines[y].split()[1].strip()
						extralocationlines=extralocationlines+1
					if '/colour=11' in lines[y].split()[1]:
						pseudo='y'
					if '/systematic_id' in lines[y].split()[1] or '/locus_tag' in lines[y].split()[1]:
						locustag=lines[y].split()[1].split('"')[1]
						#print locustag
						foundpool='y'
					y=y+1
				
				direction='f'
				feature=words[1]
				if feature=='CDS' and pseudo=='n':
					CDScount=CDScount+1
					CDSposns[CDScount]=[]
					CDSnames[CDScount]=locustag
				elif feature=='CDS':
					pseudocount=pseudocount+1
					pseudoposns[pseudocount]=[]
					pseudonames[pseudocount]=locustag
					embldata, pseudoposns=getCDSseq(embldata, location, direction, pseudocount, feature, pseudoposns, 'P')
				if pseudo=='n'and feature in ['CDS', 'rRNA', 'tRNA']:
					embldata, CDSposns=getCDSseq(embldata, location, direction, CDScount, feature, CDSposns, feature[0])
				
	
		
		snptoCDSname={}
		for i in CDSposns.keys():
			for j in CDSposns[i]:
				snptoCDSname[j[0]]=CDSnames[i]
		
		
#		print pseudoposns.keys()
#		print pseudonames
		
		snptopseudoname={}
		for i in pseudoposns.keys():
			for j in pseudoposns[i]:
				snptopseudoname[j[0]]=pseudonames[i]
			
		print "Done"
		sys.stdout.flush()
	