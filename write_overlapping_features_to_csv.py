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
	print '\nChlamydia_genomes_to_MLST.py -t <tab file containing features of interest> -k <key in query tab file> -f <comma separated list of features to save> <tab files to find overlapping features in>'
	print '\nWritten by Simon R. Harris, Wellcome Trust Sanger Institute, UK. 2009\n'

##############################
# Get command line arguments #
##############################


def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "ht:f:k:", ["primers=", "features=", "key="])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	tabfile=""
	tabfiles=[]
	features=[]
	key=""

	for opt, arg in opts:
	
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-t", "--primers"):
			tabfile=arg
		elif opt in ("-f", "--features"):
			features=arg.split(',')
		elif opt in ("-k", "--key"):
			key=arg
	
	tabfiles=args
	
	if tabfile=='':
		print 'Error: Missing reference tab file'
		Usage()
		sys.exit()
	if tabfile=='':
		print 'Error: Missing key'
		Usage()
		sys.exit()
	elif tabfiles==[]:
		print 'Error: No query tab files given'
		Usage()
		sys.exit()
	elif features==[]:
		print 'Error: No features to save given'
		Usage()
		sys.exit()
	
	
	return tabfile, tabfiles, features,key


if __name__ == "__main__":
	argv=sys.argv[1:]
	tabfile, tabfiles, features,key=getOptions(argv)
	
	
	queryfeatures=[]
	count=-1
	for line in open(tabfile,'rU'):
		if len(line.strip().split())>2 and line.strip().split()[1]!="source" and len(line.strip().split()[2].split(".."))>1:
			count=count+1
			line=line.replace("complement(","").replace("join(","").replace(")","")
			bits=line.split()[2].split(',')
			location=[]
			for bit in bits:
				start=int(bit.split('..')[0])
				end=int(bit.split('..')[1])
				location.append([start,end])
			queryfeatures.append({"number":count+1})
			queryfeatures[count]["loc"]=location
		elif count>=0 and len(line.strip().split())>1 and line.split()[0]=="FT" and line.split()[1][0]=='/' and len(line.split()[1].split("="))>1:
			queryfeatures[count][line.split()[1].split("=")[0].replace("/","")]=' '.join(line.split()[1:]).split("=")[1].replace('"',"")
	
	for f in tabfiles:
		tabfeatures=[]
		count=-1
		for line in open(f,'rU'):
			if len(line.strip().split())>2 and line.strip().split()[1]!="source" and len(line.strip().split()[2].split(".."))>1:
				count=count+1
				line=line.replace("complement(","").replace("join(","").replace(")","")
				bits=line.split()[2].split(',')
				location=[]
				for bit in bits:
					start=int(bit.split('..')[0])
					end=int(bit.split('..')[1])
					location.append([start,end])
				tabfeatures.append({})
				tabfeatures[count]["loc"]=location
			elif count>=0 and len(line.strip().split())>1 and line.split()[0]=="FT" and line.split()[1][0]=='/' and len(line.split()[1].split("="))>1:
				tabfeatures[count][line.split()[1].split("=")[0].replace("/","")]=' '.join(line.split()[1:]).split("=")[1].replace('"',"")


	print "Region",
	for l in features:
		print '\t', l,
	print '\n',
	
	for q in queryfeatures:
		for t in tabfeatures:
			for l in q['loc']:
				for l2 in t['loc']:
					if (l[0]>l2[0] and l[0]<l2[1]) or (l[1]>l2[0] and l[1]<l2[1]) or (l[1]>l2[0] and l[1]<l2[1]) or (l[0]<l2[0] and l[1]>l2[1]):
						print q[key],
						for f in features:
							if t.has_key(f):
								print '\t', t[f],
							else:
								print '\t',
						print '\n',

				
				
			
				

