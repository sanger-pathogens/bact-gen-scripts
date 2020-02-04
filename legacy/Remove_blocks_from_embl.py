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
	print '\nRemove_blocks_from_embl.py Usage:\n'
	print '\nRemove_blocks_from_embl.py -s <fasta file> -f <tab file of regions to remove> -e <tab file to be shortened> -o <output file name prefix> -k'
	print '\n-k option tells the program to keep the blocks instead of removing them'
	print '\nCopyright Simon R Harris, Wellcome Trust Sanger Institute, UK. 2009\n'





################################################################
# Recursive algorithm to reconstruct locations with new values #
################################################################

def fixlocationstring(locationstring, start, end):
	
	prefix=""
	suffix=""
#	print locationstring
	if '(' in locationstring:
		prefix=locationstring.split('(')[0]+'('
		parts='('.join(locationstring.split('(')[1:])
		suffix=')'+parts.split(')')[-1]
		parts=')'.join(parts.split(')')[:-1])
		newlocstring=fixlocationstring(parts, start, end)
		
	else:
		joinlist=locationstring.split(',')
		newjoinlist=[]
		for join in joinlist:
			if '..' in join:
				a=int(join.split('..')[0])
				b=int(join.split('..')[1])
				if a>end and b>end:
					newjoinlist.append(str(a-(end+1-start))+'..'+str(b-(end+1-start)))
				elif a<end and a>start and b>end:
					print a, b, start, end, str(start+1)+'..'+str(b-(end+1-start))
					length=(b+1-(end+1-start))-(start+1)
					newstart=start
					while float(length)/3!=length/3:
						newstart=newstart+1
						length=(b+1-(end+1-start))-(newstart+1)
					newjoinlist.append(str(newstart+1)+'..'+str(b-(end+1-start)))
				elif a<start and b<end and b>start:
					length=start-a
					newstart=start
					while float(length)/3!=length/3:
						newstart=newstart-1
						length=newstart-a
					
					#print a, b, start, end, str(a)+'..'+str(start-1)
					newjoinlist.append(str(a)+'..'+str(newstart-1))
				elif a<start and b>end:
					#print a, b, start, str(a)+'..'+str(b-(end+1-start))
					newjoinlist.append(str(a)+'..'+str(b-(end+1-start)))
				elif a<start and b<start:
					newjoinlist.append(str(a)+'..'+str(b))
			else:
				a=int(join)
				if a>end:
					newjoinlist.append(str(a-(end-start)))
				elif a<start:
					newjoinlist.append(str(a))
			
			newlocstring=','.join(newjoinlist)			

	newstring=prefix+newlocstring+suffix

	return newstring




	

##############################
# Get command line arguments #
##############################


def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "hs:f:e:o:k", ["sequence=", "features=", "embl=", "keep"])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	seq=''
	embl1=''
	embl2=''
	outfile=''
	keep='n'

	for opt, arg in opts:
	
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-s", "--sequence"):
			seq=arg
		elif opt in ("-e", "--embl"):
			embl2=arg
		elif opt in ("-f", "--features"):
			embl1=arg
		elif opt in ("-o", "--out"):
			outfile=arg
		elif opt in ("-k", "--keep"):
			keep='y'


	if seq=='' or embl1=='' or embl2=='':
		print 'Error: Missing input file'
		Usage()
		sys.exit()
	if outfile=='':
		print 'Error: Missing output file'
		Usage()
		sys.exit()
	
	
	return seq, embl1, embl2, outfile, keep




if __name__ == "__main__":
        argv=sys.argv[1:]
        seq, embl1, embl2, outfile, keep=getOptions(argv)



featurelines=open(embl1,"rU").readlines()

features=[]
for line in featurelines:
	if len(line.strip().split())>2 and len(line.strip().split()[2].split('..'))==2:
		try:
			start=int(line.strip().split()[2].split('..')[0])
		except ValueError:
			print "sorry, '%s' isn't a valid integer" % line.strip().split()[2].split('..')[0]

		try:
			end=int(line.strip().split()[2].split('..')[1])
		except ValueError:
			print "sorry, '%s' isn't a valid integer" % line.strip().split()[2].split('..')[1]

		for feature in features:
			if (start>feature[0] and start<feature[1]) or (start>feature[0] and start<feature[1]):
				print "Error, there are overlapping features in your file to be removed. Quitting!"
		features.append([start,end])
			
features.sort()
if keep=='n':
	print "Found", len(features), "blocks to remove"


else:
	oldfeatures=features
	features=[]
	for x in range(0,len(oldfeatures)):
		features.append([1,1])
	features.append([1,100000000])

	for x, feature in enumerate(oldfeatures):
		features[x][1]=feature[0]-1
		features[x+1][0]=feature[1]+1
#	print features
	print "Found", len(features), "blocks to remove"
	

print "Removing blocks from embl file:"

lines=open(embl2,"rU").readlines()


for count, feature in enumerate(features[::-1]):
	print str(count)+"\r"
	start=feature[0]
	end=feature[1]
	for x, line in enumerate(lines):
#		print line
		line=line.strip()
		if ".." in line:
			words=line.split()
			for word in words:
				if ".." in word and word[0]!='/':
					newlocstring=fixlocationstring(word,start,end)
#					word, newlocstring, start, end
					if not '..' in newlocstring:
						lines[x]=lines[x].replace(word,"REMOVETHISLINE")
					else:
						lines[x]=lines[x].replace(word,newlocstring)


output=open(outfile,'w')
present='y'
for line in lines:
	if "REMOVETHISLINE" in line:
		present='n'
	elif len(line.split())>2 and ".." in line:
		print >> output, line.strip()
		present='y'
	elif present=='y':
		print >> output, line.strip()
		

output.close()





		
		