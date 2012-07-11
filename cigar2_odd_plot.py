#!/usr/bin/env python
import string, re, gzip
import os, sys


if len(sys.argv)!=7 or '-h' in sys.argv[1:]:
	print "cigar2_odd_plot.py <reference fasta file> <cigar file (can be zipped)> <max insert length> <min insert length> <max allowed insert> <output prefix>"
	sys.exit()

print "Calculating reference length...",
sys.stdout.flush()

if sys.argv[1].split('.')[-1]=='gz':
	reflen=len(''.join(gzip.open(sys.argv[1], 'r').read().split('>')[1].split('\n')[1:]))	
else:
	reflen=len(''.join(open(sys.argv[1], 'rU').read().split('>')[1].split('\n')[1:]))
print "Done."

sys.stdout.flush()

maxinsertlen=int(sys.argv[3])
mininsertlen=int(sys.argv[4])
cutoff=int(sys.argv[5])
insert=[0]*reflen
smallinsert=[0]*reflen
bothneg=[0]*reflen
bothpos=[0]*reflen
forward=[0]*reflen
reverse=[0]*reflen
print "Reading cigar file...",
sys.stdout.flush()
if sys.argv[2].split('.')[-1]=='gz':
	lines=gzip.open(sys.argv[2], 'r').readlines()
else:
	lines=open(sys.argv[2], 'r').readlines()
print "Done."
print "Identifying reads that do not map in pairs properly..."
sys.stdout.flush()

numlines=len(lines)
count=0
total=0

for loc, line in enumerate(lines):
	count=count+1
	if count>1000:
		total=total+count
		count=0
		print str(int((float(total)/numlines)*100))+"% complete\r",
	
	words=line.split()
	if len(words)<5:
		continue
	#find reverse reads with no partner and add positions to correct list
	if words[1][-1] in ["R","2"] and (loc==0 or words[1].replace(".R","").replace("/2","")!=lines[loc-1].split()[1].replace(".F","").replace("/1","") or (words[4]=="-" and ((int(words[6])-int(lines[loc-1].split()[7]))>cutoff or (int(words[6])-int(lines[loc-1].split()[7]))<0)) or (words[4]=="+" and ((int(lines[loc-1].split()[6])-int(words[7]))>cutoff or (int(lines[loc-1].split()[6])-int(words[7]))<0))):
		for x in range(int(words[6])-1,int(words[7])):
			if words[4]=="+":
				forward[x]=forward[x]+1
			elif words[4]=="-":
				reverse[x]=reverse[x]+1
	#find forward reads with no partner and add positions to correct list
	elif words[1][-1] in ["F","1"] and (loc==reflen-1 or words[1].replace(".F","").replace("/1","")!=lines[loc+1].split()[1].replace(".R","").replace("/2","") or (words[4]=="-" and ((int(words[6])-int(lines[loc+1].split()[7]))>cutoff or (int(words[6])-int(lines[loc+1].split()[7]))<0)) or (words[4]=="+" and ((int(lines[loc+1].split()[6])-int(words[7]))>cutoff or (int(lines[loc+1].split()[6])-int(words[7]))<0))):
		for x in range(int(words[6])-1,int(words[7])):
			if words[4]=="+":
				forward[x]=forward[x]+1
			elif words[4]=="-":
				reverse[x]=reverse[x]+1
	#Identify pairs that map
	elif loc>0 and words[1].replace(".R","").replace("/2","")==lines[loc-1].split()[1].replace(".F","").replace("/1",""):
		#if insert is too big or too small add positions to the appropriate list

		if words[4]=="-" and lines[loc-1].split()[4]=="+":
			if (int(words[6])-int(lines[loc-1].split()[7]))>maxinsertlen:
				for x in range(int(lines[loc-1].split()[6])-1,int(words[7])):
					insert[x]=insert[x]+1
			elif (int(words[6])-int(lines[loc-1].split()[7]))<mininsertlen:
#
				for x in range(int(lines[loc-1].split()[6])-1,int(words[7])):
					smallinsert[x]=smallinsert[x]+1
		elif words[4]=="+" and lines[loc-1].split()[4]=="-":
			if (int(lines[loc-1].split()[6])-int(words[7]))>maxinsertlen:
				for x in range(int(words[6])-1,int(lines[loc-1].split()[7])):
					insert[x]=insert[x]+1
			elif (int(lines[loc-1].split()[6])-int(words[7]))<mininsertlen:
				for x in range(int(words[6])-1,int(lines[loc-1].split()[7])):
					smallinsert[x]=smallinsert[x]+1
		#Add reads with pairs in the same orientation to appropriate list
		elif words[4]=="+" and lines[loc-1].split()[4]=="+":
			for x in range(int(words[6])-1,int(words[7])):
				bothpos[x]=bothpos[x]+1
			for x in range(int(lines[loc-1].split()[6])-1,int(lines[loc-1].split()[7])):
				bothpos[x]=bothpos[x]+1
		elif words[4]=="-" and lines[loc-1].split()[4]=="-":
			for x in range(int(words[6])-1,int(words[7])):
				bothneg[x]=bothneg[x]+1
			for x in range(int(lines[loc-1].split()[6])-1,int(lines[loc-1].split()[7])):
				bothneg[x]=bothneg[x]+1
			
				
print "100% complete"	
print "Printing output files...",
sys.stdout.flush()	

output=open(sys.argv[6]+"_odd_coverage.plot",'w')

for x in range(0,len(forward)):
	print >> output, forward[x], reverse[x], insert[x], smallinsert[x]# smallinserts are very common, so I've left this out

output.close()

output=open(sys.argv[6]+"_inv_coverage.plot",'w')

for x in range(0,len(forward)):
	print >> output, bothpos[x], bothneg[x]

output.close()
print "Done."
print "Zipping output files...",
sys.stdout.flush()
os.system("gzip -f "+sys.argv[6]+"_odd_coverage.plot "+sys.argv[6]+"_inv_coverage.plot")
print "Done."

#The following code plots the reads from the forward and reverse fastq files
#forward=[0]*reflen
#reverse=[0]*reflen
#for line in lines:
#	words=line.split()
#	if len(words)<5:
#		continue
#	if words[1].split('.')[-1]=='F':
#		for x in range(int(words[6])-1,int(words[7])):
#			forward[x]=forward[x]+1
#	elif words[1].split('.')[-1]=='R':
#		for x in range(int(words[6])-1,int(words[7])):
#			reverse[x]=reverse[x]+1
#
#output=open("test.plot",'w')
#
#for x in range(0,len(forward)):
#	print >> output, forward[x], reverse[x]
#
#output.close()

		