#!/usr/bin/env python
import string, re, gzip
import os, sys, getopt, random, math

if len(sys.argv)!=5:
	print "Fix_alignment.py input_alignment reference_name changes_file output_file"
	
alnfile=sys.argv[1]
sequences={}
if not os.path.isfile(alnfile):
	print "Cannot find file:", alnfile
	sys.exit()

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
print "Found", len(sequences.keys()), "sequences"

refname=sys.argv[2]
if not refname in sequences.keys():
	print "Cannot find", refname, "in alignment"
	sys.exit()

fixfile=sys.argv[3]
if not os.path.isfile(fixfile):
	print "Cannot find file:", fixfile
	sys.exit()

lines=open(fixfile,"rU").readlines()


for line in lines:
	words=line.strip().split(",")
	if len(words)<1:
		continue
	elif words[0]=="DEL" and len(words)>4:
		strains = words[4].split(";")
		start=int(words[1].strip())-1
		end=int(words[2].strip())
		for strain in strains:
			strain=strain.strip()
			
			if not strain in sequences.keys():
				print "Warning, cannot find", strain, "in alignment. Skipping DEL at", start, "->", end
			else:
				#print sequences[strain][start-10:end+10]
				insert="-"*(end-start)
				sequences[strain]=sequences[strain][:start]+insert+sequences[strain][end:]
				#print sequences[strain][start-10:end+10]
	
	
	elif words[0]=="SNP_?" and len(words)>4:
		strains = words[4].split(";")
		location=int(words[1].strip())-1
		for strain in strains:
			strain=strain.strip()
			
			if not strain in sequences.keys():
				print "Warning, cannot find", strain, "in alignment. Skipping SNP at", location
			else:
				#print sequences[strain][location-10:location+10]
				sequences[strain]=sequences[strain][:location]+"N"+sequences[strain][location+1:]
				#print sequences[strain][location-10:location+10]
	
	
	
	
	
	
	elif words[0]=="RR" and len(words)>1:
		location=int(words[1].strip())-1
		for strain in sequences.keys():
			if strain==refname:
				continue
			#print sequences[strain][location-10:location+10]
			sequences[strain]=sequences[strain][:location]+"N"+sequences[strain][location+1:]
			#print sequences[strain][location-10:location+10]

	
	
	
	
	elif words[0]=="SNP" and len(words)>4:
		strains = words[4].split(";")
		location=int(words[1].strip())-1
		for strain in strains:
			strain=strain.strip()
			if strain=="":
				continue
			elif not strain in sequences.keys():
				print "Warning, cannot find", strain, "in alignment. Skipping SNP at", location
				print line
			elif sequences[strain][location]==sequences[refname][location]:
				print "Warning, ", strain, "at location", location, "listed as a SNP, but is the same as", refname
	


output=open(sys.argv[4], "w")
keys=sequences.keys()
keys.sort()
for sequence in keys:
	print >> output, ">"+sequence
	print >> output, sequences[sequence]

output.close()
	
print "Done"
	
	
	
	
	
	
	
	