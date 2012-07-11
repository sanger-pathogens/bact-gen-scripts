#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math


if (len (sys.argv)!=6) or '-h' in sys.argv[1:]:
	print "location_to_mapping.py <text files with lines to extract> <snp file> <alignment of real seq with ref> <ref name> <alignment of mapped data with ref>"
	sys.exit()

locations=open(sys.argv[1], "rU").read().split('\n')[:-1]

lines=open(sys.argv[2], "rU").readlines()



names=lines[0].split()[4:-1]


snps={}
for name in names:
	snps[name]=[]
missing={}
for name in names:
	missing[name]=[]

for line in lines[1:]:
	words=line.strip().split()
	for x, strain in enumerate(words[4:-1]):
		if strain!=words[1] and strain!='-':
			snps[names[x]].append(words[0])
		elif strain=='-':
			missing[names[x]].append(words[0])

names.sort()

falsepositives={}
falsenegatives={}
inlocmissing={}
correct={}
unmappedinTW20={}
unmappedinstrain={}
innewphage={}
innewsapi={}

for name in names:
	falsepositives[name]=0
	falsenegatives[name]=0
	inlocmissing[name]=0
	unmappedinTW20[name]=0
	unmappedinstrain[name]=0
	correct[name]=0
	innewphage[name]=0
	innewsapi[name]=0


#for name in names:
#	for snp in snps[name]:
#		if not snp in locations:
#			falsepositives[name]=falsepositives[name]+1


ref=sys.argv[4]
lines=open(sys.argv[3], "rU").read().split('>')[1:]

seqs={}


for line in lines:
	seqs[line.split('\n')[0].split()[0]]=''.join(line.split('\n')[1:])

realname=""
newseqs={}
for name in seqs.keys():
	newseqs[name]=""
	if name!=ref:
		realname=name


start=0
for x in range(0,len(seqs[ref])):
	if seqs[ref][x]=='-':
		for seq in seqs.keys():
			newseqs[seq]=newseqs[seq]+seqs[seq][start:x]
				
		start=x+1
		
for seq in seqs.keys():
	newseqs[seq]=newseqs[seq]+seqs[seq][start:]


seqs={}
lines=open(sys.argv[5], "rU").read().split('>')[1:]

for line in lines:
	seqs[line.split('\n')[0].split()[0]]=''.join(line.split('\n')[1:])







for name in names:
	for location in locations:
		if seqs[name][int(location)-1]=='-':
			unmappedinstrain[name]=unmappedinstrain[name]+1
		elif not location in snps[name]:
			falsenegatives[name]=falsenegatives[name]+1
		else:
			correct[name]=correct[name]+1


for name in names:
	for snp in snps[name]:
		if newseqs[realname][int(snp)-1]=='-':
			unmappedinTW20[name]=unmappedinTW20[name]+1
			if int(snp)>1891230 and int(snp)<1936914:
				innewphage[name]=innewphage[name]+1
			if int(snp)>376285 and int(snp)<391437:
				innewsapi[name]=innewsapi[name]+1
		elif not snp in locations:
			falsepositives[name]=falsepositives[name]+1
			
		


for name in names:
	print name, correct[name], falsepositives[name], falsenegatives[name], unmappedinstrain[name], unmappedinTW20[name], innewsapi[name], innewphage[name]

