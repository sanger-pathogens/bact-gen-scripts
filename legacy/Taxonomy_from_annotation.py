#!/usr/bin/env python
import string, re
import os, sys
import subprocess

genus_list={}
species_list={}
strain_list={}
Nonecount=0
matchcount=0

for i, line in enumerate(open(sys.argv[1])):
	if len(line)==0 or line[0]=="#":
		continue
	words=line.strip().split('\t')
	
	
	sys.stdout.flush()
	
	if len(words)<8 or words[2]!="CDS":
		continue
	
	for info in words[8].split(";"):
		print i, "\r",
		parts=info.split(":")
		if len(parts)>3 and parts[-2]=="RefSeq":
			pfetchout=subprocess.check_output( ["pfetch", "-D", parts[-1]])
			try:
				taxon_string = pfetchout.strip().split("[")[1].split("]")[0]
			except IndexError:
				Nonecount += 1
				continue
			matchcount+=1
			if len(taxon_string)==0:
				continue
			genus=taxon_string.split()[0]
			species=taxon_string.split()[1]
			strain=' '.join(taxon_string.split()[1:])
			
			if not genus in genus_list:
				genus_list[genus]=0
				
			if not ' '.join([genus, species]) in species_list:
				species_list[' '.join([genus, species])]=0
			
			if not ' '.join([genus, species, strain]) in strain_list:
				strain_list[' '.join([genus, species, strain])]=0
			
			genus_list[genus]+=1
			species_list[' '.join([genus, species])]+=1
			strain_list[' '.join([genus, species, strain])]+=1

print "\n\nTotal number of genes annotated =", matchcount
		
print "\nGenera:"	
print "Match\tNumber of matches\tPercent of matches"		
l=[]	
for genus in genus_list:
	l.append([genus_list[genus], genus, (float(genus_list[genus])/matchcount)*100])

l.sort()
l.reverse()
for i in l:
	print '\t'.join(map(str,[i[1], i[0], i[2]]))
		
print "\nSpecies:"
print "Match\tNumber of matches\tPercent of matches"
l=[]	
for species in species_list:
	l.append([species_list[species], species, (float(species_list[species])/matchcount)*100])

l.sort()
l.reverse()
for i in l:
	print '\t'.join(map(str,[i[1], i[0], i[2]]))


print "\nStrains:"
print "Match\tNumber of matches\tPercent of matches"
l=[]	
for strain in strain_list:
	l.append([strain_list[strain], strain, (float(strain_list[strain])/matchcount)*100])

l.sort()
l.reverse()
for i in l:
	print '\t'.join(map(str,[i[1], i[0], i[2]]))
