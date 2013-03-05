#!/usr/bin/env python

import os, sys

lines=open(sys.argv[3], "rU").read().split(">")[1:]
strain_seqs={}
for line in lines:
	words=line.split('\n')
	name=words[0].strip().split()[0]
	seq=''.join(words[1:])
	strain_seqs[name]=seq

lengths={}
for line in open(sys.argv[2]):
	line=line.strip()
	words=line.split()
	lengths[words[0]]=int(words[1])
	

lastsubject=''
lastquery=''
keep=False
first_query=True
first_subject=True
for line in open(sys.argv[1]):
	line=line.strip()
	words=line.split()
	
	if words[0]==words[1]:
		continue
	
	if first_query:
		match=[]
		for x in range(lengths[words[0]]):
			match.append('0')
	
	if words[1]!=lastsubject or words[0]!=lastquery:
		if lastsubject and lastquery:
			unmatched_run=0
			max_run=0
			for base in match:
				if base=='0':
					unmatched_run+=1
				else:
					if unmatched_run>max_run:
						max_run=unmatched_run
					unmatched_run=0
			
			if (lengths[lastsubject]>lengths[lastquery] or (lengths[lastsubject]==lengths[lastquery] and lastsubject>lastquery)) and max_run<200:
				keep=False
				unmatched_run_score=max_run
				failed_match=''.join(map(str,match))
			
			if max_run>longest_run:
				longest_run=max_run
			if (not keep or lengths[lastquery]<200) and lastquery in strain_seqs:	
				del strain_seqs[lastquery]
		first_subject=False
		match=[]
		for x in range(lengths[words[0]]):
			match.append('0')
	
	if words[0]!=lastquery:
		#if keep:
		#	print "Keep", lastquery, unmatched_run_score, longest_run, lengths[lastquery]
		#elif not first_query:
		#	print "Remove", lastquery, unmatched_run_score#, failed_match
		if not first_query and (not keep or lengths[lastquery]<200) and lastquery in strain_seqs:	
			del strain_seqs[lastquery]
		lastsubject=''
		first_subject=True
		first_query=False
		keep=True
		searching=True
		unmatched_run_score=0
		longest_run=0
	
		
	if keep==False:
		continue
	
	
	
	fromloc=int(words[6])
	toloc=int(words[7])
	
	if fromloc>toloc:
		start=toloc
		end=fromloc
	else:
		start=fromloc
		end=toloc
	
	for x in range(start-1, end):
		match[x]=1
	
	
	querylen=lengths[words[0]]
	subjectlen=lengths[words[1]]
	
	lastquery=words[0]
	lastsubject=words[1]
	
if not first_query and (not keep or lengths[lastquery]<200) and lastquery in strain_seqs:	
	del strain_seqs[lastquery]


for seq in strain_seqs:
	print ">"+seq
	print strain_seqs[seq]
