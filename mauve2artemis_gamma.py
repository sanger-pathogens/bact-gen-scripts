#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math

def Usage():
	print "Usage:"
	print "mauve2artemis_gamma.py <Mauve backbone file> <multifasta file> <output embl file name>"
	print "NOTE: the order of sequences in the multifasta file must be the same as in the Mauve backbone file"


if len (sys.argv)!=4 or '-h' in sys.argv[1:]:
	Usage()
	sys.exit()


bbfile=sys.argv[1]	
fastafile=sys.argv[2]
outfile=sys.argv[3]

lines=open(fastafile, "rU").read().split('>')[1:]

sequences={}
seqlist=[]

for line in lines:
	words=line.strip().split('\n')
	if sequences.has_key(words[0].split()[0]):
		print "Error: You have more than one sequence in your multifasta called "+words[0].split()[0]
		Usage()
		sys.exit()
	sequences[words[0].split()[0]]=''.join(words[1:])
	seqlist.append(words[0])

lines=open(bbfile, "rU").readlines()[1:]

blocksbystrain={}
for x in range(0, len(seqlist)):
	blocksbystrain[x]=[]

for line in lines:
	words=line.strip().split()
	seedstrain=-1
	position=[]
	instrains=[]
	for x in range(0,len(words),2):
		if words[x]!='0':# and (int(words[x+1])-int(words[x]))!=0:
			if seedstrain==-1:
				seedstrain=x/2
				if int(words[x].replace('-',''))<int(words[x+1].replace('-','')):
					position.append(int(words[x].replace('-','')))
					position.append(int(words[x+1].replace('-','')))
				else:
					position.append(int(words[x+1].replace('-','')))
					position.append(int(words[x].replace('-','')))
			position.append(x/2)
	
	if seedstrain!=-1:
		blocksbystrain[seedstrain].append(position)
	
output=open(outfile, 'w')

pan_genome_seq=sequences[seqlist[0]]
curr_base=0
colour_scaling_float=float(255)/(len(seqlist)-1)

for x in range(0, len(seqlist)):
	sortedblocks=blocksbystrain[x]
	sortedblocks.sort()
	
	if x==0:
		print >> output, "FT   source          "+str(sortedblocks[0][0])+".."+str(len(sequences[seqlist[0]]))
		print >> output, "FT                   /strain="+seqlist[x]
		curr_base=len(sequences[seqlist[0]])
	
	else:
		start_base=curr_base+1
	
	for block in sortedblocks:
		if x>0:
			#curr_base=curr_base+1
			if block[1]>len(sequences[seqlist[x]]):
				print "Error: Block out of range. Are you sure you have the sequences in your multifasta in the correct order?"
				print "They must be in the same order as the Mauve xmfa file header"
				print "Sequence "+seqlist[x]+" does not have "+str(block[1])+" bases"
				Usage()
				sys.exit()
			oldlen=len(pan_genome_seq)
			pan_genome_seq=pan_genome_seq+sequences[seqlist[x]][block[0]-1:block[1]]
			#print len(pan_genome_seq), curr_base, curr_base+((block[1]+1)-block[0]), sequences[seqlist[x]][block[0]-1:block[1]]
			if len(pan_genome_seq)!=curr_base+((block[1]+1)-block[0]):
				print len(pan_genome_seq), curr_base, curr_base+((block[1]+1)-block[0]), len(pan_genome_seq)-oldlen, block[1], block[0]
				print block
				sys.exit()
			print >> output, "FT   misc_feature    "+str(curr_base)+".."+str(curr_base+((block[1]+1)-block[0]))
			curr_base=curr_base+((block[1]+1)-block[0])
			if len(block)-2==len(seqlist):
				print >> output, "FT                   /present=Core"
			for strain in block[2:]:
				print >> output, "FT                   /present="+seqlist[strain]
			
		else:
			print >> output, "FT   misc_feature    "+str(block[0])+".."+str(block[1])
			if len(block)-2==len(seqlist):
				print >> output, "FT                   /present=Core"
			for strain in block[2:]:
				print >> output, "FT                   /present="+seqlist[strain]
		
		print >> output, "FT                   /colour="+str(int(255-((len(block)-3)*colour_scaling_float)))+" 0 "+str(int(((len(block)-3)*colour_scaling_float)))
		
	
	if x>0:
		print >> output, "FT   source          "+str(start_base)+".."+str(curr_base)
		print >> output, "FT                   /strain="+seqlist[x]


print >> output, "SQ   Sequence "+str(len(pan_genome_seq))+" BP;"
for x in range(0, len(pan_genome_seq),60):
	print >> output, "    ",
	if len(pan_genome_seq)<x+60:
		seqbit=pan_genome_seq[x:]
	else:
		seqbit=pan_genome_seq[x:x+60]
	spacecount=0
	for y in range(0,len(seqbit),10):
		spacecount=spacecount+1
		print >> output, seqbit[y:y+10],
	print >> output, " "*(9-len(str(x+60)))+" "*(66-(len(seqbit)+spacecount))+str(x+60)

output.close()