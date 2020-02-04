#!/usr/bin/env python
import string, re, sets
import os, sys, getopt, random, math


if (len (sys.argv)!=8) or '-h' in sys.argv[1:]:
	print "addsnpstoseq.py <input fasta> <number of snps to add> <number of small insertions to add> <number of small deletions to add> <phage to add> <Sapi to add> <output file>"
	sys.exit()

line=open(sys.argv[1], "rU").read().split('>')[1]

name=line.split('\n')[0].split()[0]
seq=''.join(line.split('\n')[1:])

random.seed()
done=sets.Set()
done.add(-1)
states=["A","C","G","T"]
newseq=seq
extendedseq=seq

print "Adding "+sys.argv[2]+" snps"

for f in range(0,int(sys.argv[2])):
	position=-1
	while position in done:
		position=random.randint(0, len(seq)-1)
	done.add(position)
	state=seq[position].upper()
	while state==seq[position].upper():
		state=states[random.randint(0, 3)]
	newseq=newseq[:position]+state+newseq[position+1:]

print "Adding "+sys.argv[3]+" insertions"
inssize=[1,1,1,1,1,2,2,2,2,3,3,3,4,4,5]

for f in range(0,int(sys.argv[3])):
	position=-1
	while position in done:
		position=random.randint(0, len(seq)-1)
	done.add(position)
	inslen=inssize[random.randint(0, 14)]
	ins=''
	gap=''
	for f in range(0,inslen):
		ins=ins+states[random.randint(0, 3)]
		gap=gap+"-"
	newseq=newseq[:position]+ins+newseq[position:]
	extendedseq=extendedseq[:position]+gap+extendedseq[position:]
	
print "Adding "+sys.argv[4]+" deletions"
for f in range(0,int(sys.argv[4])):
	position=-1
	dellen=inssize[random.randint(0, 14)]
	while position in done:
		position=random.randint(0, len(seq)-(dellen+1))
	done.add(position)
	gap=''
	for f in range(0,dellen):
		newseq=newseq[:position]+"-"+newseq[position+1:]
		position=position+1
		
print "Removing phage1"		
#remove phage1

gap="-"*(419728-376404)
newseq=newseq[:376405]+gap+newseq[419729:]

print "Adding "+sys.argv[5]

phage=open(sys.argv[5],"rU").read().replace('\n','')
position=-1
while position in done:
	position=random.randint(0, len(seq)-1)
gap='-'*len(phage)
newseq=newseq[:position]+phage+newseq[position:]
extendedseq=extendedseq[:position]+gap+extendedseq[position:]

print "Adding "+sys.argv[6]

sapi=open(sys.argv[6],"rU").read().replace('\n','')
position=-1
while position in done:
	position=random.randint(0, len(seq)-1)
gap='-'*len(sapi)
newseq=newseq[:position]+sapi+newseq[position:]
extendedseq=extendedseq[:position]+gap+extendedseq[position:]




output=open(sys.argv[7],'w')

print >> output, ">"+name
print >> output, extendedseq
print >> output, ">"+name+"_a"
print >> output, newseq

output.close()



