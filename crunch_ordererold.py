#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math

def revcomp(sequence):
	rev=sequence[::-1]
	revcomp=''
	d={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'g':'c', 'c':'g'}
	for i in rev:
		if d.has_key(i):
			revcomp=revcomp+d[i]
		else:
			revcomp=revcomp+i
	
	return revcomp

def rev(sequence):
	rev=sequence[::-1]
	
	return rev


if len (sys.argv)!=7 or '-h' in sys.argv[1:]:
	print "MUMmer_tiling_to_tab.py <crunch file> <contig fasta file> <reference fasta file> <id cutoff> <score cutoff> <min length>"
	sys.exit()

crunchfile=sys.argv[1]	
contigfile=sys.argv[2]
reffile=sys.argv[3]
percentid=int(sys.argv[4])
score=int(sys.argv[5])
minlength=int(sys.argv[6])

lines=open(crunchfile, "rU").readlines()

crdb={}
contigbits={}
for line in lines:
	words=line.strip().split()
	if int(words[0])<score or int(words[1])<percentid or ((int(words[6])-int(words[5]))<minlength and (int(words[5])-int(words[6]))<minlength):
		continue
	else:
		crdb[int(words[2])]=words
		if contigbits.has_key(words[7]):
			contigbits[words[7]]=contigbits[words[7]]+1
		else:
			contigbits[words[7]]=1
		
keys=crdb.keys()
keys.sort()

lines=open(contigfile, "rU").read().split('>')[1:]

contigbitcount={}
contigs={}

for line in lines:
	words=line.split('\n')
	contigs[words[0].split()[0]]=''.join(words[1:])
	contigbitcount[words[0].split()[0]]=0
	
	
lines=open(reffile, "rU").read()

words=lines.split('\n')
ref=''.join(words[1:])
ref=ref.upper()


crunchout=open(crunchfile.split('.')[0]+'_ordered.'+'.'.join(crunchfile.split('.')[1:]), 'w')
contigout=open(contigfile.split('.')[0]+'_ordered.fasta', 'w')
tabout=open(contigfile.split('.')[0]+'_ordered.tab', 'w')


curcontig=''
curstart=1
curend=1
totalend=0
curcrunch={}
lastdirn='f'
colour='11'

print >> contigout, '>ordered_contigs_all_bases'
print >>tabout, 'ID   contigs'

for i, key in enumerate(keys):
	
	if int(crdb[key][6])>int(crdb[key][5]):
		dirn='f'
	else:
		dirn='r'
		hitstart=crdb[key][5]
		crdb[key][5]=crdb[key][6]
		crdb[key][6]=hitstart

	#weird code to stop small bits being put into a big contig. Improve!
	if i+1<len(keys) and crdb[keys[i+1]][7]==curcontig and crdb[key][7]!=curcontig:
		continue
	
	offset=0
	#Identify where a contig overlaps the previous one. Needs work on rev/revcomp!!!!and crdb[key][7]!=curcontig
	if i>0  and int(crdb[keys[i-1]][3])>int(crdb[key][2]):
		print "found overlap: ", crdb[keys[i-1]][7],crdb[key][7], int(crdb[keys[i-1]][3])-int(crdb[key][2]),
		#print crdb[keys[i-1]]
		#print crdb[key]
		if lastdirn=='f':
			lastseq= contigs[curcontig][(curend-1)-(int(crdb[keys[i-1]][3])-int(crdb[key][2])):curend]
		else:
			lastseq= revcomp(contigs[curcontig][curstart-1:curstart+(int(crdb[keys[i-1]][3])-int(crdb[key][2]))])
		if dirn=='f':
			newseq= contigs[crdb[key][7]][int(crdb[key][5])-1:int(crdb[key][5])+(int(crdb[keys[i-1]][3])-int(crdb[key][2]))]
		else:
			newseq= revcomp(contigs[crdb[key][7]][int(crdb[key][6])-(int(crdb[keys[i-1]][3])-(int(crdb[key][2])-1)):int(crdb[key][6])])
		if lastseq==newseq:
			
			offset=int(crdb[keys[i-1]][3])-int(crdb[key][2])
			countdiffsa=0
			countdiffsb=0
			refseq=ref[int(crdb[key][2])-1:int(crdb[keys[i-1]][3])]
			for x, base in enumerate(lastseq):
				if base!=refseq[x]:
					countdiffsa=countdiffsa+1
			for x, base in enumerate(newseq):
				if base!=refseq[x]:
					countdiffsb=countdiffsb+1
			print "SAME!", offset, countdiffsa, countdiffsb
		else:
			countdiffs=0
			countdiffsa=0
			countdiffsb=0
			for x, base in enumerate(lastseq):
				if base!=newseq[x]:
					countdiffs=countdiffs+1
			refseq=ref[int(crdb[key][2])-1:int(crdb[keys[i-1]][3])]
			for x, base in enumerate(lastseq):
				if base!=refseq[x]:
					countdiffsa=countdiffsa+1
			for x, base in enumerate(newseq):
				if base!=refseq[x]:
					countdiffsb=countdiffsb+1
			
			print "DIFFERENT!", 100-(float(countdiffsa)/len(refseq)*100), 100-(float(countdiffsa)/len(refseq)*100), 100-(float(countdiffsb)/len(refseq)*100)
#			print lastseq
#			print newseq
#			print refseq
			
	
	
	if crdb[key][7]!=curcontig or dirn!=lastdirn:
		#print crdb[key], dirn, lastdirn
		if curcontig!='':
			if colour=='10':
				colour='11'
			else:
				colour='10'
			if lastdirn=='f':
				print >> contigout, contigs[curcontig][curstart-1:(curend-1)-offset]
				
				print >>tabout, "FT   contig          "+str(totalend+1)+".."+str(totalend+(curend-curstart))
				print >>tabout, 'FT                   /systematic_id="'+curcontig+'"'
				print >>tabout, 'FT                   /colour='+colour
			else:
				#print curcontig, curstart, curend, len(contigs[curcontig]), rev(contigs[curcontig][curstart-1:curend-1])[:20]
				print >> contigout, revcomp(contigs[curcontig][(curstart-1)+offset:curend-1])
				print >>tabout, "FT   contig          complement("+str(totalend+1)+".."+str(totalend+(curend-curstart))+')'
				print >>tabout, 'FT                   /systematic_id="'+curcontig+'"'
				print >>tabout, 'FT                   /colour='+colour
				blocklen=curend-curstart
			crunchkeys=curcrunch.keys()
			crunchkeys.sort()
			for crunchkey in crunchkeys:
				if lastdirn=='f':
					
					curcrunch[crunchkey][5]=str(int(curcrunch[crunchkey][5])+totalend+1-curstart)
					curcrunch[crunchkey][6]=str(int(curcrunch[crunchkey][6])+totalend-curstart)
				else:
					#print curcrunch[crunchkey]
					newend=int(curcrunch[crunchkey][5])-curstart
					newstart=int(curcrunch[crunchkey][6])-curstart
					curcrunch[crunchkey][5]=str((blocklen-newstart)+totalend+1)
					curcrunch[crunchkey][6]=str((blocklen-newend)+totalend)
					
					#print curcrunch[crunchkey]
				#print curcrunch[crunchkey]
				print >> crunchout, ' '.join(curcrunch[crunchkey])
			curcrunch={}
		totalend=totalend+(curend-curstart)-offset
		curcontig=crdb[key][7]
		curstart=int(crdb[key][5])
		curend=curstart
	lastdirn=dirn


	if int(crdb[key][5])<curstart:
		curstart=int(crdb[key][5])
	if int(crdb[key][6])>curend:
		curend=int(crdb[key][6])

	
	curcrunch[int(crdb[key][5])]=crdb[key][:]
	
if curcontig!='':
	if colour=='10':
		colour='11'
	else:
		colour='10'
	if lastdirn=='f':
		print >> contigout, contigs[curcontig][curstart-1:curend-1]
		print >>tabout, "FT   contig          "+str(totalend+1)+".."+str(totalend+(curend-curstart))
		print >>tabout, 'FT                   /systematic_id="'+curcontig+'"'
		print >>tabout, 'FT                   /colour='+colour
	else:
		#print curcontig, curstart, curend, len(contigs[curcontig]), rev(contigs[curcontig][curstart-1:curend-1])[:20]
		print >> contigout, revcomp(contigs[curcontig][curstart-1:curend-1])
		print >>tabout, "FT   contig          complement("+str(totalend+1)+".."+str(totalend+(curend-curstart))+')'
		print >>tabout, 'FT                   /systematic_id="'+curcontig+'"'
		print >>tabout, 'FT                   /colour='+colour
		blocklen=curend-curstart
	crunchkeys=curcrunch.keys()
	crunchkeys.sort()
	for crunchkey in crunchkeys:
		if lastdirn=='f':
			
			curcrunch[crunchkey][5]=str(int(curcrunch[crunchkey][5])+totalend+1-curstart)
			curcrunch[crunchkey][6]=str(int(curcrunch[crunchkey][6])+totalend-curstart)
		else:
			#print curcrunch[crunchkey]
			newend=int(curcrunch[crunchkey][5])-curstart
			newstart=int(curcrunch[crunchkey][6])-curstart
			curcrunch[crunchkey][5]=str((blocklen-newstart)+totalend+1)
			curcrunch[crunchkey][6]=str((blocklen-newend)+totalend)
			
			#print curcrunch[crunchkey]
		#print curcrunch[crunchkey]
		print >> crunchout, ' '.join(curcrunch[crunchkey])


print "done"
print "To see the output in Act, run the following command:"
print "act "+contigfile.split('.')[0]+'_ordered.fasta '+crunchfile.split('.')[0]+'_ordered.'+'.'.join(crunchfile.split('.')[1:])+' '+reffile

crunchout.close()
contigout.close()
	
	
	
	
	
	
	
	