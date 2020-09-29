#!/usr/bin/env python
import string, re, gzip
import os, sys

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

if (len(sys.argv)!=6 and len(sys.argv)!=7) or '-h' in sys.argv[1:]:
	print "fastq2fastq_shredder.py <fasta file (can be zipped and can be multifasta)> <output prefix> <read length> <step> <linear/circular (l/c)> <insert length>"
	
	print "E.g. to create 100bp paired reads from a 350bp insert every 3 bases along the circular genome of ref.fasta.gz"
	print "\tfasta2fastq_shredder.py ref.fasta.gz outfiles 100 3 c 350"
	print "E.g. to create 50bp paired reads from a 200bp insert every base along the linear genome of ref2.fasta"
	print "\tfasta2fastq_shredder.py ref.fasta outfiles 50 1 c 200"
	sys.exit()

if sys.argv[1].split('.')[-1]=='gz':
	refcontigs=gzip.open(sys.argv[1], 'r').read().split(">")[1:]
else:
	refcontigs=open(sys.argv[1], 'rU').read().split(">")[1:]


refs=[]
for refcontig in refcontigs:
	refs.append(''.join(refcontig.split("\n")[1:]))


outprefix=sys.argv[2]
readlen=int(sys.argv[3])
step=int(sys.argv[4])
lc=sys.argv[5]
if lc not in ['l','c']:
	print 'linear/circular option must be either "l" or "c"'
	sys.exit()

if len(sys.argv)==7:
	insert=int(sys.argv[6])

if len(sys.argv)==7:
	outputf=open(outprefix+"_1.fastq", "w")
	outputr=open(outprefix+"_2.fastq", "w")
else:
	outputf=open(outprefix+".fastq", "w")

	

for refcount, ref in enumerate(refs):
	count=0
	if lc=="c" and insert+(readlen*2)>len(ref):
		print "Warning, some of your contigs are shorter than the fragment size. This causes problems when assuming contigs are circularised. Skipping contig", str(refcount+1)
		print "You could try using the linear option"
		continue
	for x in range(0-(readlen+insert),len(ref),step):
		if x+readlen>len(ref):
			continue
		
		count=count+1
		
		if len(sys.argv)==7:
			print >> outputf, "@IL0_0000:1:1:"+str(refcount+1)+":"+str(count)+"#0/1"
		else:
			print >> outputf, "@IL0_0000:1:1:"+str(refcount+1)+":"+str(count)
			
		if x<0 and lc=="l":
			print >> outputf, "N"*readlen
			print >> outputf, "+"
			print >> outputf, "&"*readlen
				
		else:
			if (x+readlen)<0:
				print >> outputf, ref[(x+len(ref)):(x+readlen)]
			elif x<0:
				print >> outputf, ref[(x+len(ref)):]+ref[0:x+readlen]
			else:
				print >> outputf, ref[x:x+readlen]
			print >> outputf, "+"
			print >> outputf, "B"*readlen
		
		if len(sys.argv)==7:
			y=x+readlen+insert
			
			print >> outputr, "@IL0_0000:1:1:"+str(refcount+1)+":"+str(count)+"#0/2"
			
			if lc=="l":
				if y+readlen<len(ref):
					print >> outputr, revcomp(ref[y:y+readlen])
					print >> outputr, "+"
					print >> outputr, "B"*readlen
				else:
					print >> outputr, "N"*readlen
					print >> outputr, "+"
					print >> outputr, "&"*readlen
				
				
			else:	
				if y+readlen<len(ref):
					print >> outputr, revcomp(ref[y:(y+readlen)])
				elif y<len(ref):
					print >> outputr, revcomp(ref[y:]+ref[:(y+readlen)-len(ref)])
				else:
					print >> outputr, revcomp(ref[y-len(ref):(y+readlen)-len(ref)])
				print >> outputr, "+"
				print >> outputr, "B"*readlen
	
outputf.close()
#sys.exit()
if len(sys.argv)==7 and sys.argv[6]=="z":
	outputr.close()
	os.system("gzip "+outprefix+"_1.fastq")
	os.system("gzip "+outprefix+"_2.fastq")
elif sys.argv[6]=="z":
	os.system("gzip "+outprefix+".fastq")
	
		
		
