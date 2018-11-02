#!/usr/bin/env python


##################
# Import modules #
##################

import pysam, sys, os
sys.path.extend(map(os.path.abspath, ['/nfs/pathogen/sh16_scripts/modules/']))
from Si_general import *
from Si_SeqIO import *

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


def extract_fastas_from_contig(refname):

	reads={}
	
	for read in samfile.fetch(refname):
		if not read.qname in reads:
			reads[read.qname]=[]
		reads[read.qname].append(read)
	
	sout=open(prefix+"tmpsing.seq", "w")
	pout=open(prefix+"tmppair.seq","w")
	singcount=0
	paircount=0
	
	for readname in reads:
		readpair=reads[readname]
		if len(readpair)==1:
			singcount+=1
			if readpair[0].is_reverse:
				print >> sout, ">"+readpair[0].qname
				print >> sout, revcomp(readpair[0].seq)
			else:
				print >> sout, ">"+readpair[0].qname
				print >> sout, readpair[0].seq
		elif len(readpair)==2 and not readpair[0].is_proper_pair:
			if readpair[0].is_reverse:
				print >> sout, ">"+readpair[0].qname
				print >> sout, revcomp(readpair[0].seq)
			else:
				print >> sout, ">"+readpair[0].qname
				print >> sout, readpair[0].seq
			if readpair[1].is_reverse:
				print >> sout, ">"+readpair[1].qname
				print >> sout, revcomp(readpair[1].seq)
			else:
				print >> sout, ">"+readpair[1].qname
				print >> sout, readpair[1].seq	
		elif len(readpair)==2:
			paircount+=1
			if readpair[0].is_reverse:
				print >> pout, ">"+readpair[1].qname
				print >> pout, readpair[1].seq
				print >> pout, ">"+readpair[0].qname
				print >> pout, revcomp(readpair[0].seq)
			else:
				print >> pout, ">"+readpair[0].qname
				print >> pout, readpair[0].seq
				print >> pout, ">"+readpair[1].qname
				print >> pout, revcomp(readpair[1].seq)
	sout.close()
	pout.close()
	return singcount, paircount
				


samfile=pysam.Samfile(sys.argv[1], "rb")

#pindelout=open(sys.argv[2], "r")

refcontig={}
outseq=[]
refseq=""

seqrecords=read_seq_file(sys.argv[2])

for record in seqrecords:
	refcontig[record.name]=str(record.seq)
	outseq=outseq+["N"]*len(refcontig[record.name])
	


variants=[]
weird_regions=[]
count=0
prefix=sys.argv[3]
tmpout=open(prefix+"_tmpoutput","w")
Sioutput=open(prefix+"_Siout.txt","w")
covplot=open(prefix+"_Sicov.plot", "w")
print >> covplot, "#POS %ALIGNED %ID"
contigs=samfile.references
startbase=0
for contig in contigs:
	refseq+=refcontig[contig]	  
	
	singlen, pairlen=extract_fastas_from_contig(contig)

	if singlen==0 and pairlen==0:
		print >> Sioutput, contig+" 0.0 0.0"
		startbase+=len(refcontig[contig])
		continue
	elif singlen==0:
		os.system("velveth "+prefix+"tmpreads.seq_velvet 31 -shortPaired -fasta "+prefix+"tmppair.seq")
	elif pairlen==0:
		os.system("velveth "+prefix+"tmpreads.seq_velvet 31 -short -fasta "+prefix+"tmpsing.seq")
	else:
		os.system("velveth "+prefix+"tmpreads.seq_velvet 31 -short -fasta "+prefix+"tmpsing.seq -shortPaired -fasta "+prefix+"tmppair.seq")
	os.system("velvetg "+prefix+"tmpreads.seq_velvet -exp_cov auto -cov_cutoff auto")

	output=open(prefix+"tmrefregion.fasta", "w")
	print >> output, ">ref"
	print >> output, refcontig[contig]
	output.close()
	#if int(os.popen("grep -c '^>' tmpreads.seq_velvet/contigs.fa").readlines()[0].strip())==1:
		#continue
	bestpercentid=0.0
	bestunmasked=0.0
	for query in open(prefix+"tmpreads.seq_velvet/contigs.fa", "rU").read().split(">")[1:]:
		queryseq=query.split("\n")
		fgbmask=''
		gbmask=''
		output=open(prefix+"tmpref.fasta", "w")
		queryseq=query.split("\n")
		print >> output, ">ref"
		print >> output, refcontig[contig]
		output.close()
		output=open(prefix+"tmpquery.fasta", "w")
		print >> output, ">"+queryseq[0]
		print >> output, ''.join(queryseq[1:])
		output.close()
		
		output=open(prefix+"tmpreads.fasta", "w")
		print >> output, ">ref"
		print >> output, refcontig[contig]
		queryseq=query.split("\n")
		print >> output, ">"+queryseq[0]
		try:
			words=os.popen("bl2seq -i "+prefix+"tmpref.fasta -j "+prefix+"tmpquery.fasta -p blastn -D 1").readlines()[3].split()
		except IndexError:
			continue
		if int(words[8])>int(words[9]):
			print >> output, ">"+queryseq[0]
			print >> output, revcomp(''.join(queryseq[1:]))
			print "Chose reverse alignment"
		else:
			print >> output, ">"+queryseq[0]
			print >> output, ''.join(queryseq[1:])
			print "Chose forward alignment"
		output.close()
		
		os.system("muscle -in "+prefix+"tmpreads.fasta -out "+prefix+"tmpreads.aln")
		
		#os.system("seaview tmpreads.aln")
		
		
		#os.system("seaview tmpreads.aln")
		
		os.system("Gblocks "+prefix+"tmpreads.aln -k=y -t=d -b3=5 -b4=76 -b5=h -p=n -s=n")
		
		gbmask=''.join(open(prefix+"tmpreads.aln-gbMask","rU").read().split(">")[5].split("\n")[1:]).replace(" ","")
		
		
		print gbmask
		
		gbmasked=set([])
		gbgood=set([])
		for x,base in enumerate(gbmask):
			if base==".":
				gbmasked.add(x)
			elif base=="#":
				gbgood.add(x)
		lines=open(prefix+"tmpreads.aln","rU").read().split(">")[1:]
		seqs=[]
		for line in lines:
			words=line.split('\n')
			seqs.append(''.join(words[1:]).upper())
		
		refbase=0
		#endmatch=0
		#startmatch=0
		#for x in xrange(len(seqs[0])):
		#	base=seqs[0][x]
		#	if base!="-":
		#		refbase+=1
		#		
		#	if refbase>endq:
		#		if x in gbgood:
		#			endmatch+=1
		#		else:
		#			break
		#			
		#	elif refbase<startq:
		#		if x in gbgood:
		#			startmatch+=1
		#		else:
		#			startmatch=0
		#print
		#print startmatch, endmatch
		#print
		#refbase=0
		#if endmatch>100 or startmatch>100:
		percentid=0.0
		unmasked=0.0
		for x in xrange(len(seqs[0])):
			base=seqs[0][x]
			if base!="-":
				refbase+=1
				
			if refbase>len(refcontig[contig]):
				break
				
			elif refbase>0:
				if x in gbmasked:
					variants.append([refbase, base, seqs[1][x], "masked", contig])
				elif base!=seqs[1][x]:
					variants.append([refbase, base, seqs[1][x], contig])
				if base!="-":
					if x in gbmasked:
						outseq[startbase+refbase-1]="?"
					else:
						outseq[startbase+refbase-1]=seqs[1][x]
						unmasked+=1
						if base==seqs[1][x]:
							percentid+=1
					
		percentid=(percentid/len(refcontig[contig]))*100
		
		unmasked=(unmasked/len(refcontig[contig]))*100
		if unmasked>bestunmasked:
			bestunmasked=unmasked
			bestpercentid=percentid
		print contig, percentid, unmasked

	print "BEST MATCH:", contig, bestpercentid, bestunmasked	
	print >> Sioutput, contig, bestpercentid, bestunmasked
	if bestunmasked>0:
		for x in xrange(startbase, startbase+len(refcontig[contig])):
			print >> covplot, x, bestunmasked, bestpercentid
		
	#else:
		#weird_regions.append(midbase)
	count+=1
	#if count>1:
	#	break
	startbase+=len(refcontig[contig])
tmpout.close()
Sioutput.close()
covplot.close()
sys.exit()
