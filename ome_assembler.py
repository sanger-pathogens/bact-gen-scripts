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
contigs=samfile.references
startbase=0
for contig in contigs:
	refseq+=refcontig[contig]	  
	
	singlen, pairlen=extract_fastas_from_contig(contig)

	if singlen==0 and pairlen==0:
		startbase+=len(refcontig[contig])
		continue
	elif singlen==0:
		os.system("velveth "+prefix+"tmpreads.seq_velvet 31 -shortPaired -fasta "+prefix+"tmppair.seq")
	elif pairlen==0:
		os.system("velveth "+prefix+"tmpreads.seq_velvet 31 -short -fasta "+prefix+"tmpsing.seq")
	else:
		os.system("velveth "+prefix+"tmpreads.seq_velvet 31 -short -fasta "+prefix+"tmpsing.seq -shortPaired -fasta "+prefix+"tmppair.seq")
	os.system("velvetg "+prefix+"tmpreads.seq_velvet -exp_cov auto -cov_cutoff auto -min_contig_lgth 200")

	
	os.system("cat "+prefix+"tmpreads.seq_velvet/contigs.fa >> "+prefix+"all.fasta")
	#if int(os.popen("grep -c '^>' tmpreads.seq_velvet/contigs.fa").readlines()[0].strip())==1:
		#continue

sys.exit()
