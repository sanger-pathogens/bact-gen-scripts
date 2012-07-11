#!/usr/bin/env python
#/usr/bin/env python
import string
import os, sys
import pysam
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from optparse import OptionParser
from random import *

sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
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


if sys.argv[1].split(".")[-1]=="bam":
	samfile = pysam.Samfile( sys.argv[1], "rb" )
elif sys.argv[1].split(".")[-1]=="sam":
	samfile = pysam.Samfile( sys.argv[1], "r" )
else:
	print "Not a bam file"
	sys.exit()
	
refs=samfile.references
lengths=samfile.lengths

refstats={}
chromosome=	sys.argv[2]
for read in samfile:
	if not read.is_unmapped:
		indel=False
		pos=read.pos
		for x in read.cigar:
			if x[0] in [1,2]:
				indel=True
				break
			else:
				pos+=x[1]
	
		if indel:			
			
			if not chromosome in refs:
				print "argument 2 must be a contig in the samfile"
				sys.exit()
			
			outfile=open("tmp.fasta","w")
			#Read the alignment file
				
			try:
				alignment=read_alignment(sys.argv[3])
			except StandardError:
				DoError("Cannot open alignment file")
			
			for record in alignment:
				if record.id==chromosome:
					my_seq_record = SeqRecord(Seq(str(record.seq[pos-54:pos+54])))
					my_seq_record.id=record.id
					my_seq_record.description="Reference"
					SeqIO.write([my_seq_record], outfile, "fasta")
		
		
		
			iter = samfile.fetch( chromosome, pos, pos)
			for x in iter:
				
			#	if x.is_reverse:
			#		print revcomp(x.seq)
			#	else:
				print x.pos
				print >> outfile,  ">"+x.qname
				
				print >> outfile, "-"*(x.pos-(pos-54))+x.seq+"-"*((pos+54)-(x.pos+54))
			
			outfile.close()
			
			os.system("muscle -in tmp.fasta -out test.aln")
			os.system("seaview test.aln")
			
		