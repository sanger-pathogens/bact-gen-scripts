#!/usr/bin/env python
import string, re, gzip
import os, sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
import numpy


fastqfiles=sys.argv[1:]


N=100000
gc_values=[]

legend_colours=[]

for x, fastqfile in enumerate(fastqfiles):
	handle = open(fastqfile)
	for i, seq_record in enumerate(SeqIO.parse(handle, "fastq")):
		if i>100000:
			break
		gc_values.append(float(GC(seq_record.seq)))
	print fastqfile, numpy.mean(gc_values), numpy.std(gc_values)

	handle.close()
