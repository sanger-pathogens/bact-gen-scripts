#!/usr/bin/env python
import string, re, gzip
import os, sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
import pylab
import numpy


fastqfiles=sys.argv[1:]
#gfffile=sys.argv[2]

#gfflines=[]
#for line in open(gfffile,"rU"):
#	words= line.split()
#	gfflines.append( [words[-1].replace("locus_tag=",""), int(words[3])-1, int(words[4]), words[6]])
#
#
#record = SeqIO.read(open(fastafile), "fasta", IUPAC.unambiguous_dna)
#
#for line in gfflines:
#	fout=open("Fasta_files/"+line[0]+".fas","a")
#	if line[3]=="-":
#		nuc_seq= record[line[1]:line[2]].seq.reverse_complement()
#		print >>fout, ">Arch\n", nuc_seq
#	else:
#		print >>fout, ">Arch\n", record[line[1]:line[2]].seq
#	fout.close()

N=100000

gc_values=[[0.0]*len(fastqfiles) for n in xrange(N) ]

colours=["b", "g", "r", "c"]
legend_colours=[]

for x, fastqfile in enumerate(fastqfiles):
	handle = open(fastqfile)
	print fastqfile
	for i, seq_record in enumerate(SeqIO.parse(handle, "fastq")):
		if i>=N:
			break
		gc_values[i][x]=float(GC(seq_record.seq))
	#gc_values.append(gc_value)
#gc_values.sort() 
	handle.close()
	legend_colours.append(pylab.Rectangle((0, 0), 1, 1, fc=colours[x]))

#print len(gc_values),min(gc_values),max(gc_values)

#pylab.plot(gc_values)
pylab.Figure()

#print gc_values
#gc_values_transposed=numpy.transpose(gc_values)
#print gc_values_transposed
n, bins, patches = pylab.hist(gc_values, bins=20) 
#pylab.title("%i fastq sequences\nGC%% %0.1f to %0.1f" % (len(gc_values),min(gc_values),max(gc_values)))
print n, bins, patches


pylab.legend(legend_colours, fastqfiles)
#pylab.legend(legend_colours, ["PT2019 Multiplexed","PT2019_2 whole lane with PCR","PT2019_4 whole lane no PCR"], loc="upper left")
pylab.title("Fastq file %GC") 
pylab.xlabel("%GC") 
pylab.ylabel("Frequency") 
pylab.show()
