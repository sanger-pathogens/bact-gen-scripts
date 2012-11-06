#!/usr/bin/env python
import os, sys, string, numpy
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import PairedFastaQualIterator
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_SeqIO import *
from optparse import OptionParser
#from Bio.SeqIO import PairedFastaQualIterator

##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-g", "--genome_size", action="store", dest="genome_size", help="Expected size of genome", default=3000000, type="int")
	parser.add_option("-c", "--coverage", action="store", dest="coverage", help="Target coverage", default=30, type="int")
	parser.add_option("-f", "--fasta", action="store", dest="fasta", help="Input fasta file name", default="")
	parser.add_option("-q", "--qual", action="store", dest="qual", help="Input qual file name (must relate to fasta input file)", default="")
	parser.add_option("-o", "--output", action="store", dest="output", help="Output file name prefix", default="Subsample")
	parser.add_option("-t", "--type", type="choice", dest="fileformat", choices=["fasta","fastq","fastaqual"], action="store", help="output file format (fasta, fastaqual (fasta plus qual file), fastq) [default=%default]", default="fastaqual", metavar="FORMAT")

	
	
	return parser.parse_args()



def check_input_options(options, args):
	
	if not os.path.isfile(options.fasta):
		print "Cannot find input file:", options.fasta
	if options.qual!="" and not os.path.isfile(options.qual):
		print "Cannot find input file:", options.qual
	if options.fileformat=="fastq":
		print "Fastq output format not yet implemented"
		sys.exit()

(options, args) = main()
check_input_options(options, args)


if options.qual=="":
	try:
		contigs=SeqIO.parse(open(options.fasta), "fasta")
	except StandardError:
		print "Could not open file", options.fasta
		sys.exit()
else:
	contigs = PairedFastaQualIterator(open(options.fasta, "rU"), open(options.qual, "rU"))
	try:
		contigs = PairedFastaQualIterator(open(options.fasta, "rU"), open(options.qual, "rU"))
	except StandardError:
		print "Could not open files", options.fasta, "and", options.qual
		sys.exit()



lengths=[]

target_length=(options.genome_size*options.coverage)
print "Aiming for a total length of", target_length, "bases"

for contig in contigs:

	
	seq=str(contig.seq)
	name=contig.name
	length=len(seq)

	if options.qual=="":
		lengths.append([length, name, seq.upper(), "N"*len(seq)])
	else:
		qual=' '.join(map(str,contig.letter_annotations["phred_quality"]))
		lengths.append([length, name, seq.upper(), qual])

lengths.sort()
lengths.reverse()
total=0
count=0

if options.fileformat in["fasta", "fastaqual"]:
	output=open(options.output+".fasta", "w")
elif options.fileformat=="fastq":
	output=open(options.output+".fastq", "w")
if options.fileformat=="fastaqual":
	qualout=open(options.output+".fasta.qual", "w")
	

for x in lengths:
	if total>target_length:
		break
	total+=x[0]
	count+=1
	print >> output, ">"+x[1]
	print >> output, x[2]
	if options.fileformat=="fastaqual":
		print >> qualout, ">"+x[1]
		print >> qualout, x[3]
#	elif options.fileformat=="fastq":
#		print fastq qualities here
print count, "contigs retained, equating to", total, "bases ("+str(total/options.genome_size)+"x coverage)"
output.close()
if options.fileformat=="fastaqual":
	qualout.close()	
	