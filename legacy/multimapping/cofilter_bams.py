#!/usr/bin/env python3

import os, sys
import pysam
import pyfaidx
from itertools import zip_longest
from optparse import OptionParser, OptionGroup


def DoError(errorstring):
	'''Function to print an error message and exit with status 1'''
	sys.exit("\nError: "+errorstring+"\nFor help use -h or --help\n")
	


def main():
	'''Function to get command line arguments'''
	usage = "usage: %prog [options] <list of fasta assemblies>"
	parser = OptionParser(usage=usage)
	
	parser.add_option("-r", "--references", action="store", dest="reference", help="Reference sequence file containing sequences of all references used for creation of all bams.", default="", metavar="FILE")
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="More verbose output", default=False)
	parser.add_option("-s", "--score", action="store", dest="score", type="choice", choices=["mapping_quality", "cigar_score", "alignment_score"], help="Score to use for comparing read mappings to alternate references. Choose from mapping_quality, cigar_score or alignment score", default="cigar_score")
	parser.add_option("-S", "--sample", action="store", dest="sample", help="Size of read sample to use for choosing best reference. Default is to use all reads.", default=0, type="int")
	
	
	return parser.parse_args()


def check_input_validity(options, args):
	'''Function to sanity check command line arguments'''

	if len(args)<2:
		DoError("At least 2 sam or bam files must be specified")

	if not os.path.isfile(options.reference) and options.score=="alignment_score":
			DoError("Cannot find file "+options.reference+", which is required when alignment_score is chosen as the scoring method")
	
	return


def get_read_names(reads):
	'''Function to return a list of read names for all reads'''
	names=[]
	for r in reads:
		if r.is_read1:
			d="1"
		else:
			d="2"
		names.append(r.qname+"/"+d)
	return names



def calculate_flag(read):
	'''Function to calculate and return the flag value for a read'''
	flag=0
	if read.is_paired:
		flag+=1
	if read.is_proper_pair:
		flag+=2
	if read.is_unmapped:
		flag+=4	
	if read.mate_is_unmapped:
		flag+=8
	if read.is_reverse:
		flag+=16
	if read.mate_is_reverse:
		flag+=32
	if read.is_read1:
		flag+=64
	if read.is_read2:
		flag+=128
	if read.is_secondary:
		flag+=256
	if read.is_qcfail:
		flag+=512
	if read.is_duplicate:
	 	flag+=1024
	if read.is_supplementary:
		flag+=2048
	return flag

def read_samfile(filename):
	'''Reads a samfile and returns an AlifnmentFile object'''
	try:
		samfile = pysam.AlignmentFile( filename )
	except:
		DoError(filename+" is not a bam or sam file")
	return samfile


def unmap_read(read, mate, unmap_read, unmap_mate):

	#a.flag = 99
	if unmap_read:
		read.mapping_quality = 0
		read.cigartuples = None
		read.cigarstring = None
		read.is_proper_pair = False
		read.is_unmapped = True
		read.pos=0
		read.reference_id=-1
		read.template_length=0
		read.is_reverse=False
		
		mate.is_proper_pair = False
		mate.mate_is_unmapped=True
		mate.next_reference_start=0
		mate.next_reference_id=-1
		mate.template_length=0
		mate.mate_is_reverse=False
		

	if unmap_mate:
		mate.mapping_quality = 0
		mate.cigartuples = None
		mate.cigarstring = None
		mate.is_proper_pair = False
		mate.is_unmapped = True
		mate.pos=0
		mate.reference_id=-1
		mate.template_length=0
		mate.is_reverse=False
		
		read.is_proper_pair = False
		read.mate_is_unmapped=True
		read.next_reference_start=0
		read.next_reference_id=-1
		read.template_length=0
		read.mate_is_reverse=False

	read.flag=calculate_flag(read)
	mate.flag=calculate_flag(mate)
	

def calculate_alignment_score(read, reference_sequence):
	'''Function to calculate alignment score for a read against the reference'''
	alignment_score=0
	readseq=read.query_sequence
	refseq=reference_sequence[read.reference_start:read.reference_end]
	for base in read.get_aligned_pairs():
		if base[0]==None or base[1]==None:
			continue
		elif str(readseq[base[0]].upper())==str(refseq[base[1]-read.reference_start]):
			alignment_score+=1
		
	return alignment_score

def calculate_cigar_score(cigar):
	'''Function to calculate alignment score for a read against the reference'''
	alignment_score=0
	
	if cigar==None:
		return 0

	for element in cigar:
		if element[0]==0:
			alignment_score+=element[1]
		
	return alignment_score

def get_score(read, fastq=None):
	if options.score=="mapping_quality":
			score=read.mapping_quality
	elif options.score=="cigar_score":
		score=calculate_cigar_score(read.cigartuples)
	elif options.score=="alignment_score":
		if fastq==None:
			DoError("fastq must be provided if requesting alignment_score as scoring method")
		score=calculate_alignment_score(read, fastq)

	return score


#Get command line arguments

(options, args) = main()
	
#Do some checking of the input files
	
check_input_validity(options, args)

if options.score=="alignment_score":
	reference_fastas=pyfaidx.Fasta(options.reference, sequence_always_upper=True)

	refs={}
	for ref in reference_fastas:
		refs[ref.name]=str(reference_fastas[ref.name])

#print(refs)

bams=[]

for filename in args:
	bams.append(read_samfile(filename))

if len(bams)<2:
		DoError("Cannot cofilter if fewer than 2 correctly formatted sam or bamfiles are provided")

outfiles=[]
for i in range(len(bams)):
	filename=bams[i].filename.decode("utf-8") 
	outfiles.append(pysam.AlignmentFile( filename+".cofiltered.bam", mode="wb", template=bams[i] ))


best=[]
done=0
for i in range(len(bams)):
	best.append([0, i])

for count, sam_records in enumerate(zip_longest(*bams, fillvalue="-")):
	
	for x in sam_records:
		if type(x)!=pysam.calignmentfile.AlignedSegment:
				DoError(str(count)+" Number of reads in each file is not the same")

	if len(set(get_read_names(sam_records)))>1:
		DoError(str(count)+" Read names do not match: "+', '.join(get_read_names(sam_records)))
	
	
	
	if sam_records[0].is_read1:
		mates=[]
		for i, x in enumerate(sam_records):
			mates.append(next(bams[i]))
	else:
		DoError(str(count)+" Read2 first in pair. Reads not ordered correctly.")

	readscores=[]
	read_maximum=0
	for i in range(len(sam_records)):
		if options.score=="alignment_score":
			score=get_score(sam_records[i], fastq=refs[bams[i].references[sam_records[i].reference_id]])
		else:
			score=get_score(sam_records[i])
		readscores.append(score)
		if score>read_maximum:
			#read_maximum=x.mapping_quality
			read_maximum=score
	
	matescores=[]
	mate_maximum=0
	for i in range(len(mates)):

		if options.score=="alignment_score":
			score=get_score(mates[i], fastq=refs[bams[i].references[mates[i].reference_id]])
		else:
			score=get_score(mates[i])
		matescores.append(score)
		if score>mate_maximum:
			mate_maximum=score
	
	for i in range(0,len(sam_records)):
		if readscores[i]==read_maximum and matescores[i]==mate_maximum:
			best[i][0]+=2
		elif readscores[i]==read_maximum:
			unmap_read(sam_records[i], mates[i], False, True)
			best[i][0]+=1
		elif matescores[i]==mate_maximum:
			unmap_read(sam_records[i], mates[i], True, False)
			best[i][0]+=1

	if options.verbose:
		done+=1
		if done>49999:
		 	done=0
		 	print("Analysed "+str((count+1)*2)+" reads")

	if options.sample!=0:
		if count+1>options.sample:
			break


best.sort()
best.reverse()
if options.verbose:
	print("References ranked by number of reads with maximum score:")
	for x in best:
		print("\t"+bams[x[1]].filename.decode("utf-8")+": "+str(x[0]))




printed=[]
for i in range(len(bams)):
	printed.append(0)
	bams[i].reset()

done=0
for count, sam_records in enumerate(zip_longest(*bams, fillvalue="-")):
	
	for x in sam_records:
		if type(x)!=pysam.calignmentfile.AlignedSegment:
				DoError(str(count)+" Number of reads in each file is not the same")

	if len(set(get_read_names(sam_records)))>1:
		DoError(str(count)+" Read names do not match: "+', '.join(get_read_names(sam_records)))
	
	
	
	if sam_records[0].is_read1:
		mates=[]
		for i, x in enumerate(sam_records):
			mates.append(next(bams[i]))
	else:
		DoError(str(count)+" Read2 first in pair. Reads not ordered correctly.")

	readscores=[]
	read_maximum=0
	for i in range(len(best)):
		if options.score=="alignment_score":
			score=get_score(sam_records[best[i][1]], fastq=refs[bams[best[i][1]].references[sam_records[best[i][1]].reference_id]])
		else:
			score=get_score(sam_records[best[i][1]])
		readscores.append(score)
		if score>read_maximum:
			#read_maximum=x.mapping_quality
			read_maximum=score
	
	matescores=[]
	mate_maximum=0
	for i in range(len(best)):

		if options.score=="alignment_score":
			score=get_score(mates[best[i][1]], fastq=refs[bams[best[i][1]].references[mates[best[i][1]].reference_id]])
		else:
			score=get_score(mates[best[i][1]])
		matescores.append(score)
		if score>mate_maximum:
			mate_maximum=score
	
	read_printed=False
	mate_printed=False
	for i in range(0,len(sam_records)):
		if readscores[i]==read_maximum and matescores[i]==mate_maximum:
			printed[i]+=2
			if read_printed and mate_printed:
				unmap_read(sam_records[best[i][1]], mates[best[i][1]], True, True)
				printed[i]-=2
				continue
			elif read_printed:
				unmap_read(sam_records[best[i][1]], mates[best[i][1]], True, False)
				printed[i]-=1
			elif mate_printed:
				unmap_read(sam_records[best[i][1]], mates[best[i][1]], False, True)
				printed[i]-=1
			outfiles[best[i][1]].write(sam_records[best[i][1]])
			outfiles[best[i][1]].write(mates[best[i][1]])
			read_printed=True
			mate_printed=True
			
		elif readscores[i]==read_maximum and not read_printed:
			unmap_read(sam_records[best[i][1]], mates[best[i][1]], False, True)
			outfiles[best[i][1]].write(sam_records[best[i][1]])
			outfiles[best[i][1]].write(mates[best[i][1]])
			printed[i]+=1
			read_printed=True
		elif matescores[i]==mate_maximum and not mate_printed:
			unmap_read(sam_records[best[i][1]], mates[best[i][1]], True, False)
			outfiles[best[i][1]].write(sam_records[best[i][1]])
			outfiles[best[i][1]].write(mates[best[i][1]])
			printed[i]+=1
			mate_printed=True
		# else:
		# 	unmap_read(sam_records[best[i][1]], mates[best[i][1]], True, True)
		# 	outfiles[best[i][1]].write(sam_records[best[i][1]])
		# 	outfiles[best[i][1]].write(mates[best[i][1]])

	if options.verbose:
		done+=1
		if done>49999:
		 	done=0
		 	print("Analysed "+str((count+1)*2)+" reads")
		# if count>50000:
		#  	break
	
if options.verbose:
	print("Number of reads now aligned to each reference:")
	for i, x in enumerate(best):
		print("\t"+bams[x[1]].filename.decode("utf-8")+": "+str(printed[i]))

for bam in bams:
	bam.close()
for outfile in outfiles:
	outfile.close()
