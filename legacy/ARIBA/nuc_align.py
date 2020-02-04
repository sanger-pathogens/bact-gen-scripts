#!/usr/bin/env python3

import os, sys

from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalOmegaCommandline, ClustalwCommandline, MuscleCommandline, MafftCommandline, TCoffeeCommandline, PrankCommandline
from Bio.Data import CodonTable
from argparse import ArgumentParser

##########################################
# Function to Get command line arguments #
##########################################

def get_translation_table_list():
	translations=list(CodonTable.generic_by_id.keys())+list(CodonTable.generic_by_name.keys())
	return translations

def get_commandline_options():

	genetic_code_choices=get_translation_table_list()
	aligner_choices=["clustalo", "clustalw", "clustalw2", "muscle", "mafft", "t_coffee", "prank"]

	parser = ArgumentParser()
	
	parser.add_argument("-i", "--nuc_file", action="store", dest="nuc_file", help="Input nucleotide sequence file in fasta format", default="", metavar="FILE", required=True)
	parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", help="Output prefix", default="", required=True)
	parser.add_argument("-n", "--noncoding", action="store_false", dest="coding", help="Sequences are noncoding (i.e. do not translate) [default= sequences are coding and will be translated]", default=True)
	parser.add_argument("-a", "--aligner", action="store", dest="aligner", help="Aligner to use", default="clustalo", choices=aligner_choices)
	parser.add_argument("-g", "--genetic_code", action="store", dest="genetic_code", help="Genetic code to use", default=11, choices=genetic_code_choices, metavar="GENETIC_CODE")
	parser.add_argument("-G", "--print_codes", action="store_true", dest="print_codes", help="Print possible Genetic codes and exit", default=False)
	
	
	
	return parser.parse_args()


def translate(infile, prefix, table):
	'''Translate nucleotide sequences to amino acids'''
	print("Translating sequences")
	
	outfilename=prefix+"_aa.fasta"
	
	aa_sequences=[]
	for record in SeqIO.parse(nuc_file, "fasta"):
		aa_record=record.seq.translate(table=table, cds=True, to_stop=False)
		aa_sequences.append(SeqRecord(aa_record, id=record.id, description=record.description))
	
	if len(aa_sequences)==0:
		print("Failed to translate any sequences from nucleotide file")
		sys.exit()
	
	SeqIO.write(aa_sequences, outfilename, "fasta")
	
	return outfilename


def make_alignment(infilename, program, prefix, coding_sequence):
	'''Align sequences using clustalw, mafft or muscle'''	
	print("Aligning with "+program)
	if coding_sequence:
		outfilename=prefix+"_aa_"+aligner+".aln"
	else:
		outfilename=prefix+"_nuc_"+aligner+".aln"
	if program in ["clustalw", "clustalw2"]:
		cline = ClustalwCommandline(program, infile=infilename, outfile=outfilename, output="fasta")
	elif program=="t_coffee":
		cline = TCoffeeCommandline(program, infile=infilename, outfile=outfilename, output="fasta")
	elif program=="prank":
		cline = PrankCommandline(program, d=infilename, o=outfilename, f=8)
		outfilename+=".best.fas"
	elif program=="clustalo":
		cline = ClustalOmegaCommandline("clustalo", infile=infilename, outfile=outfilename, outfmt="fasta")
	elif program=="muscle":
		cline = MuscleCommandline("muscle", input=infilename, out=outfilename)
	elif program=="mafft":
		cline = MafftCommandline("mafft", input=infilename, auto=True)
	try:
		stdout, stderr = cline()
		if program=="mafft":
			with open(outfilename, "w") as handle:
				handle.write(stdout)
	except OSError:
		print("Alignment failed")
		print("Here's the command line that didn't work:")
		print(cline)
		sys.exit()
		
	align = AlignIO.read(outfilename, "fasta")

	return align

def back_translate(aa_alignment, nuc_sequence_file, table, prefix):
	'''Back translate to get aligned nucleotides'''
	print("Back translating alignment")
	outfilename=prefix+"_nuc_"+aligner+".aln"
	seq_dict=SeqIO.index(nuc_sequence_file, "fasta")
	aligned_sequences=[]
	for aa_record in aa_alignment:
		seq_record = seq_dict[aa_record.id]
		aligned_seq=[]
		i=0
		codons=Seq("")
		for aa in aa_record.seq:
			if aa!="-":
				assert aa==seq_record.seq[i:i+3].translate(table=table)
				codons+=seq_record.seq[i:i+3]
				i+=3
			else:
				codons+=Seq("---")
			
		aligned_sequences.append(SeqRecord(seq=codons, id=seq_record.id, description=seq_record.description, name=seq_record.name))
		
	SeqIO.write(aligned_sequences, outfilename, "fasta")

def print_codes():
	codes=get_translation_table_list()
	print("Possible inputs for genetic codes option (-g):")
	for code in codes:
		print('\t'+str(code))
	sys.exit()

args = get_commandline_options()

nuc_file=args.nuc_file
genetic_code=args.genetic_code
aligner=args.aligner
coding=args.coding
output_prefix=args.output_prefix
#clustalw_location="clustalw"
#muscle_location="muscle"
#mafft_location="mafft"

if args.print_codes:
	print_codes()

if coding:
	
	file_to_align=translate(nuc_file, output_prefix, genetic_code)

else:
	file_to_align=nuc_file

alignment=make_alignment(file_to_align, aligner, output_prefix, coding)

if coding:
	back_translate(alignment, nuc_file, genetic_code, output_prefix)
	
