#!/usr/bin/env python3

import os, sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalwCommandline, MuscleCommandline, MafftCommandline

nuc_file="porB.fasta"
genetic_code=11

aa_sequences=[]
for record in SeqIO.parse(nuc_file, "fasta"):
	aa_record=record.seq.translate(table=genetic_code, cds=True, to_stop=False)
	aa_sequences.append(SeqRecord(aa_record, id=record.id))
SeqIO.write(aa_sequences, "aa.fasta", "fasta")
cline = ClustalwCommandline("clustalw", infile="aa.fasta")
print(dir(ClustalwCommandline))
