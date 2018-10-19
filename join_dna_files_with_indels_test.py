#!/usr/bin/env python
import os, sys, string
from random import *
from optparse import OptionParser, OptionGroup
import pysam
import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, Gapped, generic_dna, SingleLetterAlphabet
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_SeqIO import *


#################################
# Simple Error Printing Funtion #
#################################

def DoError(ErrorString):
	print "!!!Error:", ErrorString,"!!!"
	sys.exit()


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options] <list of mfa files>"
	parser = OptionParser(usage=usage)

	parser.add_option("-r", "--reference", action="store", dest="ref", help="Reference fasta (or multifasta). Must be the one used for mapping", default="")
	
	parser.add_option("-o", "--output", action="store", dest="output", help="output file prefix", default="")
#	parser.add_option("-m", "--muscle", action="store_true", dest="muscle", help="Try to align regions with muscle [default=%default]", default=False)
	parser.add_option("-c", "--curate", action="store_true", dest="curate", help="Manually curate added insertions using seaview [default=%default]", default=False)
	parser.add_option("-t", "--textfile", action="store", dest="textfile", help="textfile containing input mfa files", default="")
	
	return parser.parse_args()



################
# Main program #
################		

if __name__ == "__main__":

	#Get command line arguments

	(options, args) = main()

	if options.ref=="":
		DoError("Reference file (-r) required")
	
	if options.output=="":
		DoError("Output file (-o) required")
	
	if options.textfile=="":
		DoError("No textfile of input mfa files selected")
	if not os.path.isfile(options.textfile):
		DoError("Cannot find file "+options.textfile)
	
	#Read the reference file
	
	if not os.path.isfile(options.ref):
		DoError("Cannot find file "+options.ref)
	
	#read the reference fasta file
	
	try:
		fasta=SeqIO.parse(open(options.ref), "fasta")
	except StandardError:
		DoError("Cannot open file "+options.ref)
		
	if options.textfile!="" and not os.path.isfile(options.textfile):
		DoError("Cannot find file "+options.textfile)
		
	
	for line in open(options.textfile, "rU"):
		line=line.strip()
		args.append(line)	
	
	
	reforder=[]
	refseq={}
	for sequence in fasta:
		if not refseq.has_key(sequence.id):
			refseq[sequence.id]={}
			reforder.append(sequence.id)
		refseq[sequence.id]=sequence
		refseq[sequence.id].id='.'.join(options.ref.split("/")[-1].split(".")[:-1])
	
	
	#Get all insertion locations
	Insertions={}
	Insertion_locations={}
	Deletions={}
	Deletion_locations={}
	
	for arg in args:
		if arg.split(".")[-1]=="mfa":
			indelfile=arg.replace(".mfa", "_indels.txt")
		elif arg.split(".")[-1]=="dna":
			indelfile=arg.replace(".dna", "_indels.txt")
		else:
			DoError("Input files must end in .dna or .mfa")
		if not os.path.isfile(indelfile):
			print "Cannot find ", indelfile
		else:
			try:
				lines=open(indelfile, "rU").readlines()
			except StandardError:
				print "Cannot open ", indelfile
			
			for line in lines:
				contig=line.split()[0]
				location=int(line.split()[1])-1
				#note that insertions must be +1, as they are AFTER the base
				indeltype=line.strip().split()[2]
				change=line.strip().split()[3]
				if indeltype=="+":
					if not Insertion_locations.has_key(contig):
						Insertion_locations[contig]=[]
						Insertions[contig]={}
					
					if not Insertions[contig].has_key(location):
						Insertion_locations[contig].append(location)
						Insertions[contig][location]={}
					Insertions[contig][location][arg]=change
				elif indeltype=="-":
					if not Deletion_locations.has_key(contig):
						Deletion_locations[contig]=[]
						Deletions[contig]={}
					
					if not Deletions[contig].has_key(location):
						Deletion_locations[contig].append(location)
						Deletions[contig][location]={}
					Deletions[contig][location][arg]=change
	

	
	
	#Read the mfa files for each isolate
	sequences={}
	
	for arg in args:
		
		try:
			fasta=SeqIO.parse(open(arg), "fasta")
		except StandardError:
			DoError("Cannot open file "+arg)
		
		for sequence in fasta:
			
			seqid=""
			
			if sequence.id in reforder:
				seqid=sequence.id
			else:
				x=0
				while x<len(sequence.id.split("_")):
					
					if "_".join(sequence.id.split("_")[x:]) in reforder:
						seqid="_".join(sequence.id.split("_")[x:])
						break
					#print reforder, sequence.id,  seqid
					x+=1
			if seqid=="":
				DoError(seqid+" and "+sequence.id+" not in reference genome")
			if not sequences.has_key(seqid):
				sequences[seqid]={}
			sequences[seqid][arg]=sequence
			sequences[seqid][arg].id='.'.join(arg.split("/")[-1].split(".")[:-1])
		
	chars = string.ascii_letters + string.digits
	tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))
	
	
	#Deal with the deletions
	
	align_size=20
	
	for contig in Deletion_locations:
		Deletion_locations[contig].sort()
		Deletion_locations[contig].reverse()
	
	for contig in sequences:
		if not contig in Deletion_locations:
			continue
		
		for location in Deletion_locations[contig]:
			#print location
			maxdellen=0
			#lens=[]
			for sequence in sequences[contig]:
				if Deletions[contig][location].has_key(sequence):
					if len(Deletions[contig][location][sequence])>maxdellen:	
						maxdellen=len(Deletions[contig][location][sequence])
						#lens.append(Deletions[contig][location][sequence])
			
			alignlen=maxdellen+align_size
			
			if location>=align_size and location<=len(refseq[contig])-(align_size+1):
				start=location-align_size
				end=location+alignlen
			elif location<align_size:
				start=0
				end=location+alignlen
			else:
				end=len(refseq[contig])
				start=location-align_size
			
#			print end-start, maxdellen, align_size, alignlen, lens 
			
			tempseqs=[refseq[contig][start:end]]
			
			for sequence in sequences[contig]:
				
				tempseqs.append(sequences[contig][sequence][start:end])
				#print sequences[contig][sequence][start:end]
				tempseqs[-1].id=sequence
				if Deletions[contig][location].has_key(sequence):
					
					deletionlen=len(Deletions[contig][location][sequence])
					#sequences[contig][sequence].seq=sequences[contig][sequence].seq[:location]+"-"*deletionlen+sequences[contig][sequence].seq[location+deletionlen:]
					if start==0:
						distfromstart=align_size-location
						
					else:
						distfromstart=align_size
					
					startbit=tempseqs[-1].seq[:distfromstart+1].upper()
					
					
					if distfromstart+deletionlen<len(tempseqs[-1].seq):
						endbit=tempseqs[-1].seq[distfromstart+1+deletionlen:].upper()
					else:
						endbit=""
					
#					print tempseqs[-1].seq
#					print location>=align_size , location<=len(refseq[contig])-(align_size+1)
#					print location, len(refseq[contig]), Deletions[contig][location][sequence]
#					print distfromstart, start, end,  startbit, endbit, distfromstart+deletionlen, len(tempseqs[-1].seq)
#					print tempseqs[-1].seq
##			
#					print startbit+"-"*deletionlen+endbit
#					sys.exit()
					tempseqs[-1].seq=startbit+"-"*deletionlen+endbit
			

			SeqIO.write(tempseqs, open(tmpname+".aln","w"), "fasta")
			if options.curate:
				os.system("seaview "+tmpname+".aln")
			
			muscleout=open(tmpname+".aln", "rU").read().split(">")[1:]
			for line in muscleout:
				words=line.split("\n")
				name=words[0].split()[0]
				if not sequences[contig].has_key(name):
					refseq[contig].seq=refseq[contig].seq[:start]+''.join(words[1:])+refseq[contig].seq[end:]
					continue
				else:
				
#					seqline=''.join(words[1:])
#					
#					x=0
#					for y, base in enumerate(seqline):
#						if base=="-":
#							continue
#						elif base!=tempseqscopy[name][x]:
#							seqline=seqline[:y]+tempseqscopy[name][x]+seqline[y+1:]
#						x+=1
				
					sequences[contig][name].seq=sequences[contig][name].seq[:start]+''.join(words[1:])+sequences[contig][name].seq[end:]
			
			refseqnewlen=len(refseq[contig].seq)
			
			for name in sequences[contig]:
				if len(sequences[contig][name].seq)!=refseqnewlen:
					print "D", start, end, location, refseqnewlen, len(sequences[contig][name].seq)
					sys.exit()
					os.system("seaview "+tmpname+".aln")
			
	
	
	#sort the inserts
	
	
	for contig in Insertion_locations:
		Insertion_locations[contig].sort()
		Insertion_locations[contig].reverse()

	#Deal with the inserts (most complicated bit)
	
	for contig in sequences:
		if not contig in Insertion_locations:
			continue
		
		for location in Insertion_locations[contig]:
			
			if location>=align_size and location<=len(refseq[contig])-(align_size+1):
				start=location-align_size
				end=location+align_size
			elif location<align_size:
				start=0
				end=align_size*2
			else:
				end=len(refseq[contig])-1
				start=len(refseq[contig])-(align_size*2)+1
			tempseqs=[refseq[contig][start:end]]
#			tempseqs[-1].id=tmpname
			
			for sequence in sequences[contig].keys():
				tempseqs.append(sequences[contig][sequence][start:end])
				tempseqs[-1].id=sequence
				if Insertions[contig][location].has_key(sequence):
					tempseqs[-1].seq=tempseqs[-1].seq[:align_size+1]+Insertions[contig][location][sequence]+tempseqs[-1].seq[align_size+1:]
					
			
#			if options.muscle:
			SeqIO.write(tempseqs, open(tmpname+".fasta","w"), "fasta")
			os.system("muscle -in "+tmpname+".fasta -out "+tmpname+".aln  >  /dev/null 2>&1 ")
#			else:
#				SeqIO.write(tempseqs, open(tmpname+".aln","w"), "fasta")
			
#			if options.curate:
#				os.system("mv "+tmpname+".aln "+tmpname+"b.aln")
#				muscleout=open(tmpname+"b.aln", "rU").read().split(">")[1:]
#				output=open(tmpname+".aln", "w")
#				for line in muscleout:
#					words=line.split("\n")
#					print >> output, ">"+words[0]
#					print >> output, words[1].replace("N-","NN").replace("-N","NN")
#
#				output.close()
#				os.system("seaview "+tmpname+".aln")
			
			
			if options.curate:
				os.system("seaview "+tmpname+".aln")
			muscleout=open(tmpname+".aln", "rU").read().split(">")[1:]
			for line in muscleout:
				words=line.split("\n")
				name=words[0].split()[0]
				if not sequences[contig].has_key(name):
					refseq[contig].seq=refseq[contig].seq[:start]+''.join(words[1:])+refseq[contig].seq[end:]
				else:
					seqline=''.join(words[1:])
					oseqline=seqline
					len1=len(seqline)
					while ("N-" in seqline) or ("-N" in seqline):
						seqline=seqline.replace("N-","NN").replace("-N","NN")
					if len1!=len(seqline):
						print oseqline, seqline
						sys.exit()
					sequences[contig][name].seq=sequences[contig][name].seq[:start]+seqline+sequences[contig][name].seq[end:]
			
			refseqnewlen=len(refseq[contig].seq)
			
			for name in sequences[contig]:
				if len(sequences[contig][name].seq)!=refseqnewlen:
					print "I", start, end, location
					print name, refseqnewlen, len(sequences[contig][name].seq)
					sys.exit()
					os.system("seaview "+tmpname+".aln")
			
	
	final_sequences=[]
	final_sequences.append(refseq[reforder[0]])
	for contig in reforder[1:]:
		final_sequences[-1].seq=final_sequences[-1].seq+refseq[contig].seq
	
	for sequence in sequences[contig].keys():
		final_sequences.append(sequences[reforder[0]][sequence])
		for contig in reforder[1:]:
			final_sequences[-1].seq=final_sequences[-1].seq+sequences[contig][sequence].seq

	
	
	SeqIO.write(final_sequences, open(options.output,"w"), "fasta")
	
	os.system("rm -f "+tmpname+"*")
	
	
	
	
		
		
	
	

				
