#!/usr/bin/env python
import string, re, gzip
import os, sys
from Bio import SeqIO
from Bio import AlignIO
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_SeqIO import *
from Si_general import *
from optparse import OptionParser


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "%prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-a", "--alignment", action="store", dest="alignment", help="Input alignment file", default="")
	parser.add_option("-r", "--reference", action="store", dest="reference", help="Name of reference sequence in alignment")
	parser.add_option("-d", "--outdir", action="store", dest="outfile", help="Output directory name", default="")

	
	return parser.parse_args()

################
# Main program #
################		

if __name__ == "__main__":

	#Get command line arguments
	
	(options, args) = main()
	
	unknowns=["N", "X", "?"]
	
	gaps=["-"]
	
	try:
		alignment=read_alignment(options.alignment)
		#alignment = AlignIO.read(open(options.alignment), "fasta")
	except ValueError:
		DoError("Cannot open alignment file "+options.alignment+". Is it in the correct format?")
		
	ref=options.reference
	refseq=False
	
	print "Finding reference"
	
	for record in alignment:
		if record.id==ref:
			refseq=str(record.seq).upper()
			print "Found", ref
			break
	
	if not refseq:
		DoError("Cannot find reference in alignment")
		sys.exit()
	
	
	
	
	
	filenames=[]
	
	pwd=os.getcwd()
	
	if not os.path.isdir(options.outfile):
		os.system("mkdir "+options.outfile)
	
	vcfdir=pwd+"/"+options.outfile
	
	for record in alignment:
		if record.id==ref:
			continue
		
		outlist=[]
		
		ininsertion=False
		indeletion=False
		inmissing=False
		
		sequence=str(record.seq).upper()
		refx=0
		for x, refbase in enumerate(refseq):
			if refbase not in gaps:
				refx+=1
			base=sequence[x]
			
			if base in unknowns:
				if ininsertion:
					ininsertion=False
					outlist.append("\t".join(map(str,[ref, insertionstart[0], "99", refseq[insertionstart[1]], sequence[insertionstart[1]]+''.join(insertion), ".", "PASS", "INDEL"])))
					#print insertion
				elif indeletion:
					indeletion=False
					outlist.append("\t".join(map(str,[ref, deletionstart[0], "99", refseq[deletionstart[1]]+''.join(deletion), sequence[deletionstart[1]], ".", "PASS", "INDEL"])))
					
				if inmissing:
					missing.append(refbase)
				else:
					missingstart=[refx-1,x-1]
					missing=[refbase]
					inmissing=True
				
			
			elif base==refbase and base not in unknowns:
				if ininsertion:
					ininsertion=False
					outlist.append("\t".join(map(str,[ref, insertionstart[0], "99", refseq[insertionstart[1]], sequence[insertionstart[1]]+''.join(insertion), ".", "PASS", "INDEL"])))
					#print insertion
				elif indeletion:
					indeletion=False
					outlist.append("\t".join(map(str,[ref, deletionstart[0], "99", refseq[deletionstart[1]]+''.join(deletion), sequence[deletionstart[1]], ".", "PASS", "INDEL"])))
					#print deletion
				elif inmissing:
					inmissing=False
					outlist.append("\t".join(map(str,[ref, missingstart[0], "0", refseq[missingstart[1]]+''.join(missing), sequence[missingstart[1]], ".", "u", "INDEL"])))
				
			
			elif base not in gaps and refbase not in gaps:
				if ininsertion:
					ininsertion=False
					outlist.append("\t".join(map(str,[ref, insertionstart[0], "99", refseq[insertionstart[1]], sequence[insertionstart[1]]+''.join(insertion), ".", "PASS", "INDEL"])))
					#print insertion
				elif indeletion:
					indeletion=False
					outlist.append("\t".join(map(str,[ref, deletionstart[0], "99", refseq[deletionstart[1]]+''.join(deletion), sequence[deletionstart[1]], ".", "PASS", "INDEL"])))
					#print deletion
				elif inmissing:
					inmissing=False
					outlist.append("\t".join(map(str,[ref, missingstart[0], "0", refseq[missingstart[1]]+''.join(missing), sequence[missingstart[1]], ".", "u", "INDEL"])))
				outlist.append("\t".join(map(str,[ref, refx, ".", refbase, base, "99", "PASS", "."])))
			
			elif base in gaps and refbase not in gaps:
				if ininsertion:
					ininsertion=False
					outlist.append("\t".join(map(str,[ref, insertionstart[0], ".", refseq[insertionstart[1]], sequence[insertionstart[1]]+''.join(insertion), "99", "PASS", "INDEL"])))
				elif inmissing:
					inmissing=False
					outlist.append("\t".join(map(str,[ref, missingstart[0], "0", refseq[missingstart[1]]+''.join(missing), sequence[missingstart[1]], ".", "u", "INDEL"])))
					#print insertion
				if indeletion:
					deletion.append(refbase)
				else:
					deletionstart=[refx-1,x-1]
					deletion=[refbase]
					indeletion=True
					
			elif base not in gaps and refbase in gaps:
				if indeletion:
					indeletion=False
					outlist.append("\t".join(map(str,[ref, deletionstart[0], ".", refseq[deletionstart[1]]+''.join(deletion), sequence[deletionstart[1]], "99", "PASS", "INDEL"])))
				elif inmissing:
					inmissing=False
					outlist.append("\t".join(map(str,[ref, missingstart[0], "0", refseq[missingstart[1]]+''.join(missing), sequence[missingstart[1]], ".", "u", "INDEL"])))
					#print deletion
				if ininsertion:
					insertion.append(base)
				else:
					insertionstart=[refx,x-1]
					insertion=[base]
					ininsertion=True
					
				
				
		#print vcf files		
		if len(outlist)>0:
			print record.id, "contains", len(outlist), "variations"
			filenames.append(vcfdir+"/"+record.id+".vcf")
			output=open(vcfdir+"/"+record.id+".vcf","w")
			print >> output, "##fileformat=VCFv4.0"
			print >> output, '##FILTER=<ID=u,Description="Missing or unmapped">'
			print >> output, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
			for line in outlist:
				print >>output, line
			output.close()
		else:
			print record.id, "shows no variation"
		
		#print x, refx, len(outlist)
		
		
	output=open(vcfdir+"/"+"filenames.list", "w")
	#sys.exit()
	for filename in filenames[::-1]:
		print >> output, filename+".gz"
		os.system("bgzip "+filename)
		os.system("tabix -p vcf "+filename+".gz")
	output.close()
		