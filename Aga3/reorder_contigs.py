#!/usr/bin/env python
# !/usr/bin/env python
import string
import os, sys
from Bio import Seq
from Bio import AlignIO
from Bio.Align import Generic
from Bio.Alphabet import IUPAC, Gapped
from optparse import OptionParser, OptionGroup
import glob
import dendropy
import shlex, subprocess

sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_general import *
from Si_SeqIO import *



####################
# Set some globals #
####################


MUMMER_DIR=""


#############################################
# Function to reverse complement a sequence #
#############################################

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


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	
	parser.add_option("-c", "--contigs", action="store", dest="contigs", help="Contig fasta file name", default="", metavar="FILE")
	parser.add_option("-o", "--output", action="store", dest="output", help="Output file name", default="", metavar="STRING")
	parser.add_option("-t", "--tree", action="store", dest="tree", help="Tree file to allow ordering by clade", default="")
	
	
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


	if options.contigs=='':
		DoError('No reference fasta file selected (-r)')
	elif not os.path.isfile(options.contigs):
		DoError('Cannot find file '+options.contigs)
	
	
	if options.output=="":
		DoError('No output prefix selected (-o)')
	
	
	return



def read_tree(treefile):
	#Try opening the tree using various schemas
	opened=False
	for treeschema in ["nexus", "newick"]:#["beast-summary-tree",  "nexus", "newick"]:
		try: 
			t = dendropy.Tree.get_from_path(treefile, schema=treeschema, as_rooted=True, preserve_underscores=True, case_insensitive_taxon_labels=False, set_node_attributes=True, extract_comment_metadata=True)
			opened=True
			t.schema=treeschema
			break
		except dendropy.utility.error.DataParseError:
			continue
		except ValueError as e:
			print "Encountered ValueError while trying to read tree file as", treeschema+":", e
			continue
			
	if not opened:
		print "Failed to open tree file"
		sys.exit()
	
	#t.deroot()
	for n in t.leaf_iter():
		print n.taxon.label
	
	print t.as_ascii_plot()
	
	
	return t




################
# Main program #
################		

if __name__ == "__main__":

	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)
	
#	read_tree(options.tree)
#	
#	sys.exit()
	
	
	try:
		refseqs=read_seq_file(options.contigs)
	except StandardError:
		DoError("Cannot open reference file")
	
	reflist=[]
	reflens=[]
	refstarts={}
	refends={}
	lens={}
	seqs={}
	totlen=0
	for seq in refseqs:
		if seq.id in reflist:
			DoError("You have two identically named sequences in your reference fasta file")
		reflist.append(seq.id)
		reflens.append(totlen)
		refstarts[seq.id]=totlen
		refends[seq.id]=totlen+len(seq.seq)
		seqs[seq.id]=str(seq.seq).upper()
		totlen+=len(seq.seq)
		lens[seq.id]=len(seq.seq)

		try:
			refstarts[seq.id]=int(seq.description.split("#")[1].split("..")[0].strip())+1
		except ValueError:
			refstarts[seq.id]=float("Inf")
		except IndexError:
				refstarts[seq.id]=float("Inf")
		
		if refstarts[seq.id]==float("Inf") or refstarts[seq.id]=="":
			try:
				refstarts[seq.id]=int(seq.description.split("#")[2].split("..")[0].strip())+1
			except ValueError:
				refstarts[seq.id]=float("Inf")
			except IndexError:
					refstarts[seq.id]=float("Inf")
		
		if refstarts[seq.id]==float("Inf") or refstarts[seq.id]=="":
			try:
				refstarts[seq.id]=int(seq.description.split("#")[1].split("..")[1].strip())+1
			except ValueError:
				refstarts[seq.id]=float("Inf")
			except IndexError:
				refstarts[seq.id]=float("Inf")
		
		if refstarts[seq.id]==float("Inf") or refstarts[seq.id]=="":
			try:
				refstarts[seq.id]=int(seq.description.split("#")[2].split("..")[1].strip())+1
			except ValueError:
				refstarts[seq.id]=float("Inf")
			except IndexError:
				refstarts[seq.id]=float("Inf")
	
	reforders=[]
	for ref in refstarts:
		reforders.append([refstarts[ref],ref])
	
	reforders.sort()
	output=open(options.output, "w")
	for ref in reforders:
		print ref[1]
		print >> output, ">"+ref[1]
		print >> output, seqs[ref[1]]
	
	
	sys.exit()
	