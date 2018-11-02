#!/usr/bin/env python
import string, re, copy
import os, sys
sys.path.extend(map(os.path.abspath, ['/nfs/pathogen/sh16_scripts/modules/']))
from Si_alignment import *
from Si_general import *
from Si_SeqIO import *
from Si_nexus import *
from optparse import OptionParser





##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-a", "--alignment", action="store", dest="alignment", help="alignment file name", default="", metavar="FILE")
	parser.add_option("-p", "--prefix", action="store", dest="prefix", help="prefix for output files", default="")
	parser.add_option("-e", "--embl", action="store", dest="embl", help="Embl/genbank annotation file for reference strain (for dN/dS etc.)", default="", metavar="FILE")
	parser.add_option("-r", "--reference", action="store", dest="reference", help="Name of reference sequence relating to the embl file", default="", metavar="string")
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


#	if options.embl!="" and options.reference=='':
#		DoError('No reference selected! If you give an embl file you must specify which sequence it is linked to (i.e. the reference)')
	if options.alignment=='':
		DoError('No alignment file selected!')
	elif not os.path.isfile(options.alignment):
		DoError('Cannot find file '+options.alignment+'!')
	elif options.embl!='' and not os.path.isfile(options.embl):
		DoError('Cannot find file '+options.embl+'!')
	if options.prefix=='':
		options.prefix=options.alignment.split("/")[-1].split(".")[0]

		
		
	return










################
# Main program #
################		

if __name__ == "__main__":


	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)


#Read the alignment file
	
	try:
		alignment=read_alignment(options.alignment)
	except StandardError:
		DoError("Cannot open alignment file")



	if options.embl!='':
		
		refnum=-1
		if options.reference!="":
			refnum=-1
			for x, taxon in enumerate(alignment):
				if taxon.id==options.reference:
					refnum=x
					break
			
			if refnum==-1:
				DoError("Reference "+options.reference+" is not in your alignment")
		else:
			refnum=0
			options.reference=alignment[refnum].id
		
		
		
			
		try:
			emblrecord=open_annotation(options.embl, remove_gaps_from_sequence(alignment[refnum].seq))
		except StandardError:
			DoError("Cannot open annotation file "+options.embl+" please check the format")
		
		
		try:
			ref_to_alignment, alignment_to_ref=get_ref_to_alignment_translations(options.reference, alignment)
		
		except StandardError:
			DoError("Reference "+options.reference+" is not in your alignment")
	
		
		#SNP_types_using_reference(SNPlocations, "JJA", SeqRecordObject, alignmentObject)
	
		#add dNdS to embl features
		
	
		emblCDSrecord=SeqRecord(emblrecord.seq, id=emblrecord.id, name=emblrecord.name, description=emblrecord.description)
		
		feature_number=0
		
		
		
		for x, feature in enumerate(emblrecord.features):
			#print emblrecord.features[x]
			#print feature
			feature=change_reference_location_to_alignment_location(feature, ref_to_alignment)
			emblrecord.features[x]=feature

			if feature.type=="CDS":
				#feature.qualifiers["dNdS"]={}
				emblCDSrecord.features.append(feature)
				feature_number+=1
				#print emblCDSrecord.features[-1]
				
		
		realignment=muscle_gene_realign(alignment, emblCDSrecord)
		
		if options.prefix!="":
			filename=options.prefix+"_realigned.aln"
		else:
			filename="Realignment.fasta"
		
		handle=open(filename, "w")
		AlignIO.write([realignment], handle, "fasta")
		handle.close()
		
		