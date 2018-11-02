#!/usr/bin/env python

#/usr/bin/env python


##################
# Import modules #
##################

import string, re, os, sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, Gapped
from Bio.SeqRecord import SeqRecord
from optparse import OptionParser, OptionGroup
sys.path.extend(map(os.path.abspath, ['/nfs/pathogen/sh16_scripts/modules/']))
from Si_general import *
from Si_SeqIO import *
from random import *



##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	
	parser.add_option("-r", "--reference", action="store", dest="ref", help="Reference embl file name", default="", metavar="FILE")
	parser.add_option("-m", "--methylation", action="store", dest="methylation", help="Methylation sites gff file", default="", metavar="FILE")
	parser.add_option("-o", "--output", action="store", dest="output", help="Output file prefix", default="", metavar="STRING")
	
	
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


#	if options.processors>cpu_count():
#		DoError('You do not have '+str(options.processors)+' processors')
#
	if options.ref=='':
		DoError('No reference embl file selected (-r)')
	elif not os.path.isfile(options.ref):
		DoError('Cannot find file '+options.ref)
	
	if options.methylation=='':
		DoError('No methylation gff file selected (-m)')
	elif not os.path.isfile(options.methylation):
		DoError('Cannot find file '+options.methylation)
	
	if options.output=="":
		DoError('No output file name selected (-o)')
		
	if options.output+".gff"==options.ref:
		DoError('Your choice of output file name will overwrite your reference file')
	if options.output+".gff"==options.methylation:
		DoError('Your choice of output file name will overwrite your methylation file')
	
	return


################
# Main program #
################		

if __name__ == "__main__":

	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)


	
	refrecord=open_annotation(options.ref, quiet=True)
	refdict={}
	reffound={}
	reforder=[]
	for feature in refrecord.features:
		if feature.type=="CDS":
			name=""
			for nametype in [ "locus_tag", "systematic_id", "gene"]:
				if nametype in feature.qualifiers:
					name=feature.qualifiers[nametype][0]
					break
			if not name in refdict:
				refdict[name]=feature
				reffound[name]=False
				reforder.append(name)
			else:
				print "Repeated name:", name, "Please note that this may cause errors in the output of the script"
				
	
	
	
	print len(refdict), "CDSs found in reference embl file"
	
	
	methylation_locations=[]
	methylation_headings=["start", "end", "type", "score", "strand"]
	for line in open(options.methylation, "rU"):
		words=line.strip().split()
		if len(words)>0 and words[0][0].strip()!="#" and len(words)==9:
			info=words[8].split(";")
			modtype=words[2]
			#if modtype in [".", "modified_base"]:
			#	continue
			if modtype==".":
				continue
			idatabase={}
			for i in info:
				idatabase[i.split("=")[0]]=i.split("=")[1]
				if not i.split("=")[0] in methylation_headings:
					methylation_headings.append(i.split("=")[0])
			if not "motif" in idatabase:
				continue
			idatabase["contig"]=words[0]
			idatabase["source"]=words[1]
			idatabase["type"]=words[2]
			idatabase["start"]=int(words[3])
			idatabase["end"]=int(words[4])
			idatabase["score"]=float(words[5])
			idatabase["strand"]=words[6]
			idatabase["frame"]=words[7]
			
			if idatabase["start"]!=idatabase["end"]:
				print idatabase
				sys.exit()
			
			methylation_locations.append(idatabase)
	
	print len(methylation_locations), "bases marked as methylated"
	
	racsvoutput=open(options.output+"_intragenic.csv", "w")
	ragffoutput=open(options.output+"_intragenic.gff", "w")
	ercsvoutput=open(options.output+"_intergenic.csv", "w")
	ergffoutput=open(options.output+"_intergenic.gff", "w")
	print >> ergffoutput, "##gff-version 3"
	print >> ragffoutput, "##gff-version 3"
	erheader=methylation_headings+["Distance from previous CDS", "Previous CDS Locus tag", "Previous CDS Systematic ID", "Previous CDS Gene", "Previous CDS strand", "Previous CDS start", "Previous CDS end", "Previous CDS product", "Distance from following CDS", "Following CDS Locus tag", "Following CDS Systematic ID", "Following CDS Gene", "Following CDS strand", "Following CDS start", "Following CDS end", "Following CDS product"]
	raheader=methylation_headings+["Distance from CDS start", "Locus tag", "Systematic ID", "Gene", "CDS strand", "CDS start", "CDS end", "CDS product"]
	print >> ercsvoutput, ",".join(erheader)
	print >> racsvoutput, ",".join(raheader)
	x=0
	prev_feature="NULL"
	first_feature="NULL"
	for feature in refrecord.features:
		if prev_feature=="NULL":
			first_feature=feature
		if feature.type=="CDS":
			qualifiers={}
			
			#print feature.location.start, feature.location.end, qualifiers, feature.strand, x, methylation_locations[x]["end"]
			
			while x<len(methylation_locations) and methylation_locations[x]["start"]<feature.location.end:
				
				if methylation_locations[x]["start"]<feature.location.end and methylation_locations[x]["start"]>feature.location.start:
					situation="Intragenic"
					if feature.strand==1:
						distance=methylation_locations[x]["start"]-feature.location.start
						strand="+"
					else:
						distance=feature.location.end-methylation_locations[x]["start"]
						strand="-"
					for qtype in [ "locus_tag", "systematic_id", "gene", "product"]:
						if qtype in feature.qualifiers:
							qualifiers[qtype]=feature.qualifiers[qtype][0]
						else:
							qualifiers[qtype]="NULL"
					feature_string=[distance, qualifiers["locus_tag"], qualifiers["systematic_id"], qualifiers["gene"], strand, feature.location.start, feature.location.end, qualifiers["product"]]
					
				elif methylation_locations[x]["start"]<feature.location.start:
					
					situation="Intergenic"
					if feature.strand==1:
						strand="+"
					
						distance=feature.location.start-methylation_locations[x]["start"]
						
						qualifiers["distance_to_follinwg_CDS"]=distance
						qualifiers["following_CDS_strand"]=strand
						for qtype in [ "locus_tag", "systematic_id", "gene", "product"]:
							if qtype in feature.qualifiers:
								qualifiers["following_CDS_"+qtype]=feature.qualifiers[qtype][0]
							else:
								qualifiers["following_CDS_"+qtype]="NULL"
						fqlist=[distance, qualifiers["following_CDS_locus_tag"], qualifiers["following_CDS_systematic_id"], qualifiers["following_CDS_gene"], strand, feature.location.start, feature.location.end, qualifiers["following_CDS_product"]]
					else:
						fqlist=["NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL"]
					
					if prev_feature.strand==1:
						pqlist=["NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL"]
					else:
						pstrand="-"
						pdistance=methylation_locations[x]["start"]-prev_feature.location.end
							
						qualifiers["previous_CDS_strand"]=pstrand
						qualifiers["distance_to_previous_CDS"]=pdistance

						for qtype in [ "locus_tag", "systematic_id", "gene", "product"]:
							if qtype in prev_feature.qualifiers:
								qualifiers["previous_CDS_"+qtype]=prev_feature.qualifiers[qtype][0]
							else:
								qualifiers["previous_CDS_"+qtype]="NULL"
						pqlist=[pdistance, qualifiers["previous_CDS_locus_tag"], qualifiers["previous_CDS_systematic_id"], qualifiers["previous_CDS_gene"], pstrand, prev_feature.location.start, prev_feature.location.end, qualifiers["previous_CDS_product"]]
					
					
					feature_string=pqlist+fqlist
					
				
				outstring=[]
				for h in methylation_headings:
					if h in methylation_locations[x]:
						outstring.append(methylation_locations[x][h])
					else:
						outstring.append("NULL")
				
				
				outstring+=feature_string
				
				gffinfo=[]
				for m in methylation_locations[x]:
					gffinfo.append(m+"="+str(methylation_locations[x][m]))
				for q in qualifiers:
					gffinfo.append(q+"="+str(qualifiers[q]))
				if methylation_locations[x]["start"]<feature.location.end and methylation_locations[x]["start"]>feature.location.start:
					print >> racsvoutput, ",".join(map(str,outstring))
					print >> ragffoutput, '\t'.join(map(str,[methylation_locations[x]["contig"], methylation_locations[x]["source"], methylation_locations[x]["type"], methylation_locations[x]["start"], methylation_locations[x]["end"], methylation_locations[x]["score"], methylation_locations[x]["strand"], methylation_locations[x]["frame"], ";".join(map(str,gffinfo))]))
				elif  methylation_locations[x]["start"]<feature.location.start:
					print >> ercsvoutput, ",".join(map(str,outstring))
					print >> ergffoutput, '\t'.join(map(str,[methylation_locations[x]["contig"], methylation_locations[x]["source"], methylation_locations[x]["type"], methylation_locations[x]["start"], methylation_locations[x]["end"], methylation_locations[x]["score"], methylation_locations[x]["strand"], methylation_locations[x]["frame"], ";".join(map(str,gffinfo))]))
				
				x+=1
			
				#sys.exit()
			prev_feature=feature
	
	qualifiers={}
		
	while x<len(methylation_locations):
		
		
		if first_feature.strand==1:
			strand="+"
		
			distance=first_feature.location.start+(len(refrecord.seq)-methylation_locations[x]["start"])
			
			qualifiers["distance_to_follinwg_CDS"]=distance
			qualifiers["following_CDS_strand"]=strand
			for qtype in [ "locus_tag", "systematic_id", "gene", "product"]:
				if qtype in first_feature.qualifiers:
					qualifiers["following_CDS_"+qtype]=first_feature.qualifiers[qtype][0]
				else:
					qualifiers["following_CDS_"+qtype]="NULL"
			fqlist=[distance, qualifiers["following_CDS_locus_tag"], qualifiers["following_CDS_systematic_id"], qualifiers["following_CDS_gene"], strand, feature.location.start, feature.location.end, qualifiers["following_CDS_product"]]
		else:
			fqlist=["NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL"]
		
		if prev_feature.strand==1:
			pqlist=["NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL"]
		else:
			pstrand="-"
			pdistance=methylation_locations[x]["start"]-prev_feature.location.end
				
			qualifiers["previous_CDS_strand"]=pstrand
			qualifiers["distance_to_previous_CDS"]=pdistance
			for qtype in [ "locus_tag", "systematic_id", "gene", "product"]:
				if qtype in prev_feature.qualifiers:
					qualifiers["previous_CDS_"+qtype]=prev_feature.qualifiers[qtype][0]
				else:
					qualifiers["previous_CDS_"+qtype]="NULL"
			pqlist=[pdistance, qualifiers["previous_CDS_locus_tag"], qualifiers["previous_CDS_systematic_id"], qualifiers["previous_CDS_gene"], pstrand, prev_feature.location.start, prev_feature.location.end, qualifiers["previous_CDS_product"]]
		
		
		feature_string=pqlist+fqlist	
		
		outstring=[]
		for h in methylation_headings:
			if h in methylation_locations[x]:
				outstring.append(methylation_locations[x][h])
			else:
				outstring.append("NULL")
		
		outstring.append(situation)
		
		outstring+=feature_string
		
		gffinfo=[]
		gffinfo.append("situation="+situation)
		for m in methylation_locations[x]:
			gffinfo.append(m+"="+str(methylation_locations[x][m]))
		for q in qualifiers:
			gffinfo.append(q+"="+str(qualifiers[q]))
		
		print >> ercsvoutput, ",".join(map(str,outstring))
		print >> ergffoutput, '\t'.join(map(str,[methylation_locations[x]["contig"], methylation_locations[x]["source"], methylation_locations[x]["type"], methylation_locations[x]["start"], methylation_locations[x]["end"], methylation_locations[x]["score"], methylation_locations[x]["strand"], methylation_locations[x]["frame"], ";".join(map(str,gffinfo))]))
		
		x+=1
	
	
	racsvoutput.close()
	ragffoutput.close()
	ercsvoutput.close()
	ergffoutput.close()
	print "Done"
