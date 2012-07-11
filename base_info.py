#!/usr/bin/env python
import string, re
import os, sys
import math
from Bio import SeqIO
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.SeqUtils import GC
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import GenBank
from Bio.Data import CodonTable
from optparse import OptionParser
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_general import *
from Si_SeqIO import *


##########################################
# Function to Get command line arguments #
##########################################

def get_user_options():
	usage = "usage: %prog [options] args"
	parser = OptionParser(usage=usage)
	
	parser.add_option("-o", "--output", action="store", dest="outputfile", help="output file name", type="string", default="", metavar="FILE")
	parser.add_option("-e", "--embl", action="store", dest="embl", help="input embl file name", type="string", default="", metavar="FILE")
	parser.add_option("-s", "--sites", action="store", dest="sites", help="input site list file name", type="string", default="", metavar="FILE")
	
	return parser.parse_args()



################################
# Check command line arguments #
################################

def check_input_validity(options, args):

	

	if options.embl=='':
		DoError('No embl file (-e) selected')
	elif not os.path.isfile(options.embl):
		DoError('Cannot find embl file')
	
	if options.sites=='':
		DoError('No site list file (-s) selected')
	elif not os.path.isfile(options.sites):
		DoError('Cannot find site list file')
		
	if options.outputfile=='':
		DoError('No output file (-o) selected')
		
	overwrite=True
	while not overwrite and os.path.isfile(options.outputfile):
		outopt=""
		outopt=raw_input('\nOutput file '+options.outputfile+' already exists.\n\nWould you like to overwrite (o), choose a new output file name (n) or quit (Q): ')
		if outopt=='Q':
			sys.exit()
		elif outopt=="o":
			overwrite=True
		elif outopt=="n":
			options.outputfile=raw_input('Enter a new output file name: ')
			
	

def revcomp(string):
	dict={"A":"T", "T":"A", "C":"G", "G":"C"}
	outlist=[]
	for base in string:
		outlist.append(dict[base])
	
	return ''.join(outlist)


			
################
# Main program #
################		

if __name__ == "__main__":


	#Get command line arguments

	(options, args) = get_user_options()
	check_input_validity(options, args)
	
	try:
		sitelines=open(options.sites, "rU")
	except StandardError:
		DoError("Failed to open site list file")
		
	sites=[]
	for line in sitelines:
		try:
			sites.append([int(line.strip().split()[0])]+line.strip().split()[1:])
		except StandardError:
			DoError("Sites file must contain integers in the first column")
	
	sites.sort()
	siteinfo={}
	for site in sites:
		siteinfo[site[0]]={}
	
	try:
		emblrecord=open_annotation(options.embl)
	except (StandardError, SimonError):
		DoError("Cannot open annotation file "+options.embl+" please check the format")
	
	added=[]
	x=0
	for feature in emblrecord.features:
		if feature.type=="CDS":
			
			#print feature.location.nofuzzy_start, feature.location.nofuzzy_end, feature.location.start, feature.location.end
			while x<len(sites) and int(sites[x][0])<feature.location.nofuzzy_start:
				
				site=int(sites[x][0])
				if sites[x][2]=="SNP":
					x+=1
					continue
				siteinfo[site]["location"]=str(site)+":"+sites[x][1]
				siteinfo[site]["start"]=site
				siteinfo[site]["end"]=int(sites[x][1])
				siteinfo[site]["type"]=sites[x][2]
				siteinfo[site]["SNP_type"]="Intergenic"
				if len(sites[x])>4:
					for y in xrange(4,len(sites[x])):
						name="additional_"+str(y)
						siteinfo[site][name]=sites[x][y]
						if not name in added:
							added.append(name)
				if siteinfo[site]["type"]!="SNP":
					#print siteinfo[site]
					if siteinfo[site]["type"]=="insertion":
						siteinfo[site]["start"]+=1
						siteinfo[site]["strain_base"]=sites[x][3].upper()
					elif siteinfo[site]["type"]=="deletion":
						siteinfo[site]["ref_base"]=sites[x][3].upper()
					
					x+=1
					continue
				siteinfo[site]["ref_base"]=str(emblrecord.seq[site-1:int(sites[x][1])])
				siteinfo[site]["strain_base"]=sites[x][3].upper()
#				siteinfo[site]["gene"]=None
#				siteinfo[site]["aa"]=None
#				siteinfo[site]["strand"]=None
				#print siteinfo[site]
				x+=1
			while x<len(sites) and int(sites[x][0])<=feature.location.nofuzzy_end:
				site=int(sites[x][0])
				if sites[x][2]=="SNP":
					x+=1
					continue
				siteinfo[site]["location"]=str(site)+":"+sites[x][1]
				siteinfo[site]["start"]=site
				siteinfo[site]["end"]=int(sites[x][1])
				siteinfo[site]["type"]=sites[x][2]
				if len(sites[x])>4:
					for y in xrange(4,len(sites[x])):
						name="additional_"+str(y)
						siteinfo[site][name]=sites[x][y]
						if not name in added:
							added.append(name)
				
				siteinfo[site]["CDS"]="Unknown"
				for f in ["CDS", "systematic_id", "locus_tag"]:
					if f in feature.qualifiers:
						siteinfo[site]["CDS"]=feature.qualifiers[f]
						break
				
				siteinfo[site]["strand"]=feature.strand
				if "product" in feature.qualifiers:
					siteinfo[site]["product"]=feature.qualifiers["product"]
				if siteinfo[site]["type"]!="SNP":
					#print siteinfo[site]
					if siteinfo[site]["type"]=="insertion":
						siteinfo[site]["start"]+=1
						siteinfo[site]["strain_base"]=sites[x][3].upper()
					elif siteinfo[site]["type"]=="deletion":
						siteinfo[site]["ref_base"]=sites[x][3].upper()
						
					x+=1
					continue
				
				if "pseudo" in feature.qualifiers:
					siteinfo[site]["pseudo"]=True
				else:
					siteinfo[site]["pseudo"]=False
				
				
				
				if feature.strand==-1:
					nucseq=str(emblrecord.seq[feature.location.nofuzzy_start:feature.location.nofuzzy_end].reverse_complement())
					if not siteinfo[site]["pseudo"]:
						aaseq=str(emblrecord.seq[feature.location.nofuzzy_start:feature.location.nofuzzy_end].reverse_complement().translate(table="Bacterial", cds=True, to_stop=False))
				else:
					nucseq=str(emblrecord.seq[feature.location.nofuzzy_start:feature.location.nofuzzy_end])
					if not siteinfo[site]["pseudo"]:
						try:
							aaseq=str(emblrecord.seq[feature.location.nofuzzy_start:feature.location.nofuzzy_end].translate(table="Bacterial", cds=True, to_stop=False))
						except CodonTable.TranslationError:
							x+=1
							continue
							
				
				siteinfo[site]["ref_base"]=str(emblrecord.seq[site-1:int(sites[x][1])])
				if siteinfo[site]["strand"]==-1:
					siteinfo[site]["CDS_location"]=feature.location.nofuzzy_end-site
				else:
					siteinfo[site]["CDS_location"]=site-1-feature.location.nofuzzy_start
				if not siteinfo[site]["pseudo"]:
#					print nucseq
#					print aaseq
					if siteinfo[site]["strand"]==-1:
						siteinfo[site]["CDS_aa_location"]=(feature.location.nofuzzy_end-site)/3
						siteinfo[site]["aa_base"]=(int(((float(feature.location.nofuzzy_end-site)/3)-((feature.location.nofuzzy_end-site)/3))*3.1))+1
						try:
							siteinfo[site]["ref_aa"]=aaseq[((feature.location.nofuzzy_end-site)/3)]
						except StandardError:
							siteinfo[site]["ref_aa"]="*"
						siteinfo[site]["ref_nucs"]=nucseq[(feature.location.nofuzzy_end-site)-(siteinfo[site]["aa_base"]-1):((feature.location.nofuzzy_end-site)-(siteinfo[site]["aa_base"]-1))+3]
					else:
						siteinfo[site]["CDS_aa_location"]=(site-1-feature.location.nofuzzy_start)/3
						siteinfo[site]["aa_base"]=(int(((float(site-1-feature.location.nofuzzy_start)/3)-((site-1-feature.location.nofuzzy_start)/3))*3.1))+1
						try:
							siteinfo[site]["ref_aa"]=aaseq[((site-1-feature.location.nofuzzy_start)/3)]
						except StandardError:
							siteinfo[site]["ref_aa"]="*"
						siteinfo[site]["ref_nucs"]=nucseq[(site-1-feature.location.nofuzzy_start)-(siteinfo[site]["aa_base"]-1):((site-1-feature.location.nofuzzy_start)-(siteinfo[site]["aa_base"]-1))+3]
					if siteinfo[site]["strand"]==-1 and nucseq[siteinfo[site]["CDS_location"]]!=revcomp(siteinfo[site]["ref_base"]):
						print Seq(nucseq[siteinfo[site]["CDS_location"]:siteinfo[site]["CDS_location"]+3]).translate(table="Bacterial", cds=False), Seq(nucseq[siteinfo[site]["CDS_location"]-1:siteinfo[site]["CDS_location"]+2]).translate(table="Bacterial", cds=False), Seq(nucseq[siteinfo[site]["CDS_location"]-2:siteinfo[site]["CDS_location"]+1]).translate(table="Bacterial", cds=False)
						print siteinfo[site]["CDS_aa_location"], siteinfo[site]["CDS_location"], site, feature.location.nofuzzy_start, feature.location.nofuzzy_end
						print siteinfo[site]["ref_base"], nucseq[siteinfo[site]["CDS_location"]]
						print aaseq[siteinfo[site]["CDS_aa_location"]]
						print siteinfo[site]["aa_base"]
						sys.exit()
						
					
					
					
#					print (site-feature.location.nofuzzy_start)-(siteinfo[site]["aa_base"]-1), ((site-feature.location.nofuzzy_start)-(siteinfo[site]["aa_base"]-1))+3
					
				else:
					siteinfo[site]["SNP_type"]="Pseudogene"
					
				if siteinfo[site]["type"]=="SNP":
					sitebase=sites[x][3].upper()
					if not sitebase in ["A","C","G","T"]:
						break
					siteinfo[site]["strain_base"]=sitebase
					if not siteinfo[site]["pseudo"]:
						try:
							strain_nucs=[siteinfo[site]["ref_nucs"][0].lower(),siteinfo[site]["ref_nucs"][1].lower(),siteinfo[site]["ref_nucs"][2].lower()]
						except IndexError:
							print siteinfo[site]["ref_nucs"], nucseq, aaseq
							
						if siteinfo[site]["strand"]==-1:
							strain_nucs[siteinfo[site]["aa_base"]-1]=revcomp(sitebase)
						else:
							strain_nucs[siteinfo[site]["aa_base"]-1]=sitebase
						siteinfo[site]["strain_nucs"]=''.join(strain_nucs)
						
						
						if siteinfo[site]["CDS_location"]<3:
							gene = Seq(siteinfo[site]["strain_nucs"])
							try:
								siteinfo[site]["strain_aa"]=str(gene.translate(table="Bacterial", cds=True))
							except CodonTable.TranslationError:
								gene = Seq("ATG"+siteinfo[site]["strain_nucs"])
								siteinfo[site]["strain_aa"]=str(gene.translate(table="Bacterial", cds=False))[1]
						else:
							gene = Seq("ATG"+siteinfo[site]["strain_nucs"])
							siteinfo[site]["strain_aa"]=str(gene.translate(table="Bacterial", cds=False))[1]
							
						
						if siteinfo[site]["strain_aa"]!=siteinfo[site]["ref_aa"]:
							siteinfo[site]["SNP_type"]="Nonsynonymous"
						else:
							siteinfo[site]["SNP_type"]="Synonymous"
				
				#print siteinfo[site]
				
				
				x+=1
				
			
	while x<len(sites) and int(sites[x][0])<feature.location.nofuzzy_start:
		
		site=int(sites[x][0])
		siteinfo[site]["location"]=str(site)+":"+sites[x][1]
		siteinfo[site]["type"]=sites[x][2]
		siteinfo[site]["SNP_type"]="Intergenic"
#				siteinfo[site]["gene"]=None
#				siteinfo[site]["aa"]=None
#				siteinfo[site]["strand"]=None
		print siteinfo[site]
		x+=1
	
	output=open(options.outputfile, "w")
	
	for x in xrange(len(sites)):
		site=int(sites[x][0])
		outlist=[]
		for datum in ["start", "type","ref_base","strain_base", "CDS", "strand", "product", "CDS_location", "SNP_type", "ref_nucs", "strain_nucs", "ref_aa", "strain_aa"]+added:
			if datum in siteinfo[site]:
				outlist.append(str(siteinfo[site][datum]))
			else:
				outlist.append("-")
		print >> output,  '\t'.join(outlist)
