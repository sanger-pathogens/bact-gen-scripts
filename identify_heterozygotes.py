#!/usr/bin/env python

#Create plot of dissimilarity between pairs of strains in an alignment



#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------
import string, re
import os, sys
import dendropy
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/lib/python2.7/site-packages/fisher-0.1.4-py2.7-linux-x86_64.egg/', '/nfs/users/nfs_s/sh16/lib/python2.5/site-packages/fisher-0.1.4-py2.7-linux-x86_64.egg/', '/nfs/users/nfs_s/sh16/scripts/modules/']))
from fisher import pvalue
from optparse import OptionParser
from Bio.Align import AlignInfo
#sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_general import *
from Si_SeqIO import *

##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options] args"
	parser = OptionParser(usage=usage)

	parser.add_option("-v", "--vcf", action="store", dest="vcf", help="SNP_sites output vcf file name.", default="", metavar="FILE")
	parser.add_option("-b", "--bcf", action="store", dest="bcf", help="bcf file name of test isolate.", default="", metavar="FILE")
	parser.add_option("-t", "--tree", action="store", dest="tree", help="tree file name", default="", metavar="FILE")
	parser.add_option("-i", "--id", action="store", dest="id", help="Name of test isolate.", default="", metavar="STRING")
	parser.add_option("-n", "--Nproportion", action="store", dest="Nproportion", help="maximum proportion of Ns to allow in isolates at a site for it to be included in the test (i.e. ignore any sites with > than this proportion of Ns) [default=%default]", default=0.05, type="float")
	parser.add_option("-p", "--proportion", action="store", dest="proportion", help="minimum proportion of mapped reads to allow to include an allele in the test [default=%default]", default=0.05, type="float")
	parser.add_option("-c", "--count", action="store", dest="count", help="minimum number of mapped reads to allow to include an allele in the test [default=%default]", default=5, type="int")
	parser.add_option("-V", "--verbose", action="store_true", dest="verbose", help="Be verbose", default=False)
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


#	if options.embl!="" and options.reference=='':
#		DoError('No reference selected! If you give an embl file you must specify which sequence it is linked to (i.e. the reference)')
	if options.vcf=='':
		DoError('No vcf file selected!')
	elif not os.path.isfile(options.vcf):
		DoError('Cannot find file '+options.vcf+'!')
	if options.tree=='':
		DoError('No tree file selected!')
	elif not os.path.isfile(options.tree):
		DoError('Cannot find file '+options.tree+'!')
	if options.Nproportion>1 or options.Nproportion<0:
		DoError('Maximum proportion of Ns must be between 0 and 1!')
	if options.proportion>1 or options.proportion<0:
		DoError('Minimum proportion of mapped reads must be between 0 and 1!')
	if options.count<0:
		DoError('Minimum number of mapped reads must be greater than or equal to 0!')
	
	if options.bcf=="":
		DoError('No bcf file selected!')
	
	if not os.path.isfile(options.bcf):
		DoError('Cannot find file '+arg+'!')
	
	if options.id=="":
		DoError("No taxon name provided")
		

gap_and_missing=set(["-", "N", "?"])
missing_set=([ "N", "?"])

def chisquarepvalue(taxon,taxonb,taxonc):
	
	taxonseq=alignment[taxon]
	taxonbseq=alignment[taxonb]
	taxoncseq=alignment[taxonc]
	
	diffN=0.0
	sameN=0.0
	diffnotN=0.0
	samenotN=0.0
	Ncount=0.0
	diffNpositions=set([])
	
	for x in xrange(len(taxonseq)):
		if x in toremove:
			continue
		if taxonseq[x] in missing_set:
			Ncount+=1
		if taxonbseq[x] in gap_and_missing or taxoncseq[x] in gap_and_missing:
			continue
		elif taxonbseq[x]==taxoncseq[x]:
			if taxonseq[x] in missing_set:
				sameN+=1
			elif taxonseq[x] not in missing_set:
				samenotN+=1
		elif taxonbseq[x]!=taxoncseq[x]:
			if taxonseq[x] in missing_set:
				diffN+=1
				diffNpositions.add(x)
			elif taxonseq[x] not in missing_set:
				diffnotN+=1
	
	#matrix=[[sameN, samenotN],[diffN, diffnotN]]
	p = pvalue(diffN, diffnotN, sameN, samenotN)
#	if p.right_tail<(0.05/tot):
#		print taxon, taxonb, taxonc, Ncount, sameN, diffN, samenotN, diffnotN, p
#	print p.right_tail
	return p.right_tail	, diffNpositions

def parse_snp_sites_vcf(handle):
	#Read the vcf file
	if options.verbose:
		print "Reading snp_sites vcf file"
	headings=[]
	data={}
	taxon_ids={}
	taxa=[]
	poscount=0
	percent_id={}
	for line in open(handle, "rU"):
		line=line.strip()
		if len(line)<2:
			continue
		elif line[0:2]=="##":
			words=line[2:].split("=")
			if words[0]=="contig":
				data[words[2].split(",")[0]]={}
		elif line[0]=="#":
			headings=line.split()
			for i, name in enumerate(headings[9:]):
				taxon_ids[i]=name
				taxa.append(name)
				percent_id[i]={}
				for j, nameb in enumerate(headings[9+i+1:]):
					percent_id[i][i+j]=0.0
		else:
			words=line.split()
			
			CHROM=words[0]
			POS=int(words[1])
#			if POS>10000:
#				break
			if not CHROM in data:
				DoError("CHROM ("+CHROM+") missing from data dictionary")
			data[CHROM][POS]={}
			ref_alleles=words[3].split(",")
			alt_alleles=words[4].split(",")
			alleles=ref_alleles+alt_alleles
			poscount+=1
			for allele in alleles:
				data[CHROM][POS][allele]=[]
			
			for i, isolate in enumerate(words[9:]):
				data[CHROM][POS][alleles[int(isolate)]].append(i)
				for j, isolateb in enumerate(words[9+i+1:]):
					if alleles[int(isolate)]==alleles[int(isolateb)] and alleles[int(isolate)] not in [".", "*"]:
						percent_id[i][i+j]+=1
			
	
	
	for i in percent_id:
		for j in percent_id[i]:
			percent_id[i][j]=(percent_id[i][j]/poscount)*100
			if not j in percent_id:
				percent_id[j]={}
			percent_id[j][i]=(percent_id[i][j]/poscount)*100
			
	if options.verbose:
		print "Found "+str(poscount)+" sites"
	return taxa, taxon_ids, data, percent_id
	

def remove_dodgy_sites(data):
	#"Remove sites with more than the specified proportion of Ns"
	if options.verbose:
		print ("Removing sites with a proportion of Ns greater than "+str(options.Nproportion))
	toremove=[]
	for chromosome in data:
		for position in data[chromosome]:
			total_taxa=0.0
			Ns=0.0
			for allele in data[chromosome][position]:
				total_taxa+=len(data[chromosome][position][allele])
			if "*" in data[chromosome][position]:
				Ns=float(len(data[chromosome][position]["*"]))
			N_proportion=Ns/total_taxa
			if N_proportion>options.Nproportion:
				toremove.append([chromosome, position])
	
	for chromosome, position in toremove:
		del data[chromosome][position]
	
	if options.verbose:
		print "Removed "+str(len(toremove))+" sites"
	
	return data


def open_bcf(filename):
	if options.verbose:
		print "Trying to open bcf file"
	fileopen=False
#	try:
#		bcffile=open(filename, "rU")
#		fileopen=True
#	except StandardError:
#		print "Cannot open bcf file as vcf"
#	if not fileopen:
#		try:
#			bcffile=os.popen("bcftools view "+filename)
#			fileopen=True
#		except StandardError:
#			print "Cannot open bcf file with default bcftools"
#	if not fileopen:
	try:
		bcffile=os.popen("bcftools-1.2 view "+filename)
		fileopen=True
	except StandardError:
		if options.verbose:
			print "Cannot open bcf file with bcftools-1.2"
	if not fileopen:
		DoError("Failed to open bcf file "+filename)
	
	return bcffile

def parse_bcf(handle, vcf_data):
	if options.verbose:
		print "Parsing bcf file"
	data={}
	for line in handle:
		line=line.strip()
		if len(line)<2:
			continue
		elif line[0:2]=="##":
			words=line[2:].split("=")
			if words[0]=="contig":
				data[words[2].split(",")[0]]=[]
		elif line[0]=="#":
			continue
		else:
			words=line.split()
			CHROM=words[0]
			if CHROM in vcf_data:
				POS=int(words[1])
				if POS in vcf_data[CHROM]:
					REF=words[3]
					ALT=words[4]
					if len(REF)>1 or len(ALT)>1:
						continue
					INFO=words[7].split(";")
					if "INDEL" in INFO:
						#print "skipping INDEL", POS
						continue
					founddp4=False
					for i in INFO:
						var=i.split("=")[0]
						val=i.split("=")[1]
						if var=="DP4":
							dp4vals=map(int,val.split(","))
							REFcount=dp4vals[0]+dp4vals[1]
							ALTcount=dp4vals[2]+dp4vals[3]
							founddp4=True
					if founddp4 and REFcount+ALTcount>0:
						if REFcount>=ALTcount:
							maj_allele=REF
							maj_allele_count=REFcount
							min_allele=ALT
							min_allele_count=ALTcount
						else:
							min_allele=REF
							min_allele_count=REFcount
							maj_allele=ALT
							maj_allele_count=ALTcount
						data[CHROM].append([POS, maj_allele, min_allele, maj_allele_count, min_allele_count])
	
	return data
					
			
	


################
# Main program #
################		

if __name__ == "__main__":
	

	
	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)
	
	
	#Read the vcf file
	
	taxon_list, taxon_id_dictionary, vcf_data, percent_id=parse_snp_sites_vcf(options.vcf)
	tid=-1
	for t, taxon in enumerate(taxon_list):
		if taxon==options.id:
			tid=t
	if tid==-1:
		DoError("Specified name not in vcf")
	
	vcf_data=remove_dodgy_sites(vcf_data)
	
	
	
	bcf_handle=open_bcf(options.bcf)
	bcf_data=parse_bcf(bcf_handle, vcf_data)
	
	maj_allele_identities={}
	min_allele_identities={}
	maj_allele_match_count=0.0
	min_allele_match_count=0.0
	for x in xrange(len(taxon_list)):
		maj_allele_identities[x]=0.0
		min_allele_identities[x]=0.0
	
	for chromosome in bcf_data:
		for i in bcf_data[chromosome]:
			if i[3]>options.count:
				if i[1] in vcf_data[chromosome][i[0]]:
					maj_allele_match_count+=1
					for taxon in vcf_data[chromosome][i[0]][i[1]]:
						maj_allele_identities[taxon]+=1
			if float(i[4])/(i[3]+i[4])>options.proportion and i[4]>options.count:
				if i[2] in vcf_data[chromosome][i[0]]:
					min_allele_match_count+=1
				#print i, float(i[4])/(i[3]+i[4])
					for taxon in vcf_data[chromosome][i[0]][i[2]]:
						min_allele_identities[taxon]+=1
			#elif i[3]>options.count:
			#	min_allele_match_count+=1
			#	for taxon in vcf_data[chromosome][i[0]][i[1]]:
			#		min_allele_identities[taxon]+=1
	
	maxmax=0
	topmaj=""
	topmajid=""
	for taxon in maj_allele_identities:
		if maj_allele_identities[taxon]>maxmax and taxon_list[taxon]!=options.id:
			maxmax=maj_allele_identities[taxon]
			topmaj=taxon_list[taxon]
			topmajid=taxon
	maxmin=0
	topmin=""
	topminid=""
	for taxon in min_allele_identities:
		if min_allele_identities[taxon]>maxmin and taxon_list[taxon]!=options.id:
			maxmin=min_allele_identities[taxon]
			topmin=taxon_list[taxon]
			topminid=taxon
	
	
	if maj_allele_match_count==0:
		majval=0
	else:
		majval=maxmax/maj_allele_match_count
		
	if min_allele_match_count==0:
		minval=0
	else:
		minval=maxmax/min_allele_match_count
		
		
	if topmaj=="" or topmin=="":
		majminperc="-"
	else:
		majminperc=percent_id[topmajid][topminid]
		
	if topmaj=="" or tid=="":
		majtperc="-"
	else:
		majtperc=percent_id[topmajid][tid]
	if topmin=="" or tid=="":
		mintperc="-"
	else:
		mintperc=percent_id[topminid][tid]
	
	print ','.join(map(str, [options.id, topmaj, maxmax, majval, topmin, maxmin, minval, majtperc, mintperc, majminperc]))
	sys.exit()
	
	
		
		
		

	
	
	
#	for x in toremove[::-1]:
#		for seq in seqnames:
#			alignment[seq]=alignment[seq][:x]+alignment[seq][x+1:]
#	alnlen=len(alignment[name])
#	
	
	for x, taxon in enumerate(seqnames):
		pairwise_distances[taxon]={}
		taxonseq=""
		print x+1, "\r",
		sys.stdout.flush()
		try: taxonseq=alignment[taxon]
		except StandardError:
			print "Cannot find ", taxon
			continue
#		print x
		for taxonb in seqnames[x+1:]:
						
			if taxon==taxonb:
				continue
			
			try: taxonbseq=alignment[taxonb]
			except StandardError:
				print "Cannot find ", taxonb
				continue
			
			count=0
			ncount=0			
			
			if len(taxonseq)==len(taxonbseq):
			
				for y, base in enumerate(taxonseq):
					if y in toremove:
						continue
					if base!=taxonbseq[y] and base not in gap_and_missing and taxonbseq[y] not in gap_and_missing:
						count+=1
						
						#print y+1, base, taxonbseq[y]
					elif base=="N":
						ncount+=1
						
			else:
				print taxon, "and", taxonb, "are not the same length"
			pairwise_distances[taxon][taxonb]=count
			pairwise_distances[taxonb][taxon]=count
			tot+=1
			#print taxon, taxonb, count, ncount

#	print "\n"
	
	significant={}
	Nlocsets={}
	
	for x, taxon in enumerate(seqnames):
		print x+1, "\r",
		sys.stdout.flush()
		for y, taxonb in enumerate(seqnames[x+1:]):
			
			sys.stdout.flush()
			for taxonc in seqnames[x+y+2:]:
				if pairwise_distances[taxon][taxonb]+pairwise_distances[taxon][taxonc]<pairwise_distances[taxonb][taxonc]:
					cspvalue, Nlocset=chisquarepvalue(taxon,taxonb,taxonc)
					if cspvalue<(0.05/tot):
						if not taxon in Nlocsets:
							Nlocsets[taxon]={}
						if not taxonb in Nlocsets[taxon]:
							Nlocsets[taxon][taxonb]={}
						if not taxonc in Nlocsets[taxon][taxonb]:
							Nlocsets[taxon][taxonb][taxonc]=Nlocset
					
						if not taxon in significant:
							significant[taxon]={}
						if not taxonb in significant[taxon]:
							significant[taxon][taxonb]=0
						significant[taxon][taxonb]+=1
						if not taxonc in significant[taxon]:
							significant[taxon][taxonc]=0
						significant[taxon][taxonc]+=1
				elif pairwise_distances[taxon][taxonb]+pairwise_distances[taxonb][taxonc]<pairwise_distances[taxon][taxonc]:
					cspvalue, Nlocset=chisquarepvalue(taxonb,taxon,taxonc)
					if cspvalue<(0.05/tot):
						if not taxonb in Nlocsets:
							Nlocsets[taxonb]={}
						if not taxon in Nlocsets[taxonb]:
							Nlocsets[taxonb][taxon]={}
						if not taxonc in Nlocsets[taxonb][taxon]:
							Nlocsets[taxonb][taxon][taxonc]=Nlocset
							
						if not taxonb in significant:
							significant[taxonb]={}
						if not taxon in significant[taxonb]:
							significant[taxonb][taxon]=0
						significant[taxonb][taxon]+=1
						if not taxonc in significant[taxonb]:
							significant[taxonb][taxonc]=0
						significant[taxonb][taxonc]+=1
				elif pairwise_distances[taxon][taxonc]+pairwise_distances[taxonb][taxonc]<pairwise_distances[taxon][taxonb]:
					cspvalue, Nlocset=chisquarepvalue(taxonc,taxonb,taxon)
					if cspvalue<(0.05/tot):
						if not taxonc in Nlocsets:
							Nlocsets[taxonc]={}
						if not taxonb in Nlocsets[taxonc]:
							Nlocsets[taxonc][taxon]={}
						if not taxonc in Nlocsets[taxonc][taxon]:
							Nlocsets[taxonc][taxon][taxonb]=Nlocset
							
						if not taxonc in significant:
							significant[taxonc]={}
						if not taxonb in significant[taxonc]:
							significant[taxonc][taxonb]=0
						significant[taxonc][taxonb]+=1
						if not taxon in significant[taxonc]:
							significant[taxonc][taxon]=0
						significant[taxonc][taxon]+=1
#				else:
#					print >> outfile, taxonb+","+"2"
#					print >> outfile, taxonc+","+"2"
	print
	logout=open("heterozygotes.log","w")
	print >> logout, "Taxon", "Number of significant comparisons"
	for taxon in seqnames:
		if taxon in significant:
			
			outfile=open("test.csv","w")
			print >> outfile, "name,type,number:c:1:"+str(len(seqnames)-1)
			print >> outfile, taxon+","+"1,"+"-"
			maxnum=0
			for taxonb in significant[taxon]:
				print >> outfile, taxonb+","+"-,"+str(significant[taxon][taxonb])
				if significant[taxon][taxonb]>maxnum:
					maxnum=significant[taxon][taxonb]

				
			outfile.close()
			print >> logout, taxon, float(maxnum)/(len(seqnames)-1)
			os.system("/nfs/users/nfs_s/sh16/scripts/reportlabtest.py -t "+options.tree+" -m test.csv -C 2,3 -O portrait -M -L left -a 2 -o "+taxon+"_"+str(maxnum)+".pdf")
			os.system("rm -f  test.csv")
		else:
			print >> logout, taxon, 0.0
	
	
	for taxon in Nlocsets:
		Nlocs={}
		for taxonb in Nlocsets[taxon]:
			for taxonc in Nlocsets[taxon][taxonb]:
				for base in Nlocsets[taxon][taxonb][taxonc]:
					if not base in Nlocs:
						Nlocs[base]=[]
					Nlocs[base].append('-'.join([taxonb,taxonc]))
		tabout=open(taxon+"_Ns.tab", "w")
		for Nloc in Nlocs:
			print >> tabout, "FT   misc_feature    "+str(refpos[Nloc])
			print >> tabout, 'FT                   /taxa="'+', '.join(Nlocs[Nloc])+'"'
		tabout.close()
				
	
	sys.exit()
                                    
