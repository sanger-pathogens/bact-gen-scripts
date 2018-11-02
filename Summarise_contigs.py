#!/usr/bin/env python

import os, sys
sys.path.extend(map(os.path.abspath, ['/nfs/pathogen/sh16_scripts/modules/']))
from Si_SeqIO import *
from Si_general import *
from Bio import SeqIO
from Bio.SeqUtils import GC
import gzip





try:
	emblrecord=open_annotation(sys.argv[1])
except StandardError:
	DoError("Cannot open annotation file "+sys.argv[1]+" please check the format")


contig_locations={}
contig_hits={}
All_hits=[]
contignumber=0

for x, seqFeature in enumerate(emblrecord.features):
	if seqFeature.type=="fasta_record":
		contignumber+=1
		seqFeature.qualifiers["label"]="contig_"+str(contignumber)
		del seqFeature.qualifiers["note"]
		seqFeature.qualifiers["gc"]= GC(emblrecord.seq[int(seqFeature.location.nofuzzy_start):int(seqFeature.location.nofuzzy_end)])
		contig_locations["contig_"+str(contignumber)]=[int(seqFeature.location.nofuzzy_start),int(seqFeature.location.nofuzzy_end)]
		contig_hits["contig_"+str(contignumber)]={"total":0}


for x, seqFeature in enumerate(emblrecord.features):
	
	if seqFeature.type=="CDS":
		start=int(seqFeature.location.nofuzzy_start)
		end=int(seqFeature.location.nofuzzy_end)
		incontig="?"
		for contig in contig_locations.keys():
			if start>=contig_locations[contig][0] and end<=contig_locations[contig][1]:
				incontig=contig
				break
		
		if seqFeature.qualifiers.has_key("fasta_file"):
			filename=seqFeature.qualifiers["fasta_file"][0].split(":")[1]+".gz"
			fastafile=gzip.open(filename,"r")
			for line in fastafile:
				if len(line.split(":")[0]) and line.split(":")[0]=="The best scores are":
					line=fastafile.next()
					while len(line.split())>0 and float(line.strip().split()[-1])<1e-50:
	#						print line.strip()
	#						
						#if float(line.strip().split()[-1])<1e-5:
						mfetchout=os.popen('mfetch -f "tax org" '+line.split()[0])
						taxonomy=[]
						top_hit=[]
						for line in mfetchout:
							if line.split()[0]=="OS":
								top_hit.append(' '.join(line.split()[1:]))
							elif line.split()[0]=="OC":
								taxonomy.append(' '.join(line.split()[1:]))
						seqFeature.qualifiers["top_hit_taxonomy"]=' '.join(taxonomy)
						seqFeature.qualifiers["top_hit"]=' '.join(top_hit)
						if contig_hits.has_key(incontig):
							contig_hits[incontig]["total"]+=1
							for taxon in ' '.join(taxonomy).replace(".","").split(";")[:3]:
								taxon=taxon.strip()
								if contig_hits[incontig].has_key(taxon):
									contig_hits[incontig][taxon]+=1
								else:
									contig_hits[incontig][taxon]=1
									if not taxon in All_hits:
										All_hits.append(taxon)
						line=fastafile.next()
		
		seqFeature.qualifiers["contig"]=incontig
#		print seqFeature
#		print contig_hits[incontig]
#		print contig_hit_list[incontig]

#print All_hits


handle =open('.'.join(sys.argv[1].split(".")[:-1])+"_summary.txt","w")

print >> handle, "\t".join(["contig", "length", "GC%", "Total_hits"]+All_hits)

for seqFeature in emblrecord.features:
	outlist=[]
	if seqFeature.type=="fasta_record":
		outlist.append(seqFeature.qualifiers["label"])
		outlist.append(str(int(seqFeature.location.nofuzzy_end)-int(seqFeature.location.nofuzzy_start)))
		outlist.append(str(seqFeature.qualifiers["gc"]))
		if contig_hits.has_key(seqFeature.qualifiers["label"]):
			outlist.append(str(contig_hits[seqFeature.qualifiers["label"]]["total"]))
		else:
			outlist.append("0")
		for taxon in All_hits:
			if contig_hits.has_key(seqFeature.qualifiers["label"]) and contig_hits[seqFeature.qualifiers["label"]].has_key(taxon):
				outlist.append(str(100*(float(contig_hits[seqFeature.qualifiers["label"]][taxon])/float(contig_hits[seqFeature.qualifiers["label"]]["total"]))))
			else:
				outlist.append("0")
		print >> handle, "\t".join(outlist)
handle.close()
handle =open('.'.join(sys.argv[1].split(".")[:-1])+"_summary.gb","w")
SeqIO.write([emblrecord],handle, "gb")
handle.close()