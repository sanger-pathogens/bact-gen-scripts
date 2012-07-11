#!/usr/bin/env python

##################
# Import modules #
##################

import os, sys, string
from optparse import OptionParser
import pysam

##############################
# Get command line arguments #
##############################

def get_user_options():
	usage = "usage: %prog [options]"
	version="%prog 1.0. Written by Simon Harris, Wellcome Trust Sanger Institute, 2012"
	parser = OptionParser(usage=usage, version=version)
	
	#do not allow arguments to be interspersed. e.g. -a -b arg1 agr2. MUST be -a arg1 -b arg2.
	parser.disable_interspersed_args()
	

	parser.add_option("-c", "--contigs", action="store", dest="contigs", help="multifasta containing contigs to search in", default="", metavar="FILE")
	parser.add_option("-g", "--genes", action="store", dest="genes", help="multifasta containing genes to search for", default="", metavar="FILE")
	parser.add_option("-b", "--bamfile", action="store", dest="bamfile", help="bamfile of mapped genes", default="", metavar="FILE")
	parser.add_option("-o", "--output", action="store", dest="output", help="prefix for output files", default="", metavar="FILE")
	parser.add_option("-s", "--sts", action="store", dest="sts", help="st file", default="", metavar="FILE")
	

	return parser.parse_args()

(options, args)=get_user_options()

if not os.path.isfile(options.contigs):
	print "Could not find contigs file"
	sys.exit()	
if not os.path.isfile(options.bamfile):
	print "Could not find bam file"
	sys.exit()	
if not os.path.isfile(options.genes):
	print "Could not find genes file"
	sys.exit()
try:
	contigsfile=open(options.contigs, "rU").read()
except StandardError:
	print "Could not open contigs file"
	sys.exit()

contigs={}
genes_present=[]
for line in contigsfile.split(">")[1:]:
	bits=line.split("\n")
	contigs[bits[0].split()[0]]=''.join(bits[1:])

filename=options.bamfile


try:
	genesfile=open(options.genes, "rU").read()
except StandardError:
	print "Could not open contigs file"
	sys.exit()
	
genes={}
geneorder=[]
for line in genesfile.split(">")[1:]:
	bits=line.split("\n")
	geneorder.append(bits[0].split()[0])
	genes[bits[0].split()[0]]=''.join(bits[1:])
refstarts={}
refends={}
count=0
for x, ref in enumerate(geneorder):
	refstarts[ref]=count
	count+=len(genes[ref])
	refends[ref]=count


if filename.split(".")[-1]=="bam":
	samfile = pysam.Samfile( filename, "rb" )
elif filename.split(".")[-1]=="sam":
        samfile = pysam.Samfile( filename, "r" )
else:
	print filename+" not a bam file"
        sys.exit() 

refs=samfile.references
lengths=samfile.lengths

for read in samfile:
	if read.is_unmapped:
            continue
        start=read.pos
        readpos=0
        refpos=start
	insertions=0
        inslength=0
        deletions=0
        dellength=0
        SNPs=0
        clipped=0
        cliplength=0
        for cig in read.cigar:
            
            if cig[0]==0:
                for x in range(0,cig[1]):
                    if read.seq[readpos].upper()!=contigs[refs[read.tid]][refpos].upper():
                        SNPs+=1
			print read.qname, read.tid, readpos, refpos, read.seq[readpos].upper(), contigs[refs[read.tid]][refpos].upper()
                    readpos+=1
                    refpos+=1
            elif cig[0]==1:
                insertions+=1
                inslength+=cig[1]
                readpos+=cig[1]
            elif cig[0]==2:
                deletions+=1
                dellength+=cig[1]
                refpos+=cig[1]
            elif cig[0]==4:
                clipped+=1
                cliplength+=cig[1]
                readpos+=cig[1]
            else:
                print cig
        end=refpos
        genes_present.append([read.qname, refs[read.tid], start, end, len(read.seq), SNPs, insertions, deletions, clipped, inslength, dellength, cliplength])
        
        #print filename, read.qname, len(read.seq), SNPs, insertions, deletions, clipped, inslength, dellength, cliplength



bestgene=genes_present[0]
secondary_genes=[]
outputlines=[]
for gene in genes_present[1:]:
    
    #print gene

    if gene[1]==bestgene[1]:
        if gene[2]<bestgene[3]:
            #print "overlap"
            if float(gene[5]+gene[6]+gene[7])/gene[4] < float(bestgene[5]+bestgene[6]+bestgene[7])/bestgene[4]:
                secondary_genes.append(bestgene[0])
                bestgene=gene
            else:
                secondary_genes.append(gene[0])
        else:
            outputlines.append(bestgene+[', '.join(secondary_genes)])
	    #print outputlines
            bestgene=gene
            secondary_genes=[]
    else:
        outputlines.append(bestgene+[', '.join(secondary_genes)])
        bestgene=gene
        secondary_genes=[]

outputlines.append(bestgene+[', '.join(secondary_genes)])


stlines=open(options.sts).readlines()

stgenes=map(string.lower,stlines[0].strip().split(","))
if not "st" in stgenes:
	print "Couldn't find ST column in ST file"
	sys.exit()
column_to_gene={}
gene_to_column={}
gene_order=[]
stcolumn=-1
profile=[]
for x,st in enumerate(stgenes):
	column_to_gene[x]=st
	if st=="st":
		stcolumn=x
	else:
		gene_order.append(st)
		profile.append(0)
		gene_to_column[st]=len(gene_order)-1
		

if stcolumn!=0:
	print "Couldn't find ST column in first column of ST file"
	sys.exit()


profiles={}
for line in stlines[1:]:
	words=line.strip().split(",")[1:]
	if len(words)!=len(gene_order):
		continue
	else:
		profiles[','.join(words)]=line.strip().split(",")[0]

notes=[]
LV=0
if options.output!="":
	output=open(options.output+"_MLST.txt","w")
	#plotlines=[]
	#print >> output, "\t".join(["gene", "contig", "start", "end", "length", "SNPs", "No. insertions", "No. deletions", "No. clipped regions", "total insertion length", "total deletion length", "clipped length", "Secondary gene hits in same region"])
	
	for line in outputlines:
		#print >> output, '\t'.join(map(str,line))

		#plotlines.append([refstarts[line[0]],refends[line[0]],((float(line[4])-float(line[5]+line[6]+line[7]))/line[4])*100])

		if line[0][:4].lower()!="st" and  line[0][:4].lower() in stgenes:
			profile[gene_to_column[line[0][:4].lower()]]=line[0][4:]
			if int(line[5])!=0:
				notes.append(str(line[5])+" SNPs in "+line[0][:4])
				LV+=1
			elif int(line[6])!=0:
				notes.append(str(line[6])+" deletions in "+line[0][:4])
			elif int(line[7])!=0:
				notes.append(str(line[7])+" insertions in "+line[0][:4])
				 
	if LV==0 and len(notes)==0:
		LVout="perfect match"
	elif LV==1:
		LVout="single locus variant"
	elif LV==2:
		LVout="double locus variant"
	elif LV==3:
		LVout="triple locus variant"
	elif LV==4:
		LVout="4 locus variant"
	elif LV==5:
		LVout="5 locus variant"
	elif LV==6:
		LVout="6 locus variant"
	elif LV==7:
		LVout="7 locus variant"
	else:
		LVout="incorrect gene length"		
		
	if ','.join(map(str,profile)) in profiles:						
		print >> output, options.output+","+','.join(map(str,profile)), ",ST"+profiles[','.join(map(str,profile))], LVout+",", ";".join(notes)
	else:
		print >> output, options.output+","+','.join(map(str,profile))+",No ST match,"+";".join(notes)
		       
			
			

	output.close()
	
	
    
else:
	#print "\t".join(["gene", "contig", "start", "end", "length", "SNPs", "No. insertions", "No. deletions", "No. clipped regions", "total insertion length", "total deletion length", "clipped length", "Secondary gene hits in same region"])
	for line in outputlines:
		#print '\t'.join(map(str,line))
		if line[0][:4].lower()!="st" and  line[0][:4].lower() in stgenes:
			profile[gene_to_column[line[0][:4].lower()]]=line[0][4:]
			if int(line[5])!=0:
				notes.append(str(line[5])+" SNPs in "+line[0][:4])
				LV+=1
			elif int(line[6])!=0:
				notes.append(str(line[6])+" deletions in "+line[0][:4])
			elif int(line[7])!=0:
				notes.append(str(line[7])+" insertiontions in "+line[0][:4])			 
			
	if LV==0 and len(notes)==0:
		LVout="perfect match"
	elif LV==1:
		LVout="single locus variant"
	elif LV==2:
		LVout="double locus variant"
	elif LV==3:
		LVout="triple locus variant"
	elif LV==4:
		LVout="4 locus variant"
	elif LV==5:
		LVout="5 locus variant"
	elif LV==6:
		LVout="6 locus variant"
	elif LV==7:
		LVout="7 locus variant"
	else:
		LVout="incorrect gene length"
		
	if ','.join(map(str,profile)) in profiles:
		print options.output, ','.join(map(str,profile)), "ST"+profiles[','.join(map(str,profile))], LVout+".", ";".join(notes)			
	else:
		
		print options.output, ','.join(map(str,profile)), "No ST match", LVout+".", ";".join(notes)	

