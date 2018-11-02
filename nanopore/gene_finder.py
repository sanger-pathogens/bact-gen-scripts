#!/usr/bin/env python

import os, sys
from optparse import OptionParser, OptionGroup
sys.path.extend(map(os.path.abspath, ['/nfs/pathogen/sh16_scripts/modules/']))
from Si_general import *
import string
from random import *

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
	
	parser.add_option("-r", "--reads", action="store", dest="reads", help="Fasta file containing reads", default="", metavar="FILE")
	parser.add_option("-d", "--database", action="store", dest="database", help="Fasta file containing database of genes", default="", metavar="FILE")
	parser.add_option("-c", "--cluster_seqs", action="store", dest="clusterseqs", help="cd-hid-est cluster representative sequence fasta file", default="", metavar="STRING")
#	parser.add_option("-C", "--clusters", action="store", dest="clusters", help="cd-hid-est cluster file", default=False)
	parser.add_option("-o", "--output", action="store", dest="output", help="Prefix for output files", default="", metavar="STRING")
	parser.add_option("-t", "--threads", action="store", dest="threads", help="Number of threads for running glsearch", default=1, metavar="INT", type="int")
	parser.add_option("-e", "--blastevalue", action="store", dest="bevalue", help="evalue cutoff for preliminary BLAST search", default=0.001, metavar="FLOAT", type="float")
	parser.add_option("-E", "--glseqrchevale", action="store", dest="gevalue", help="evalue cutoff for preliminary glsearch search", default=1e-100, metavar="FLOAT", type="float")
	parser.add_option("-i", "--ID", action="store", dest="percentid", help="Percent ID cutoff for preliminary glsearch search. Default = 0. Default for assembly mode = 70", default=0, metavar="FLOAT", type="float")
	parser.add_option("-C", "--consensus", action="store", dest="consensus", help="Percent of reads possessing bae required to call consensus. Must be between 50 and 100.", default=50, metavar="FLOAT", type="float")
	parser.add_option("-s", "--strict", action="store_true", dest="strict", help="Stricter search where only the best match from each read for each gene is included in consensus", default=False)
	parser.add_option("-a", "--assembly", action="store_true", dest="assembly", help="Assembly mode when using the script on assemblies rather than reads, where multiple copies of genes should be treated independently", default=False)
	parser.add_option("-m", "--measure", action="store", type="choice", dest="measure", choices=["percentid", "bitscore", "evalue"], help="Measure to use to rank glsearch hits. Choose from bitscore, evalue and percentid [default= %default]", default="bitscore")
	parser.add_option("-R", "--restart", action="store_true", dest="restart", help="Attempt to restart previous analysis. Requires other command line options to be the same as previously", default=False)
#	parser.add_option("-t", "--tab", action="store_true", dest="tab", help="Create tab file showing aligned blocks", default=False)
#	parser.add_option("-s", "--supress", action="store_true", dest="supress", help="Supress creation of fasta of novel regions", default=False)
	
	
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


#	if options.processors>cpu_count():
#		DoError('You do not have '+str(options.processors)+' processors')
#
	if options.reads=='':
		DoError('No reads fasta file selected (-r)')
	elif not os.path.isfile(options.reads):
		DoError('Cannot find file '+options.reads)
	
	if options.database=='':
		DoError('No database fasta file selected (-d)')
	elif not os.path.isfile(options.database):
		DoError('Cannot find file '+options.database)
	
	if options.clusterseqs=='':
		DoError('No cd-hit cluster sequence fasta file selected (-c)')
	elif not os.path.isfile(options.clusterseqs):
		DoError('Cannot find file '+options.clusterseqs)
	options.clusters=options.clusterseqs+".clstr"
	if not os.path.isfile(options.clusters):
		DoError('Cannot find cluster file '+options.clusters+" from cd-hit-est")
	
	if options.output=="":
		DoError('No output prefix selected (-o)')
	
	if options.threads>32 or options.threads<1:
		DoError('Threads must be between 1 and 32')
	
	if options.percentid>100 or options.percentid<0:
		DoError('ID cutoff must be between 0 and 100')
	
	if options.assembly and not "-p" in sys.argv:
		print "Assembly mode. Setting default %ID cutoff to 70%"
		options.percentid=70
	
	if options.bevalue>10 or options.bevalue<0:
		DoError('BLAST evalue (-e) must be between 0 and 10')
	
	if options.gevalue>10 or options.gevalue<0:
		DoError('glsearch evalue (-E) must be between 0 and 10')
	
	if options.consensus>100 or options.consensus<50:
		DoError('Consensus percentage (-C) must be between 50 and 100.')
	
	return



##################################################
# Function to find a consensus from an alignment #
##################################################

def print_consensus(alignment_file, percent=50, min_matches=0, ref=""):
	
	seqs=[]
	inread=False
	for line in open(alignment_file):
		words=line.strip().split()
		if len(words)==0:
			continue
		if words[0][0]==">":
			if inread:
				seqs.append(''.join(sequence))
			name=words[0][1:]
			sequence=[]
			inread=True
		else:
			sequence.append(words[0].upper())
	if inread:
		seqs.append(''.join(sequence))
	
	#check they're all the same length
	seqnum=0
	for seq in seqs:
		if ref!="" and seq==ref:
			continue
		seqnum+=1
		if len(seq)!=len(seqs[0]):
			print "Expecting sequences in alignment to all be of the same length"
			return
	
	
	consensus=[]
	
	for basenum in xrange(0, len(seqs[0])):
		bases={'A':0, 'C':0, 'G':0, 'T':0, '-':0}
		
		for seq in seqs:
			if ref!="" and seq==ref:
				continue
			if seq[basenum] in bases:
				bases[seq[basenum]]+=1
				
		n=True
		
		for base in bases:
			if min_matches!=0:
				if base in ['A', 'C', 'G', 'T'] and (bases[base]>=float(seqnum)/(100.0/percent) or (bases[base]>=float(seqnum-bases["-"])/(100.0/percent) and bases[base]>min_matches)):
					consensus.append(base)
					n=False
			elif base in ['A', 'C', 'G', 'T'] and bases[base]>=float(seqnum)/(100.0/percent):
				consensus.append(base)
				n=False
			elif base=='-' and bases[base]>=float(seqnum)/(100.0/percent):
				n=False
				
		if n:
			consensus.append("N")
	
	
	return ''.join(consensus)
				


#######################################################################################
# Function to find best match to query sequence in an alignment and report some stats #
#######################################################################################

def alignment_stats(alignment_file, query="Consensus"):
	
	seqs={}
	inread=False
	for line in open(alignment_file):
		words=line.strip().split()
		if len(words)==0:
			continue
		if words[0][0]==">":
			if inread:
				seqs[name]=''.join(sequence)
			name=words[0][1:]
			sequence=[]
			inread=True
		else:
			sequence.append(words[0].upper())
	if inread:
		seqs[name]=''.join(sequence)
	
	if not query in seqs:
		print "Could not find query sequence in alignment file"
		print "Query =", query
		print "Sequences in aligment =", seqs.keys()
		return 1
	
	
	#check they're all the same length and create stats dictionary
	seq_stats={}
	for seq in seqs:
		if len(seqs[seq])!=len(seqs[query]):
			print "Expecting sequences in alignment to all be of the same length"
			return 1
		if seq!=query:
			seq_stats[seq]=0
	
	consensus=[]
	queryseq=seqs[query]
	for seq in seqs:
		if not seq==query:
			for basenum in xrange(0, len(seqs[query])):
				sb=seqs[seq][basenum]
				qb=seqs[query][basenum]
				if sb in ['A', 'C', 'G', 'T'] and sb==qb:
					seq_stats[seq]+=1
	
	best_match=0
	best_match_seq=""
	for seq in seq_stats:
		if seq_stats[seq]>best_match:
			best_match_seq=seq
			best_match=seq_stats[seq]

	SNPs=0
	base_matches=0
	gaps=0
	Ns=0
	insertions=0
	deletions=0
	indel=False
	inins=False
	for basenum in xrange(0, len(seqs[query])):
		sb=seqs[best_match_seq][basenum]
		qb=seqs[query][basenum]
		if qb in ['A', 'C', 'G', 'T'] and sb==qb:
			base_matches+=1
		elif qb in ['A', 'C', 'G', 'T'] and sb in ['A', 'C', 'G', 'T'] and sb!=qb:
			SNPs+=1
		elif qb=='-' and sb in ['A', 'C', 'G', 'T']:
			gaps+=1
			if not indel:
				indel=True
				deletions+=1
		elif sb=='-' and qb in ['A', 'C', 'G', 'T']:
			if not inins:
				inins=True
				insertions+=1
		elif qb not in ['A', 'C', 'G', 'T', '-']:
			Ns+=1
		
		
		
		if indel and qb!='-':
			indel=False
		elif inins and sb!='-':
			inins=False
	
	percentID=100*(float(base_matches)/(base_matches+SNPs))
	
	return [best_match_seq, len(seqs[query]), base_matches, SNPs, gaps, Ns, insertions, deletions, percentID]




################
# Main program #
################		

if __name__ == "__main__":

	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)
	
	
	chars = string.ascii_letters + string.digits
	tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))
	#tmpname="tmpZ5irP0sy"
	
	#could add code here to run ch-dit-est
	
	if options.restart and os.path.isfile(options.output+"_BLAST_matches.txt"):
		reads_passing=set([])
		genes_passing=set([])
		for line in open(options.output+"_BLAST_matches.txt"):
			words=line.strip().split()
			if len(words)==0:
				continue
			if words[0]!="#" and len(words)==12:
				query=words[0]
				subject=words[1]
				evalue=float(words[10])
				if evalue<options.bevalue:
					reads_passing.add(subject)
					genes_passing.add(query)
	else:
		print "Creating BLAST database...",
		sys.stdout.flush()
		os.system("formatdb -p F -i "+options.reads)
		print "Done"
		print "Running preliminary BLAST search...",
		sys.stdout.flush()
		os.system("blastall -p blastn -i "+options.clusterseqs+" -d "+options.reads+" -m 8 -e "+str(options.bevalue)+" -o "+tmpname+"_blast.out")
		print "Done"
	
		print "Parsing BLAST output...",
		sys.stdout.flush()
		#Read blast output and make list of reads and genes to include in fasta search
		reads_passing=set([])
		genes_passing=set([])
		output=open(options.output+"_BLAST_matches.txt","w")
		for line in open(tmpname+"_blast.out"):
			words=line.strip().split()
			if len(words)==0:
				continue
			if words[0]!="#" and len(words)==12:
				query=words[0]
				subject=words[1]
				evalue=float(words[10])
				if evalue<options.bevalue:
					reads_passing.add(subject)
					genes_passing.add(query)
					print >> output, line.strip()
		output.close()

	print "Found", len(reads_passing), "reads and", len(genes_passing), "clusters with matches"
	
	print "Creating files for glsearch...",
	sys.stdout.flush()
	#Extract reads which have been found with matches
	
	output=open(tmpname+"_reads.fasta", "w")
	inread=False
#	printed_warning=False
	new_2_old_name={}
	old_2_new_name={}
	seqcount=0
	for line in open(options.reads):
		words=line.strip().split()
		if len(words)==0:
			continue
		if words[0][0]==">":
			
			if inread and name in reads_passing:
				print >> output, ">"+old_2_new_name[name]
				print >> output, ''.join(sequence)
			name=words[0][1:]
			#to restart more recent runs
			
			
			seqcount+=1
			old_2_new_name[name]="Seq"+str(seqcount)
			new_2_old_name["Seq"+str(seqcount)]=name
			#to restart older runs
			
#			old_2_new_name[name]=name
#			new_2_old_name[name]=name

			
			sequence=[]
			inread=True
		else:
			sequence.append(line.strip())
	if inread and name in reads_passing:
		print >> output, ">"+old_2_new_name[name]
		print >> output, ''.join(sequence)
	
	output.close()
	
	#Extract genes which have been found with matches
	
	output=open(tmpname+"_genes.fasta", "w")
	inread=False
	gene_comments={}
	new_2_old_gene={}
	old_2_new_gene={}
	seqcount=0
	for line in open(options.database):
		words=line.strip().split()
		if len(words)==0:
			continue
		if words[0][0]==">":
			if inread and name in genes_passing:
				print >> output, ">"+old_2_new_gene[name]
				print >> output, ''.join(sequence)
			name=words[0][1:]
			
			if len(words)>1:
				gene_comments[name]=' '.join(words[1:])
			else:
				gene_comments[name]=''
#			if len(name)>50:
			seqcount+=1
			old_2_new_gene[name]="Gene"+str(seqcount)
			new_2_old_gene["Gene"+str(seqcount)]=name
			#to restart older runs
			
#			old_2_new_gene[name]=name
#			new_2_old_gene[name]=name

			
			sequence=[]
			inread=True
		else:
			sequence.append(words[0].upper())
	if inread and name in genes_passing:
		print >> output, ">"+old_2_new_gene[name]
		print >> output, ''.join(sequence)
	
	output.close()
	print "Done"
	
	
	if options.restart and os.path.isfile(options.output+"_glsearch_matches.txt"):
		print "Parsing glsearch output"
		sys.stdout.flush()
		matchdb={}
		for line in open(options.output+"_glsearch_matches.txt"):
			words=line.strip().split()
			if words[0]!="#" and len(words)==12:
				query=new_2_old_gene[words[0]]
				subject=new_2_old_name[words[1]]
				qstart=int(words[6])
				qend=int(words[7])
				sstart=int(words[8])
				send=int(words[9])
				evalue=float(words[10])
				bitscore=float(words[11])
				percentid=float(words[2])
				if qstart>qend:
					fr="r"
				else:
					fr="f"
				if evalue<options.gevalue and percentid>=options.percentid:
					if not subject in matchdb:
						matchdb[subject]={}
					if not query in matchdb[subject]:
						matchdb[subject][query]=[]
					if options.measure=="bitscore":
						matchdb[subject][query].append([bitscore, sstart, send, fr, evalue])
					elif options.measure=="percentid":
						matchdb[subject][query].append([percentid, sstart, send, fr, evalue])
					elif options.measure=="evalue":
						matchdb[subject][query].append([evalue, sstart, send, fr, evalue])
	else:
		print "Running glsearch to find global gene matches...",
		sys.stdout.flush()
		#run global search to find gene matches in reads
		os.system("glsearch -E "+str(options.gevalue)+" -T "+str(options.threads)+" -m 8C "+tmpname+"_genes.fasta "+tmpname+"_reads.fasta > "+tmpname+"_glsearch.out")
		print "Done"
		
		
		
		#sys.exit()
		
		print "Parsing glsearch output"
		sys.stdout.flush()
		#Read the fasta output and save matches for each query to a dictionary
		
		matchdb={}
		output=open(options.output+"_glsearch_matches.txt","w")
		for line in open(tmpname+"_glsearch.out"):
			words=line.strip().split()
			if words[0]!="#" and len(words)==12:
				query=new_2_old_gene[words[0]]
				subject=new_2_old_name[words[1]]
				qstart=int(words[6])
				qend=int(words[7])
				sstart=int(words[8])
				send=int(words[9])
				evalue=float(words[10])
				bitscore=float(words[11])
				percentid=float(words[2])
				if qstart>qend:
					fr="r"
				else:
					fr="f"
				if evalue<options.gevalue and percentid>=options.percentid:
					if not subject in matchdb:
						matchdb[subject]={}
					if not query in matchdb[subject]:
						matchdb[subject][query]=[]
					if options.measure=="bitscore":
						matchdb[subject][query].append([bitscore, sstart, send, fr, evalue])
					elif options.measure=="percentid":
						matchdb[subject][query].append([percentid, sstart, send, fr, evalue])
					elif options.measure=="evalue":
						matchdb[subject][query].append([evalue, sstart, send, fr, evalue])
					print >> output, line.strip()
		output.close()
	
	
	print "Removing overlapping matches"
	sys.stdout.flush()
	#Reduce matches that overlap to the best match decided by bitscore/%id/evalue
	
	querymatchdb={}
	
	for s in matchdb:
		matches=[]
		for q in matchdb[s]:
			matchdb[s][q].sort()
			matchdb[s][q].reverse()
			for x in xrange(len(matchdb[s][q])):
				matches.append([matchdb[s][q][x][0],q, matchdb[s][q][x][1], matchdb[s][q][x][2], matchdb[s][q][x][3], matchdb[s][q][x][4]])
		matches.sort()
		matches.reverse()
		
		unique_matches=[]
		
		for m in matches:
			keep=True
			for um in unique_matches:
				 if m[2]<um[3] and m[3]>um[2]:
				 	keep=False
				 	break
			if keep:
				unique_matches.append(m)
				
		for um in unique_matches:
			if not um[1] in querymatchdb:	
				querymatchdb[um[1]]={}
			if not s in querymatchdb[um[1]]:
				querymatchdb[um[1]][s]=[]
			querymatchdb[um[1]][s].append([um[2], um[3], um[4], um[5], um[0]])
	
	
	
	print "Creating cluster dictionary"
	sys.stdout.flush()
	#Read the cd-hit cluster file to create a cluster dictionary. Also save a dictionary linking the reference gene for each cluster to the cluster number
	
	clusters={}
	ref_to_cluster={}
	for line in open(options.clusters):
		words=line.strip().split()
		if words[0][0]==">":
			cluster=int(words[1])
			clusters[cluster]=[]
		else:
			gene=words[2][1:-3]
			clusters[cluster].append(gene)
			if words[-1]=="*" and gene in genes_passing:
				ref_to_cluster[gene]=cluster
	
	print "Extracting matching reads"
	sys.stdout.flush()
	#Extract reads which have been found with matches
	
	reads={}
	sequence=[]
	inread=False
	for line in open(options.reads):
		words=line.strip().split()
		if words[0][0]==">":
			if inread and name in matchdb:
				reads[name]=''.join(sequence)
			name=words[0][1:]
			sequence=[]
			inread=True
		else:
			sequence.append(line.strip())
	if inread and name in matchdb:
		reads[name]=''.join(sequence)
	
	
	print "Extracting matching gene sequences...",
	sys.stdout.flush()
	#Create sequence database for clusters
	cluster_ref_seqs={}
	for line in open(options.clusterseqs):
		words=line.strip().split()
		if words[0][0]==">":
			if inread and name in querymatchdb:
				cluster_ref_seqs[name]=''.join(sequence)
			name=words[0][1:]
			sequence=[]
			inread=True
		else:
			sequence.append(words[0].upper())
	if inread and name in querymatchdb:
		cluster_ref_seqs[name]=''.join(sequence)
	
	print "Found", len(cluster_ref_seqs), "clusters with matches"
	
	print "Aligning matches and creating consensuses"
	if options.assembly:
		print "Using assembly option whereby all matches will be treated independently"
	elif options.strict:
		print "Using strict option whereby only best match from each contig for each gene will be used in consensus"
	sys.stdout.flush()
	
	options.limit=50
	if not options.assembly:
		if options.limit>0:
			print "Limiting to", options.limit, "best matches in consensus alignment"
		sys.stdout.flush()
	
	if options.assembly:
		consensus_sequences={}
		ref_matches={}
		for ref in cluster_ref_seqs:
			consensus_sequences[ref]=[]
			ref_matches[ref]=[0,0]
			for match in querymatchdb[ref]:
				ref_matches[ref][0]+=1
				for x in xrange(len(querymatchdb[ref][match])):
					ref_matches[ref][1]+=1
					mymatch=[querymatchdb[ref][match][x][4], match, x, querymatchdb[ref][match][x][0], querymatchdb[ref][match][x][1], querymatchdb[ref][match][x][2]]
					if mymatch[5]=="f":
						consensus_sequences[ref].append(reads[mymatch[1]][mymatch[3]-1:mymatch[4]])
					else:
						consensus_sequences[ref].append(revcomp(reads[mymatch[1]][mymatch[3]-1:mymatch[4]]))
				
	else:
		#align each set of matched sequences and create 50% consensus
		consensus_sequences={}
		ref_matches={}
		for ref in cluster_ref_seqs:
			musclefile=open(tmpname+".fasta", "w")
			print "\t"+ref
			sys.stdout.flush()
	#		print >>  musclefile, ">"+ref
	#		print >>  musclefile, cluster_ref_seqs[ref]
			matchlist=[]
			ref_matches[ref]=[0,0]
			for match in querymatchdb[ref]:
				ref_matches[ref][0]+=1
				for x in xrange(len(querymatchdb[ref][match])):
					ref_matches[ref][1]+=1
					if x==0 or not options.strict:
						matchlist.append([querymatchdb[ref][match][x][4], match, x, querymatchdb[ref][match][x][0], querymatchdb[ref][match][x][1], querymatchdb[ref][match][x][2]])
	
			
			matchlist.sort()
			matchlist.reverse()
			
			
			for x, match in enumerate(matchlist):
				if x>=options.limit:
					break
				print >>  musclefile, ">"+match[1]
				if match[5]=="f":
					print >>  musclefile, reads[match[1]][match[3]-1:match[4]]
				else:
					print >>  musclefile, revcomp(reads[match[1]][match[3]-1:match[4]])
			
	#		for match in querymatchdb[ref]:
	#			ref_matches[ref][0]+=1
	#			for x in xrange(len(querymatchdb[ref][match])):
	#				ref_matches[ref][1]+=1
	#				if x==0 or not options.strict:#This if makes the script only use the top match for each read against each cluster for consensus calculation
	#					print >>  musclefile, ">"+match
	#					if querymatchdb[ref][match][x][2]=="f":
	#						print >>  musclefile, reads[match][querymatchdb[ref][match][x][0]:querymatchdb[ref][match][x][1]]
	#					else:
	#						print >>  musclefile, revcomp(reads[match][querymatchdb[ref][match][x][0]:querymatchdb[ref][match][x][1]])
			
			musclefile.close()
	#		os.system("seaview "+tmpname+".fasta")
			os.system("muscle -in "+tmpname+".fasta -out "+tmpname+".aln -gapopen -6 &> /dev/null")
	#		os.system("seaview "+tmpname+".aln")
			consensus_sequences[ref]=[print_consensus(tmpname+".aln", percent=options.consensus)]
			#sys.exit()
	
	csout=open(options.output+"_consensus_sequences.mfa", "w")
	output=open(options.output+"_hits.txt", "w")
	print >> output, '\t'.join(map(str,["Cluster number", "Cluster reference", "Comments", "# matching reads", "# matches in reads", "Best match gene in cluster", "Alignment length", "Matches", "SNPs", "Gaps", "Ns", "Insertions", "Deletions", "Percent ID of called bases"]))
	print "Aligning consensuses with clusters"
	#Align consensus sequence with cluster
	for ref in consensus_sequences:
		for x in xrange(len(consensus_sequences[ref])):
			tmpoutput=open(tmpname+".fasta","w")
			print "\t"+ref
			sys.stdout.flush()
			print >> tmpoutput, ">Consensus"
			print >> tmpoutput, consensus_sequences[ref][x]
			print >> csout, ">"+ref+"_match_consensus"
			print >> csout, consensus_sequences[ref][x]
			
			inread=False
			for line in open(options.database):
				words=line.strip().split()
				if len(words)==0:
					continue
				if words[0][0]==">":
					if inread and name in clusters[ref_to_cluster[ref]]:
						print >> tmpoutput, ">"+name
						print >> tmpoutput, ''.join(sequence)
					name=words[0][1:]
					sequence=[]
					inread=True
				else:
					sequence.append(words[0].upper())
			if inread and name in clusters[ref_to_cluster[ref]]:
				print >> tmpoutput, ">"+name
				print >> tmpoutput, ''.join(sequence)
			tmpoutput.close()
			
			os.system("muscle -in "+tmpname+".fasta -out "+tmpname+".aln &> /dev/null")
	#		os.system("seaview "+tmpname+".aln")
			stats=alignment_stats(tmpname+".aln")
	#		print stats
			print >> output, '\t'.join(map(str,[ref_to_cluster[ref]]+[ref]+[gene_comments[ref]]+ref_matches[ref]+stats))
		
		#read the alignment and print some stats
	output.close()
	csout.close()
	print "Cleaning up"
	os.system("rm -f "+tmpname+"*")
	print "All Done"
		
