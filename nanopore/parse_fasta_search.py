#!/usr/bin/env python

import os, sys


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


def print_consensus(alignment_file, percent=50):
	
	seqs=[]
	inread=False
	for line in open(alignment_file):
		words=line.strip().split()
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
	for seq in seqs:
		if len(seq)!=len(seqs[0]):
			print "Expecting sequences in alignment to all be of the same length"
			return
	
	
	consensus=[]
	for basenum in xrange(0, len(seqs[0])):
		bases={'A':0, 'C':0, 'G':0, 'T':0, '-':0}
		for seq in seqs:
			if seq[basenum] in bases:
				bases[seq[basenum]]+=1
				
		n=True
		for base in bases:
			if base in ['A', 'C', 'G', 'T'] and bases[base]>=float(len(seqs))/(100.0/percent):
				consensus.append(base)
				n=False
				
			elif base=='-' and bases[base]>=float(len(seqs))/(100.0/percent):
				n=False
				
		if n:
			consensus.append("N")
	
	return ''.join(consensus)
				




def alignment_stats(alignment_file, query="Consensus"):
	
	seqs={}
	inread=False
	for line in open(alignment_file):
		words=line.strip().split()
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

	
	return best_match_seq, (float(best_match)/len(seqs[query]))*100





print "Creating BLAST database"
os.system("formatdb -p F -i "+sys.argv[3])
print "Running preliminary BLAST search"
os.system("blastall -p blastn -i Curated_Database_20131022_clusters -d "+sys.argv[3]+" -m 8 -e 1e-5 -o tmp_blast.out")


print "Parsing BLAST output"
#Read blast output and make list of reads and genes to include in fasta search
reads_passing=set([])
genes_passing=set([])

for line in open("tmp_blast.out"):
	words=line.strip().split()
	if words[0]!="#" and len(words)==12:
		query=words[0]
		subject=words[1]
		evalue=float(words[10])
		if evalue<1e-10:
			reads_passing.add(subject)
			genes_passing.add(query)


print "Creating files for glsearch"
#Extract reads which have been found with matches

output=open("tmp_reads.fasta", "w")
inread=False
for line in open(sys.argv[3]):
	words=line.strip().split()
	if words[0][0]==">":
		if inread and name in reads_passing:
			print >> output, ">"+name
			print >> output, ''.join(sequence)
		name=words[0][1:]
		sequence=[]
		inread=True
	else:
		sequence.append(line)
if inread and name in reads_passing:
	print >> output, ">"+name
	print >> output, ''.join(sequence)

output.close()

#Extract genes which have been found with matches

output=open("tmp_genes.fasta", "w")
inread=False
for line in open(sys.argv[5]):
	words=line.strip().split()
	if words[0][0]==">":
		if inread and name in genes_passing:
			print >> output, ">"+name
			print >> output, ''.join(sequence)
		name=words[0][1:]
		sequence=[]
		inread=True
	else:
		sequence.append(words[0].upper())
if inread and name in genes_passing:
	print >> output, ">"+name
	print >> output, ''.join(sequence)

output.close()

print "Running glsearch to find global gene matches"
#run global search to find gene matches in reads
os.system("glsearch -E 1e-10 -T 1 -m 8C tmp_genes.fasta tmp_reads.fasta > tmp_glsearch.out")


#sys.exit()

print "Parsing glsearch output"
#Read the fasta output and save matches for each query to a dictionary

matchdb={}

for line in open("tmp_glsearch.out"):
	words=line.strip().split()
	if words[0]!="#" and len(words)==12:
		query=words[0]
		subject=words[1]
		qstart=int(words[6])
		qend=int(words[7])
		sstart=int(words[8])
		send=int(words[9])
		evalue=float(words[10])
		bitscore=float(words[11])
		if qstart>qend:
			fr="r"
		else:
			fr="f"
		if evalue<1e-50:
			if not subject in matchdb:
				matchdb[subject]={}
			if not query in matchdb[subject]:
				matchdb[subject][query]=[]
			matchdb[subject][query].append([bitscore, sstart, send, fr, evalue])


print "Removing overlapping matches"
#Reduce matches that overlap to the best match decided by bitscore	

querymatchdb={}

for s in matchdb:
	matches=[]
	for q in matchdb[s]:
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
		querymatchdb[um[1]][s].append([um[2], um[3], um[4], um[5]])


print "Creating cluster dictionary"
#Read the cd-hit cluster file to create a cluster dictionary. Also save a dictionary linking the reference gene for each cluster to the cluster number

clusters={}
ref_to_cluster={}
for line in open(sys.argv[2]):
	words=line.strip().split()
	if words[0][0]==">":
		cluster=int(words[1])
		clusters[cluster]=[]
	else:
		gene=words[2][1:-3]
		clusters[cluster].append(gene)
		if words[-1]=="*":
			ref_to_cluster[gene]=cluster

print "Extracting matching reads"
#Extract reads which have been found with matches

reads={}
inread=False
for line in open(sys.argv[3]):
	words=line.strip().split()
	if words[0][0]==">":
		if inread and name in matchdb:
			reads[name]=''.join(sequence)
		name=words[0][1:]
		sequence=[]
		inread=True
	else:
		sequence.append(line)
if inread and name in matchdb:
	reads[name]=''.join(sequence)

print "Extracting matching gene sequences"
#Create sequence database for clusters
cluster_ref_seqs={}
for line in open(sys.argv[4]):
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

print "Aligning matches and creating consensus's"
#align each set of matched sequences and create 50% consensus
consensus_sequences={}
for ref in cluster_ref_seqs:
	musclefile=open("tmp.fasta", "w")
#	print >>  musclefile, ">"+ref
#	print >>  musclefile, cluster_ref_seqs[ref]
	for match in querymatchdb[ref]:
		for x in xrange(len(querymatchdb[ref][match])):
			print >>  musclefile, ">"+match
			if querymatchdb[ref][match][x][2]=="f":
				print >>  musclefile, reads[match][querymatchdb[ref][match][x][0]:querymatchdb[ref][match][x][1]]
			else:
				print >>  musclefile, revcomp(reads[match][querymatchdb[ref][match][x][0]:querymatchdb[ref][match][x][1]])
	
	musclefile.close()
	
	os.system("muscle -in tmp.fasta -out tmp.aln -gapopen -6")
	consensus_sequences[ref]=print_consensus("tmp.aln")
		#sys.exit()


print "Aligning consensus's with clusters"
#Align consensus sequence with cluster
for ref in consensus_sequences:
	output=open("tmp.fasta","w")
	
	print >> output, ">Consensus"
	print >> output, consensus_sequences[ref]
	
	inread=False
	for line in open(sys.argv[5]):
		words=line.strip().split()
		if words[0][0]==">":
			if inread and name in clusters[ref_to_cluster[ref]]:
				print >> output, ">"+name
				print >> output, ''.join(sequence)
			name=words[0][1:]
			sequence=[]
			inread=True
		else:
			sequence.append(words[0].upper())
	if inread and name in clusters[ref_to_cluster[ref]]:
		print >> output, ">"+name
		print >> output, ''.join(sequence)
	output.close()
	
	os.system("muscle -in tmp.fasta -out tmp.aln")
#	os.system("seaview tmp.aln")
	print alignment_stats("tmp.aln")
	
	#read the alignment and print some stats
	
	
	