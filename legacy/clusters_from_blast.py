#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math


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



if (len (sys.argv)!=3) or '-h' in sys.argv[1:]:
	print "clusters_from_blast.py <blast file> <contigs file>"
	sys.exit()
	
lines=open(sys.argv[2], "rU").read().split('>')[1:]

contigseqs={}

for line in lines:
	words=line.strip().split('\n')
	contigseqs[words[0].split()[0]]=''.join(words[1:])


blastlines=open(sys.argv[1],'rU').readlines()
linesb=[]
for line in blastlines:
	words=line.split()
	if words[0]!=words[1]:
		linesb.append(line)
lasthit=''
lastquery=''
output=open('newblastdb.blastb','w')
for line in linesb:
	words=line.split()
	#if lasthit!=words[1]:
	#	print >> output, lastquery, words[0], lasthit, words[1]
	if lastquery!=words[0] or lasthit==words[1]:
		print >> output, line.strip()
		lastquery=words[0]
		lasthit=words[1]

output.close()

os.system('mv newblastdb.blastb newblastdb.blast')
blastlines=open('newblastdb.blast','rU').readlines()
hits={}

conversiondict={}

for line in blastlines:
	words=line.split()

	if not hits.has_key(words[1]):
		hits[words[1]]={}
	if not hits[words[1]].has_key(words[0]):
		if int(words[8])<int(words[9]):
			hits[words[1]][words[0]]=['f', int(words[8]), int(words[9]), int(words[6]), int(words[7])]
		else:
			hits[words[1]][words[0]]=['r', int(words[9]), int(words[8]), int(words[6]), int(words[7])]
	else:
		if int(words[8])<hits[words[1]][words[0]][1] and hits[words[1]][words[0]][0]=='f' and int(words[8])<int(words[9]):
			hits[words[1]][words[0]][1]=int(words[8])
			hits[words[1]][words[0]][3]=int(words[6])
		elif int(words[9])<hits[words[1]][words[0]][1] and hits[words[1]][words[0]][0]=='r' and int(words[9])<int(words[8]):
			hits[words[1]][words[0]][1]=int(words[9])
			hits[words[1]][words[0]][4]=int(words[7])
		if int(words[8])>hits[words[1]][words[0]][2] and hits[words[1]][words[0]][0]=='r' and int(words[9])<int(words[8]):
			hits[words[1]][words[0]][2]=int(words[8])
			hits[words[1]][words[0]][3]=int(words[6])
		elif int(words[9])>hits[words[1]][words[0]][2] and hits[words[1]][words[0]][0]=='f' and int(words[8])<int(words[9]):
			hits[words[1]][words[0]][2]=int(words[9])
			hits[words[1]][words[0]][4]=int(words[7])


oldhitlen=100000000

iteration=0
while len(hits)<oldhitlen:
	iteration=iteration+1
	print "Iteration:", iteration
	oldhitlen=len(hits)
#	print "Making new fasta file"
#	output=open('newblastdb.fasta','w')
#	for hit in hits.keys():
#		print >> output, '>'+hit
#		print >> output, contigseqs[hit]
#	
#	output.close()
#	
#	print "blasting against new fasta file"
#	os.system('formatdb -i newblastdb.fasta -p F')
#	os.system('blastall -a 4 -m 8 -e 1e-300 -i '+sys.argv[2]+' -d newblastdb.fasta -p blastn -b 2 -v 2 -o newblastdb.blast')
#	lines=open('newblastdb.blast','rU').readlines()
#	linesb=[]
#	for line in lines:
#		words=line.split()
#		if words[0]!=words[1]:
#			linesb.append(line)
#	lasthit=''
#	lastquery=''
#	output=open('newblastdb.blastb','w')
#	for line in linesb:
#		words=line.split()
#		#if lasthit!=words[1]:
#		#	print >> output, lastquery, words[0], lasthit, words[1]
#		if lastquery!=words[0] or lasthit==words[1]:
#			print >> output, line.strip()
#			lastquery=words[0]
#			lasthit=words[1]
#	
#	output.close()
#	
#	os.system('mv newblastdb.blastb newblastdb.blast')
	
	blastlines=open('newblastdb.blast','rU').readlines()
	
	hits={}
	
	conversiondict={}
	
	for line in blastlines:
		words=line.split()
	
		if not hits.has_key(words[1]):
			hits[words[1]]={}
		if not hits[words[1]].has_key(words[0]):
			if int(words[8])<int(words[9]):
				hits[words[1]][words[0]]=['f', int(words[8]), int(words[9]), int(words[6]), int(words[7])]
			else:
				hits[words[1]][words[0]]=['r', int(words[9]), int(words[8]), int(words[6]), int(words[7])]
				#direction, start of subject, end of subject, start of query, end of query
		else:
			#if words[0]=='pool4_2A8_79_length_7291_cov_13.062406':
			#	print  int(words[8]), hits[words[1]][words[0]][2], hits[words[1]][words[0]][0], int(words[9]), int(words[8])
			#	print hits[words[1]][words[0]][2], int(words[8])
			#	print hits[words[1]][words[0]][3], int(words[6])
			#	print hits[words[1]][words[0]]
			if int(words[8])<hits[words[1]][words[0]][1] and hits[words[1]][words[0]][0]=='f' and int(words[8])<int(words[9]):
				hits[words[1]][words[0]][1]=int(words[8])
				hits[words[1]][words[0]][3]=int(words[6])
			elif int(words[9])<hits[words[1]][words[0]][1] and hits[words[1]][words[0]][0]=='r' and int(words[9])<int(words[8]):
				hits[words[1]][words[0]][1]=int(words[9])
				hits[words[1]][words[0]][4]=int(words[7])
			if int(words[8])>hits[words[1]][words[0]][2] and hits[words[1]][words[0]][0]=='r' and int(words[9])<int(words[8]):
				hits[words[1]][words[0]][2]=int(words[8])
				hits[words[1]][words[0]][3]=int(words[6])
			elif int(words[9])>hits[words[1]][words[0]][2] and hits[words[1]][words[0]][0]=='f' and int(words[8])<int(words[9]):
				hits[words[1]][words[0]][2]=int(words[9])
				hits[words[1]][words[0]][4]=int(words[7])
			#if words[0]=='pool4_2A8_79_length_7291_cov_13.062406':
			#	print hits[words[1]][words[0]]
	
	for hit in hits.keys():
		if not hits.has_key(hit):
			continue
		
		for subject in hits[hit]:
			if subject in hits.keys() and hit in hits[subject]:
				if len(contigseqs[subject])<len(contigseqs[hit]) and hits.has_key(subject):
					del hits[subject]
				elif hits.has_key(hit):
					del hits[hit]
	

for count, hit in enumerate(hits.keys()):
	fastafiles=str(count)+'.fasta'
	output=open(fastafiles, 'w')
	outputb=open(str(count)+'.txt', 'w')
	print >> outputb, hit.split('_')[1], hit, 'f'
	print >> output, '>'+hit.split('_')[1]
	print >> output, contigseqs[hit]
	contigs=hits[hit].keys()
	contigs.sort()
	contigsbystrain={}
	for contig in contigs:
		if not contigsbystrain.has_key(contig.split('_')[1]):
			contigsbystrain[contig.split('_')[1]]=[]
		contigsbystrain[contig.split('_')[1]].append([hits[hit][contig][1], hits[hit][contig][2], hits[hit][contig][0], hits[hit][contig][3], hits[hit][contig][4], contig])

	
	#contig order = ref start, ref end, direction, query (contig) start query end, name of contig
	
	strainnames=[hit]
	for strain in contigsbystrain.keys():
		strainbits={}
		contigs=contigsbystrain[strain]
		contigs.sort()
		
		countb=0
		straincontents={}
		straincontents[countb]=''
		for contig in contigs:
			added='n'
			#print contig, len(contigseqs[contig[5]]), len(contigseqs[hit])
			if countb==0:
				straincontents[countb]=straincontents[countb]+', '+contig[5]+' '+contig[2]
				if contig[2]=='f':
					strainbits[countb]=[(contig[0]-contig[3])+1, len(contigseqs[contig[5]])+(contig[0]-contig[3]), contigseqs[contig[5]]]
				else:
					strainbits[countb]=[(contig[0]-contig[3])+1, len(contigseqs[contig[5]])+(contig[0]-contig[3]), revcomp(contigseqs[contig[5]])]
					#strainbits[countb]=[(contig[0]-contig[3])+1, len(contigseqs[contig[5]])+(contig[0]-contig[3]), contigseqs[contig[5]][::-1]]
				countb=countb+1
				straincontents[countb]=''
				continue
			for strainbitkey in strainbits.keys():
				
				if (contig[0]-contig[3])+1 > strainbits[strainbitkey][1] and contig[2]=='f':
					strainbits[strainbitkey][2]=strainbits[strainbitkey][2]+'NNN'+contigseqs[contig[5]]
					strainbits[strainbitkey][1]=(contig[0]-contig[3])+len(contigseqs[contig[5]])
					added='y'
				elif (contig[0]-(len(contigseqs[contig[5]])-contig[4]))+1 > strainbits[strainbitkey][1] and contig[2]=='r':
					#strainbits[strainbitkey][2]=strainbits[strainbitkey][2]+'NNN'+contigseqs[contig[5]][::-1]
					strainbits[strainbitkey][2]=strainbits[strainbitkey][2]+'NNN'+revcomp(contigseqs[contig[5]])
					strainbits[strainbitkey][1]=(contig[1]+contig[3])+1
					
					added='y'
					

				elif len(contigseqs[contig[5]])+(contig[0]-contig[3]) < strainbits[strainbitkey][0] and contig[2]=='f':
					strainbits[strainbitkey][2]=contigseqs[contig[5]]+'NNN'+strainbits[strainbitkey][2]
					strainbits[strainbitkey][0]=(contig[0]-contig[3])+1
					added='y'
				elif contig[1]+contig[3] < strainbits[strainbitkey][0] and contig[2]=='r':
						#print strainbits[strainbitkey]
						#strainbits[strainbitkey][2]=contigseqs[contig[5]][::-1]+'NNN'+strainbits[strainbitkey][2]
					strainbits[strainbitkey][2]=revcomp(contigseqs[contig[5]])+'NNN'+strainbits[strainbitkey][2]
					strainbits[strainbitkey][0]=(contig[1]+contig[3])-len(contigseqs[contig[5]])
					added='y'
				if added=='y':
					break

			if added=='n':
				straincontents[countb]=straincontents[countb]+', '+contig[5]+' '+contig[2]
				if contig[2]=='f':
					strainbits[countb]=[(contig[0]-contig[3])+1, len(contigseqs[contig[5]])+(contig[0]-contig[3]), contigseqs[contig[5]]]
				else:
					strainbits[countb]=[(contig[0]-contig[3])+1, len(contigseqs[contig[5]])+(contig[0]-contig[3]), revcomp(contigseqs[contig[5]])]
					#strainbits[countb]=[(contig[0]-contig[3])+1, len(contigseqs[contig[5]])+(contig[0]-contig[3]), contigseqs[contig[5]][::-1]]
				countb=countb+1
				straincontents[countb]=''
				added='y'
			else:
				straincontents[strainbitkey]=straincontents[strainbitkey]+', '+contig[5]+' '+contig[2]
				
		#if count==62:
		#	print strain, strainbits
		
		for strainbit in strainbits.keys():
			print >> outputb, strain+"_"+str(strainbit+1)+straincontents[strainbit]
			print >> output, '>'+strain+"_"+str(strainbit+1)
			print >> output, strainbits[strainbit][2]
			strainnames.append(strain+"_"+str(strainbit+1))
		
	#print "muscle -in "+fastafiles+" -out "+str(count)+".aln"
	output.close()
	outputb.close()
	#os.system("progressiveMauve --collinear --output="+str(count)+".xmfa "+fastafiles)
	#os.system("rm -f *.sml")
	
	
	xmfalines=open(str(count)+".xmfa", "rU").read().split("=")[0].split(">")[1:]
	
	alnout=open(str(count)+".aln","w")
	for xmfaline in xmfalines:
		header=xmfaline.split("\n")[0]
		print (int(header.split(":")[0].strip()))-1
		print >>alnout, ">"+strainnames[(int(header.split(":")[0].strip()))-1]
		print >>alnout, '\n'.join(xmfaline.split("\n")[1:])
	alnout.close()
	
	#os.system("muscle -in "+fastafiles+" -out "+str(count)+".aln")

			

	
	
	
print len(hits)