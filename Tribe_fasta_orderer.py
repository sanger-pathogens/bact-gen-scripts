#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math


if (len (sys.argv)!=2 and len (sys.argv)!=4) or '-h' in sys.argv[1:]:
	print "MUMmer_tiling_to_tab.py <fastafile from tribemcl output>"
	sys.exit()

tribefile=sys.argv[1]

lines=open(tribefile, "rU").read().split('>')[1:]

contigs={}
straincount={}
strains=[]

ref=lines[0].strip().split('_')[1]

for line in lines:
	words=line.strip().split('\n')
	
	if not straincount.has_key(words[0].split('_')[1]):
		straincount[words[0].split('_')[1]]=[0]
		contigs[words[0].split('_')[1]]=[''.join(words[1:])]
		strains.append(words[0].split('_')[1])
	else:
		contigs[words[0].split('_')[1]].append(''.join(words[1:]))
		
output=open('ref.fasta','w')
print >> output, '>'+ref
print >> output, ''.join(contigs[ref])
output.close()

for strain in strains:
	if strain!=ref:
		output=open('contigs.fasta', 'w')
		count=1
		for contig in contigs[strain]:
			print >> output, '>'+strain+'_'+str(count)
			print >> output, contig
			count=count+1
		output.close()
		
		os.system("formatdb -i contigs.fasta -p F")
		os.system("bigger_blast.pl -q normal contigs.fasta ref.fasta")
		os.system("~/scripts/crunch_orderer.py contigs.fasta-qnormalref.fasta.crunch contigs.fasta ref.fasta 90 1000 100")
		os.system("progressiveMauve --output=contig.aln ref.fasta contigs_ordered.fasta")
		
		alines=open('contig.aln', 'rU').read().split('=')
		sequfragments={}
		print "Reading Mauve output..."
		for aline in alines[:-1]:
			posn=-1
			
			lines=aline.split('>')
			
			if int(lines[1].strip().split(':')[0])!=1:
				continue
			
			dirn='+'
			
			for line in lines[1:]:
			
				if len(line)==0:
					continue
				number=int(line.strip().split(':')[0])-1
				
				if number==0:
					posn=int(line.strip().split(':')[1].split('-')[0])
					sequfragments[posn]=[]
					dirn=line.strip().split()[1]
					for i in range(len(subjectnames)+1):
							sequfragments[posn].append('')
				
				if dirn=='+':
					sequfragments[posn][number]=sequfragments[posn][number]+''.join(line.strip().split('\n')[1:])
				elif dirn=='-':
					sequfragments[posn][number]=sequfragments[posn][number]+rev(''.join(line.strip().split('\n')[1:]))
				
		
		
		keys=sequfragments.keys()
		keys.sort()
		
		sequences={ref.split('.')[0]:''}
		seqnames=[ref.split('.')[0]]
		for name in subjectnames:
			sequences[name]=''
			seqnames.append(name)
		
		print "Creating fasta alignment..."
		
		for key in keys:
			fraglen=0
			for x, fragment in enumerate(sequfragments[key]):
				if x==0:
					sequences[ref.split('.')[0]]=sequences[ref.split('.')[0]]+fragment
					fraglen=len(fragment)
				elif len(fragment)==fraglen:
					sequences[subjectnames[x-1]]=sequences[subjectnames[x-1]]+fragment
				else:
					sequences[subjectnames[x-1]]=sequences[subjectnames[x-1]]+'-'*fraglen
		
		
		
		
		
		
		alnout=open(strain+'.aln','w')
		
		for key in seqnames:
			if key!=ref:
				print >> alnout, '>'+key
				print >> alnout, sequences[key]
			
		
		alnout.close()
		
		
		sys.exit()

