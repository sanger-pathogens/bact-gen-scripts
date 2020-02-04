#!/usr/bin/env python

#Compares mummer snp output files.  NOTE: must be from the same reference sequence



#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------

import string, re, math
import os, sys, getopt


def Usage():
	print '\nmauve_iteration_align_multi.py Usage:\n'
	print '\t/nfs/users/nfs_s/mauve_iteration_align_multi.py -r reference.fasta input fasta files'
	print '\nWritten by: Simon R Harris, Wellcome Trust Sanger Institute, UK. 2009\n'



#------------------------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------------------------

def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "hr:", ["help", "ref="])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	ref=''
	contigs=[]

	
	for opt, arg in opts:
		
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-r", "--ref"):
			ref=arg
	
	contigs=args

	if contigs==[]:
		print 'no contig files selected'
		Usage()
		sys.exit()
	elif ref=='':
		print 'no reference file selected'
		Usage()
		sys.exit()

	return ref, contigs

def rev(sequence):
	rev=sequence[::-1]
	
	return rev
	
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

if __name__ == "__main__":
	argv=sys.argv[1:]
	alnfiles=''
	ref, contigs=getOptions(argv)
	refdata=open(ref,'rU').read().split('\n')[1:]
	maxlen=len(''.join(refdata))
	for contig in contigs:
		#newcontig=open(contig.split('.')[0]+'_rearranged.fasta', 'w')
		contigname=contig.split('.')[0]
		
		print "Creating pairwise alignment with Mauve for "+ref.split('.')[0]+" and "+contigname+"..."
		sys.stdout.flush()
		#os.system('/nfs/users/nfs_s/sh16/mauve_2.3.0/progressiveMauve --output='+ref+'.mauvealn1 '+ref+' '+contig+' > tempmauveout.tmp')
		os.system('progressiveMauve --output='+ref+'.mauvealn1 '+ref+' '+contig+' > tempmauveout.tmp')
		alines=open(ref+'.mauvealn1', 'rU').read().split('=')
		sequfragments={}
		
		newstart=-1
		dirn='+'
		
		for aline in alines[:-1]:
			posn=-1
			
			lines=aline.split('>')
			
			dirn='+'
			for line in lines[1:]:
				if len(line)==0:
					continue
				number=int(line.strip().split(':')[0])-1
				if number==0:
					if newstart==-1:
						newstart=int(line.strip().split(':')[1].split('-')[0])
					posn=int(line.strip().split(':')[1].split('-')[0])

					if posn<newstart:
						posn=posn+maxlen
					
					dirn=line.strip().split()[1]
				elif number==1:
					if posn==-1:
						posn=int(line.strip().split(':')[1].split('-')[0])
						posn=posn+maxlen
					#print posn, newstart, dirn, maxlen
					#sequfragments[posn]=sequfragments[posn]+''.join(line.strip().split('\n')[1:])
					if not sequfragments.has_key(posn):
						sequfragments[posn]=''
					if dirn=='+':
						sequfragments[posn]=sequfragments[posn]+''.join(line.strip().split('\n')[1:])
					elif dirn=='-':
						sequfragments[posn]=sequfragments[posn]+revcomp(''.join(line.strip().split('\n')[1:]))
				
		
		
		keys=sequfragments.keys()
		keys.sort()
		
		#print keys
		
		alnseq=''
		for key in keys:
			alnseq=alnseq+sequfragments[key].replace('-','')
			#print sequfragments[key][-20:]
		
		alnfile=open(ref+'_'+contig+'.mauvefna1','w')
		print >> alnfile, ">"+contigname
		print >> alnfile, alnseq
		alnfile.close()
		
		#sys.exit()
		os.system('rm -f *.sml')
			
		alnfiles=alnfiles+' '+ref+'_'+contig+'.mauvefna1'
	#snpfiles.append(prefix+'.snps')

	print "Creating multiple alignment with Mauve..."
	sys.stdout.flush()
	#os.system('/nfs/users/nfs_s/sh16/mauve_2.3.0/progressiveMauve --output='+ref+'.mauvealn --backbone-output='+ref+'.backbone '+ref+' '+alnfiles+' > tempmauveout.tmp')
	os.system('progressiveMauve --output='+ref+'.mauvealn --backbone-output='+ref+'.backbone '+ref+' '+alnfiles+' > tempmauveout.tmp')
	#os.system('progressiveMauve --collinear --output='+ref+'.mauvealn '+ref+' '+alnfiles)
	alines=open(ref+'.mauvealn', 'rU').read().split('=')
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
				for i in range(len(contigs)+1):
						sequfragments[posn].append('')
			
			if dirn=='+':
				sequfragments[posn][number]=sequfragments[posn][number]+''.join(line.strip().split('\n')[1:])
			elif dirn=='-':
				sequfragments[posn][number]=sequfragments[posn][number]+rev(''.join(line.strip().split('\n')[1:]))
			
	
	
	keys=sequfragments.keys()
	keys.sort()
	
	sequences={ref.split('.')[0]:''}
	seqnames=[ref.split('.')[0]]
	for name in contigs:
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
				sequences[contigs[x-1]]=sequences[contigs[x-1]]+fragment
			else:
				sequences[contigs[x-1]]=sequences[contigs[x-1]]+'-'*fraglen
	
	
	
	
	
	
	alnout=open(ref+'_mauve_iteration.aln','w')
	
	for key in seqnames:
		print >> alnout, '>'+key
		print >> alnout, sequences[key]
		
	
	alnout.close()

	os.system('rm -f *.sml tempmauveout.tmp *.mauvefna1 *.mauvealn1*')
	print "Done\n"