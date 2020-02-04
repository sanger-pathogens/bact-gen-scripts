#!/usr/bin/env python
#/usr/bin/env python
import string
import os, sys
import pysam
from optparse import OptionParser
from numpy import mean, std, max
sys.path.extend(map(os.path.abspath, ['/nfs/pathogen/sh16_scripts/modules/']))
from Si_general import *
from Si_SeqIO import *

#################################
# Simple Error Printing Funtion #
#################################

def DoError(ErrorString):
	print "!!!Error:", ErrorString,"!!!"
	sys.exit()

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

	parser.add_option("-b", "--bam", action="store", dest="bam", help="input file in sam or bam format (must have suffix .sam or .bam respectively)", default="", metavar="FILE")
	parser.add_option("-o", "--output", action="store", dest="output", help="output prefix for file name(s)", default="", metavar="FILE")
	parser.add_option("-r", "--reference", action="store", dest="ref", help="input reference fasta file", default="", metavar="FILE")
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):

	if options.bam=='':
		DoError('No bam/sam input file selected')
	elif options.bam.split('.')[-1] not in ['sam', 'bam']:
		DoError('Input file '+options.bam+' must have the suffix .sam or .bam')
	elif not os.path.isfile(options.bam):
		DoError('Cannot find file '+options.bam)
	if options.output=="":
		DoError('No output file name (-o) selected')
	elif not os.path.isfile(options.ref):
		DoError('Cannot find file '+options.ref)

	
	return



########
# Main #
########

if __name__ == "__main__":
	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)


	# read the reference fasta file

	refrecords=read_seq_file(options.ref)

	
	#open the bam/sam file
	
	if options.bam.split(".")[-1]=="bam":
#		print "Reading bam file"
		format="bam"
		try:
			samfile = pysam.Samfile( options.bam, "rb" )
		except StandardError:
			DoError('Failed to open '+options.bam+'. Is it in bam format?')
	elif options.bam.split(".")[-1]=="sam":
#		print "Reading sam file"
		format="sam"
		try:
			samfile = pysam.Samfile( options.bam, "r" )
		except StandardError:
			DoError('Failed to open '+options.bam+'. Is it in sam format?')
	
	#get reference names and lengths from the sam file header
	refs=samfile.references
	lengths=samfile.lengths

	reflengths={}
	for x, ref in enumerate(refs):
		reflengths[ref]=lengths[x]

	refseqs={}
	for record in refrecords:
		if record.name not in reflengths or len(str(record.seq))!=reflengths[record.name]:
			DoError("Sequence "+record.name+" not in bam file header")
		refseqs[record.name]=str(record.seq)


	frefreps={}
	rrefreps={}
	for record in refrecords:
		frefreps[record.name]={}
		rrefreps[record.name]={}
		for replen in xrange(1,3):
			basenum=0
			while basenum<reflengths[record.name]-replen:
				
				
				baseseq=refseqs[record.name][basenum:basenum+replen]
					
				currpos=basenum+replen

				while reflengths[record.name]>currpos+replen and refseqs[record.name][currpos:currpos+replen]==baseseq:
					currpos+=replen

				numreps=((currpos-basenum)/replen)-1
				totalreplen=(numreps+1)*replen
				if numreps>0:
					totalreplen=(numreps+1)*replen
					#print record.name, baseseq, basenum, currpos, numreps, (numreps+1)*baseseq
					x=0
					while x<totalreplen:
						if not basenum+x in frefreps[record.name] or frefreps[record.name][basenum+x]<totalreplen-x:
							frefreps[record.name][basenum+x]=totalreplen-x
						if not basenum+x in rrefreps[record.name] or rrefreps[record.name][basenum+x]<x:
							rrefreps[record.name][basenum+x]=x
						x+=1
					basenum+=numreps*replen
				basenum+=1

	

	
	#open the output file

	if format=="sam":
		output=pysam.Samfile(options.output+".sam", mode='wh', referencenames=refs, referencelengths=lengths)
	elif format=="bam":
		output=pysam.Samfile(options.output+".bam", mode='wb', referencenames=refs, referencelengths=lengths)
	
	
	
	#Iterate through input file and fix appropriate lines
	samfile.reset()
	count=0
	rpos=0
	rnum=1
	for read in samfile:
		#if read.pos<820000 or read.pos>830000:
			#continue
		
		if not read.is_unmapped:
			start=read.pos
			readpos=0
			refpos=start
			newcigarpos=0
			contig=samfile.getrname(read.rname)
			
			#trim reads that run into repeats
			mapposs=[]
			maptoread=[]
			maptoref=[]
			for cig in read.cigar:
				
				if cig[0]==0:
					for x in xrange(cig[1]):
						mapposs.append(refpos+x)
						maptoread.append(readpos+x+1)
						maptoref.append(refpos+x+1)
					readpos+=cig[1]
					refpos+=cig[1]
				elif cig[0]==4:
					readpos+=cig[1]
				elif cig[0]==1:
					readpos+=cig[1]
				elif cig[0]==2:
					refpos+=cig[1]

			
			ffound=False
			rfound=False
			ftrim=float("-Inf")
			rtrim=float("Inf")
			
			for x, base in enumerate(mapposs):
				
				if not ffound and base in frefreps[contig]:
					#sys.exit()
					if (len(mapposs)-x)<(frefreps[contig][base]+1):
						#print "ftrim", frefreps[contig][base], base, maptoread[x]
						ftrim=maptoread[x]-1
						ffound=True
				if base in rrefreps[contig]:
					#print "r", base, rrefreps[contig][base]
					if (x+1)<(rrefreps[contig][base]+1):
						#print "rtrim", rrefreps[contig][base], base, maptoread[x]
						rtrim=maptoread[x]
						newstart=maptoref[x]
						rfound=True

			readpos=0
			newcigar=[]
			if ffound or rfound:
					#sys.exit()
				fcut=False
				if not rfound:
					rcut=True
				else:
					rcut=False
					read.pos=newstart

					
				for cig in read.cigar:
					oldreadpos=readpos
					
					if cig[0]==0:
						readpos+=cig[1]
					elif cig[0]==4:
						readpos+=cig[1]
					elif cig[0]==1:
						readpos+=cig[1]

					if (ffound and ftrim>=oldreadpos and ftrim<=readpos) and (rfound and rtrim<=readpos):
						ffound=False
						fcut=True
						rfound=False
						rcut=True
						newcigar.append((4,rtrim))
						diff=ftrim-rtrim
						if diff>0:
							newcigar.append((cig[0],diff))
						newcigar.append((4,len(read.seq)-ftrim))
						
					elif ffound and ftrim>=oldreadpos and ftrim<=readpos:
						ffound=False
						fcut=True
						diff=ftrim-oldreadpos
						if diff>0:
							newcigar.append((cig[0],diff))
						newcigar.append((4,len(read.seq)-ftrim))
						
					elif rfound and rtrim<=readpos:
						rfound=False
						rcut=True
						newcigar.append((4,rtrim))
						diff=readpos-rtrim
						if diff>0:
							newcigar.append((cig[0],diff))
							
					elif not fcut and rcut:
						newcigar.append(cig)
						
				
				
							
						
			else:
				newcigar=read.cigar

			
			
					
			read.cigar=newcigar
			
			start=read.pos
			readpos=0
			refpos=start
			newcigar=read.cigar
			newcigarpos=0
			for cig in read.cigar:
				
				if cig[0]==0:
					readpos+=cig[1]
					refpos+=cig[1]
				elif cig[0]==4:
					readpos+=cig[1]
				elif cig[0]==1:
					indellength=cig[1]
					indelpos=refpos
					
					if read.is_reverse:
						revcompseq=revcomp(read.seq)
						indel=revcompseq[readpos:readpos+indellength]
					else:
						indel=read.seq[readpos:readpos+indellength]

					
					shortestrepeat=indellength
					for x in xrange(indellength,0,-1):
						if indellength%x==0:
							repeat=True
							repeatseq=refseqs[contig][refpos:refpos+x]
							for y in xrange(0,indellength,x):
								if refseqs[contig][refpos+y:refpos+y+x]!=repeatseq:
									repeat=False
							if repeat:
								shortestrepeat=x

					#print indellength, indel, readpos, indelpos, shortestrepeat, read.cigar			       
					
					while refseqs[contig][indelpos-shortestrepeat:indelpos+(indellength-shortestrepeat)]==indel:
						#print read.pos, indelpos, indel, refseqs[samfile.getrname(read.rname)][indelpos-indellength:indelpos]
						indelpos-=shortestrepeat
						count+=1
						#print count
					
					
					diff=refpos-indelpos
					if refpos!=indelpos:
						print indel, refpos, indelpos
						sys.exit()
					readpos+=cig[1]
					
				elif cig[0]==2:
					#print read.pos, newcigar
					indellength=cig[1]
					indelpos=refpos
					indel=refseqs[contig][refpos:refpos+indellength]
					shortestrepeat=indellength
					for x in xrange(indellength,0,-1):
						if indellength%x==0:
							repeat=True
							repeatseq=refseqs[contig][refpos:refpos+x]
							for y in xrange(0,indellength,x):
								if refseqs[contig][refpos+y:refpos+y+x]!=repeatseq:
									repeat=False
							if repeat:
								shortestrepeat=x

								       
					
					while refseqs[contig][indelpos-shortestrepeat:indelpos+(indellength-shortestrepeat)]==indel:
						#print read.pos, indelpos, indel, refseqs[samfile.getrname(read.rname)][indelpos-indellength:indelpos]
						indelpos-=shortestrepeat
						count+=1
						#print count
					
					
					diff=refpos-indelpos
					

					if diff!=0:
						if newcigarpos==0:
							if len(newcigar)>newcigarpos+1:
								if newcigar[newcigarpos+1][0]==0:
									newcigar[newcigarpos+1]=(newcigar[newcigarpos+1][0],newcigar[newcigarpos+1][1]+diff)
								else:
									newcigar=newcigar[:newcigarpos]+[(0,diff)]+newcigar[newcigarpos:]
									newcigarpos+=1
						elif len(newcigar)>newcigarpos+1:
							if newcigar[newcigarpos-1][0]!=0:

								print "weird", newcigar[newcigarpos-1]
								
							if newcigar[newcigarpos-1][1]<=diff:
								newcigar.pop(newcigarpos-1)
								newcigarpos-=1
							else:
								newcigar[newcigarpos-1]=(newcigar[newcigarpos-1][0],newcigar[newcigarpos-1][1]-diff)
							
								
							if newcigar[newcigarpos+1][0]==0:
								newcigar[newcigarpos+1]=(newcigar[newcigarpos+1][0],newcigar[newcigarpos+1][1]+diff)
							else:
								newcigar=newcigar[:newcigarpos]+[(0,diff)]+newcigar[newcigarpos:]
								newcigarpos+=1
						else:
							if newcigar[newcigarpos-1][0]!=0:

								print "weird", newcigar[newcigarpos-1]
								
							if newcigar[newcigarpos-1][1]<=diff:
								newcigar.pop(newcigarpos-1)
								newcigarpos-=1
							else:
								newcigar[newcigarpos-1]=(newcigar[newcigarpos-1][0],newcigar[newcigarpos-1][1]-diff)
							
								
							newcigar=newcigar+[(0,diff)]
							newcigarpos+=1
					#refpos=(refpos-indelpos)+indellength
					#readpos=(readpos-indelpos)+indellength
					refpos=indelpos+indellength
					readpos=(readpos-diff)+indellength
						
					if diff!=0:
						ciglen=0
						newciglen=0
						for cig in read.cigar:
							if cig[0] in [0,1,3,4,5,6]:
								ciglen+=cig[1]
						for cig in newcigar:
							if cig[0] in [0,1,3,4,5,6]:
								newciglen+=cig[1]
						if ciglen!=newciglen:
							print read.pos, read.cigar, newcigar, ciglen, newciglen
							newcigar=read.cigar #NEED TO FIX THESE READS!!!!!
					#read.cigar=newcigar
					
					
					
				else:
					print cig
				newcigarpos+=1
				
			read.cigar=newcigar
				
			output.write(read)

			rpos=float(read.pos)/rnum
			
			if rpos>=10000:
				print 10000*rnum
				rnum+=1
			if read.pos>100000:
				break
	
	output.close()
	if format=="bam":
		
		os.system("samtools index "+options.output+".bam")

		    	
			
