#!/usr/bin/env python
import string, re, gzip
import os, sys, getopt, random, math


###############################################################
# Recursive algorithm to reconstruct CDSs from EMBL CDS lines #
###############################################################

def getCDSseq(embldata, CDSstring, direction, CDScount, feature, CDSposns, CorP):
	
	if '(' in CDSstring:
		if CDSstring.split('(')[0]=='complement':
			if direction=='f':
				direction='r'
			else:
				direction='f'
		parts='('.join(CDSstring.split('(')[1:])
		parts=')'.join(parts.split(')')[:-1])
		getCDSseq(embldata, parts, direction, CDScount, feature, CDSposns, CorP)
		
	else:
		joinlist=CDSstring.split(',')
		for join in joinlist:
			a=int(join.split('..')[0])
			b=int(join.split('..')[1])
			if direction=='f':
				start=a
				end=b+1
				y=1
			else:
				start=b
				end=a-1
				y=-1
			for x in range(start,end,y):
				
				embldata[donelen+x]=CorP
				if words[1]=='CDS':
					CDSposns[CDScount].append((x, direction))


	return embldata, CDSposns



def Usage():
	print 'output_covered_genes.py Usage:'
	print 'output_covered_genes.py'
	print 'Written by: Simon R Harris, Wellcome Trust Sanger Institute, UK. 2009'


#------------------------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------------------------


def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "he:r:d:o:", ["help", "embl=", "ref=", "out=", "depth="])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	outfile=''
	plots=[]
	embl=''
	depth=1
	ref=''

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-e", "--embl"):
			embl=arg
		elif opt in ("-o", "--out"):
			outfile=arg
		elif opt in ("-d", "--depth"):
			depth=int(arg)
		elif opt in ("-r", "--ref"):
			ref=arg

	plots=args

	if plots==[] or embl=='' or ref=='':
		print 'Error: Missing input file!'
		Usage()
		sys.exit()
	if outfile=='':
		print 'Error: Missing output file!'
		Usage()
		sys.exit()

	return outfile, plots, embl, depth, ref



if __name__ == "__main__":
        argv=sys.argv[1:]
        outfile, plots, embl, depth, ref=getOptions(argv)

refin=open(ref, "rU").read().split('\n')[1:]
reflen=len(''.join(refin))


print "\nReading EMBL file...",
sys.stdout.flush()
CDSposns={}
CDSnames={}
CDSproducts={}
CDSgenes={}
embldata=[]
donelen=0
CDScount=0
pseudoposns={}
pseudonames={}
pseudocount=0
for i in range(reflen+1):
	embldata.append('I')

lines=open(embl,'rU').readlines()

for x, line in enumerate(lines):
	words=line.strip().split()
	if len(words)>2 and words[1] in ('CDS', 'rRNA', 'tRNA') and '..' in words[2]:
		#print words
		location=words[2].replace('<','').replace('>','')
		extralocationlines=1
		foundpool='n'
		foundpool2='n'
		
		foundpool3='n'
		y=x+1
		locustag=''
		pseudo='n'
		product=''
		gene=''
		while (foundpool=='n' or foundpool2=='n' or foundpool3=='n') and len(lines[y].split())>1 and lines[y].split()[1] not in ('CDS', 'rRNA', 'tRNA', 'gene', 'misc_feature', 'sig_peptide', 'repeat_unit', 'repeat_region', 'misc_RNA'):
			if lines[y].split()[1][0]!='/' and y==x+extralocationlines:
				location=location+lines[y].split()[1].strip()
				extralocationlines=extralocationlines+1
			if '/colour=11' in lines[y].split()[1]:
				pseudo='y'
			if '/systematic_id' in lines[y].split()[1] or '/locus_tag' in lines[y].split()[1]:
				locustag=lines[y].split()[1].split('"')[1]
				#print locustag
				foundpool='y'
			if '/product' in lines[y].split()[1]:
				product=' '.join(lines[y].split()[1:]).split('"')[1]
				#print locustag
				foundpool2='y'
			if '/gene' in lines[y].split()[1]:
				gene=' '.join(lines[y].split()[1:]).split('"')[1]
				#print locustag
				foundpool3='y'
			y=y+1
		
		direction='f'
		feature=words[1]
		if feature=='CDS' and pseudo=='n':
			CDScount=CDScount+1
			CDSposns[CDScount]=[]
			CDSnames[CDScount]=locustag
			CDSproducts[CDScount]=product
			CDSgenes[CDScount]=gene
		elif feature=='CDS':
			pseudocount=pseudocount+1
			pseudoposns[pseudocount]=[]
			pseudonames[pseudocount]=locustag
			embldata, pseudoposns=getCDSseq(embldata, location, direction, pseudocount, feature, pseudoposns, 'P')
		if pseudo=='n'and feature in ['CDS', 'rRNA', 'tRNA']:
			embldata, CDSposns=getCDSseq(embldata, location, direction, CDScount, feature, CDSposns, feature[0])
		
	
print "Done"
sys.stdout.flush()

all={}
allinfo={}
strainnumber=0
for plot in plots:
	coverage=gzip.open(plot,'r').read().split('\n')[:-1]
	
	present={}
	
	for CDS in range(1,CDScount+1):
	
		totalcov=0.0
		totalwithcov=0.0
	
		for x in CDSposns[CDS]:
			if x[0]<len(coverage) and int(coverage[x[0]])>=depth:
				totalwithcov=totalwithcov+1
				totalcov=totalcov+int(coverage[x[0]])
		
		if (totalwithcov/len(CDSposns[CDS]))*100>95:
			if CDSgenes[CDS]!='' and not present.has_key(CDSgenes[CDS]):
				present[CDSgenes[CDS]]=[CDSgenes[CDS], CDSproducts[CDS], str(totalcov/len(CDSposns[CDS])), str((totalwithcov/len(CDSposns[CDS]))*100), str(1)]
			elif CDSgenes[CDS]!='':
				present[CDSgenes[CDS]][4]=str(int(present[CDSgenes[CDS]][4])+1)
			elif not present.has_key(CDSproducts[CDS]):
				present[CDSproducts[CDS]]=[CDSgenes[CDS], CDSproducts[CDS], str(totalcov/len(CDSposns[CDS])), str((totalwithcov/len(CDSposns[CDS]))*100), str(1)]
			else:
				present[CDSproducts[CDS]][4]=str(int(present[CDSproducts[CDS]][4])+1)
		
		if CDSgenes[CDS]!='':
			if not all.has_key(CDSgenes[CDS]):
				all[CDSgenes[CDS]]=[]
				allinfo[CDSgenes[CDS]]=[CDSgenes[CDS], CDSproducts[CDS]]
				if (totalwithcov/len(CDSposns[CDS]))*100>95:
					all[CDSgenes[CDS]].append(str(totalcov/len(CDSposns[CDS])))
				else:
					all[CDSgenes[CDS]].append('0.0')
			if len(all[CDSgenes[CDS]])==(strainnumber+1):
				if(totalwithcov/len(CDSposns[CDS]))*100>95 and totalcov/len(CDSposns[CDS])>float(all[CDSgenes[CDS]][strainnumber]):
					all[CDSgenes[CDS]][strainnumber]=str(totalcov/len(CDSposns[CDS]))
			else:
				if (totalwithcov/len(CDSposns[CDS]))*100>95:
					all[CDSgenes[CDS]].append(str(totalcov/len(CDSposns[CDS])))
				else:
					all[CDSgenes[CDS]].append('0.0')
		
		elif CDSproducts[CDS]!='':
			if not all.has_key(CDSproducts[CDS]):
				all[CDSproducts[CDS]]=[]
				allinfo[CDSproducts[CDS]]=[CDSgenes[CDS], CDSproducts[CDS]]
				if (totalwithcov/len(CDSposns[CDS]))*100>95:
					all[CDSproducts[CDS]].append(str(totalcov/len(CDSposns[CDS])))
				else:
					all[CDSproducts[CDS]].append('0.0')
			if len(all[CDSproducts[CDS]])==(strainnumber+1):
				if (totalwithcov/len(CDSposns[CDS]))*100>95 and totalcov/len(CDSposns[CDS])>float(all[CDSproducts[CDS]][strainnumber]):
					all[CDSproducts[CDS]][strainnumber]=str(totalcov/len(CDSposns[CDS]))
			else:
				if (totalwithcov/len(CDSposns[CDS]))*100>95:
					all[CDSproducts[CDS]].append(str(totalcov/len(CDSposns[CDS])))
				else:
					all[CDSproducts[CDS]].append('0.0')

			
	strainnumber=strainnumber+1

output=open(outfile, 'w')

print >> output, '\t'.join(["gene","product"]+plots)
for key in all.keys():
	print >> output, '\t'.join(allinfo[key]+all[key])

output.close()












	
