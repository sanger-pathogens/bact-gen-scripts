#!/usr/bin/env python

#Compares mummer snp output files.  NOTE: must be from the same reference sequence



#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------

import string, re, math
import os, sys, getopt


def Usage():
	print '\nget_MUMmer_out.py Usage:\n'
	print '\t/nfs/team81/sh16/scripts/get_MUMmer_out.py [options] input files\nOptions:\n\t-a\tAlignment output file\n\t-c\tReference sequence is circular\n\t-d\tDirty mode: Does not remove MUMmer or Mauve files created during the analysis\n\t-e\tEMBL file(s) for reference (comma seperated list)\n\t-i\tInclude indels\n\t-m\tUse previous MUMmer SNP analysis output files\n\t-o\tOutput file name\n\t-t\tArtemis tab output file name\n\t-r\tReference sequence fasta file\n\t-h\tPrint this help'
	print '\nWritten by: Simon R Harris, Wellcome Trust Sanger Institute, UK. 2008\n'



#------------------------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------------------------

def getOptions(arg):

        try:
                opts, args = getopt.getopt(argv, "he:r:o:t:a:imcdCA", ["help", "dirty" "embl=" "ref=", "out=", "align=", "nomummer", "tabfile=", "indels", "circular", "contigs", "abacus"])
        except getopt.GetoptError:
		print "Option Error!", argv
                Usage()
                sys.exit(2)

        ref=''
        outfile=''
        infiles=[]
	i=0
	runmummer='y'
	tabfile=''
	indels='n'
	align=''
	circular='n'
	embl=''
	dirty='n'
	showcontigs='n'
	abacus='n'
	
        for opt, arg in opts:
		
                if opt in ("-h", "--help"):
                        Usage()
                        sys.exit()
                elif opt in ("-r", "--ref"):
                        ref=arg
                elif opt in ("-o", "--out"):
                        outfile=arg
		elif opt in ("-m", "--nomummer"):
                        runmummer='n'
		elif opt in ("-t", "--tabfile"):
                        tabfile=arg
		elif opt in ("-i", "--indels"):
                        indels='y'
		elif opt in ("-c", "--circular"):
                        circular='y'
		elif opt in ("-d", "--dirty"):
                        dirty='y'
		elif opt in ("-a", "--align"):
                        align=arg
		elif opt in ("-e", "--embl"):
                        embl=arg
		elif opt in ("-C", "--contigs"):
                        showcontigs='y'
		elif opt in ("-A", "--abacus"):
                        abacus='y'
	infiles=args
	

        if infiles==[]:
                print 'no input files selected'
                Usage()
                sys.exit()
        if outfile=='':
                outfile="MUMmersnps.out"
        if abacus=='y':
        		runmummer='a'

        return ref, infiles, outfile, runmummer, tabfile, indels, align, circular, embl, dirty, showcontigs


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

def rev(sequence):
	rev=sequence[::-1]
	
	return rev
	
def getCDSseq(embldata, CDSstring, direction, CDScount, feature):
	
	if '(' in CDSstring:
		if CDSstring.split('(')[0]=='complement':
			if direction=='f':
				direction='r'
			else:
				direction='f'
		parts='('.join(CDSstring.split('(')[1:])
		parts=')'.join(parts.split(')')[:-1])
		getCDSseq(embldata, parts, direction, CDScount, feature)
		
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
				
				embldata[donelen+x]=words[1][0]
				if words[1]=='CDS':
					CDSposns[CDScount].append((x, direction))


	return embldata, CDSposns


if __name__ == "__main__":
        argv=sys.argv[1:]
        ref, infiles, outfile, runmummer, tabfile, indels, align, circular, embl, dirty, showcontigs=getOptions(argv)


snps={}
snps['all']={}
snpfiles=[]
subjectnames=[]
alnfiles=''
prefixes=[]

for inputfile in infiles:
	if ref!='':
		prefix=ref.split('.')[0]+'_'+inputfile.split('.')[0]
	else:
		prefix=inputfile.split('.')[0]
	prefixes.append(prefix)
	subjectnames.append(inputfile.split('.')[0])
	if runmummer=='y':
		print "Running MUMmer..."
		os.system('nucmer -l 12 --prefix='+prefix+' '+ref+' '+inputfile)
		os.system('delta-filter -r -q '+prefix+'.delta > '+prefix+'.filter')
#		os.system('delta-filter -r -q -u 50 '+prefix+'.delta > '+prefix+'.filter')
		if indels=='n':
			os.system('show-snps -ClrHTI '+prefix+'.filter > '+prefix+'.snps')
		else:
			os.system('show-snps -ClrHT '+prefix+'.filter > '+prefix+'.snps')
#		if indels=='n':
#			os.system('show-snps -ClrHTI '+prefix+'.delta > '+prefix+'.snps')
#		else:
#			os.system('show-snps -ClrHT '+prefix+'.delta > '+prefix+'.snps')
		if align!='':
			if circular=='y':
				#os.system('show-tiling -c -v 1 -V 1 -p '+prefix+'.mummerfna '+prefix+'.filter > '+prefix+'.tiling')
				print 'show-tiling -i 60 -v 40 -V 1 -u '+prefix+'.unmapped.tab -c -p '+prefix+'.mummerfna '+prefix+'.filter > '+prefix+'.tiling'
				os.system('show-tiling -i 60 -v 40 -V 1 -u '+prefix+'.unmapped.tab -c -p '+prefix+'.mummerfna '+prefix+'.filter > '+prefix+'.tiling')
			else:
				os.system('show-tiling -i 60 -v 40 -V 1 -u '+prefix+'.unmapped.tab '+prefix+'.mummerfna '+prefix+'.filter > '+prefix+'.tiling')
			alnseq=open(prefix+'.mummerfna', 'rU').read().replace('N','')
			alnfile=open(prefix+'.mummerfna','w')
			print >> alnfile, alnseq,
			alnfile.close()
			
			maxlen=int(os.popen('wc -c '+ref).read().split()[0])
			
			
			print "Creating pairwise alignment with Mauve..."
			os.system('progressiveMauve --output='+ref+'.mummeraln1 '+ref+' '+prefix+'.mummerfna')
			alines=open(ref+'.mummeraln1', 'rU').read().split('=')
			sequfragments={}
			print "Reading Mauve output..."
			
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
						dirn=line.strip().split()[1]
					elif number==1:
						if newstart==-1:
							newstart=int(line.strip().split(':')[1].split('-')[0])
						posn=int(line.strip().split(':')[1].split('-')[0])

						if posn<newstart:
							posn=posn+maxlen
						if not sequfragments.has_key(posn):
							sequfragments[posn]=''
						#print posn, newstart, dirn
#						sequfragments[posn]=sequfragments[posn]+''.join(line.strip().split('\n')[1:])
						if dirn=='+':
							sequfragments[posn]=sequfragments[posn]+''.join(line.strip().split('\n')[1:])
						elif dirn=='-':
							sequfragments[posn]=sequfragments[posn]+rev(''.join(line.strip().split('\n')[1:]))
					
			
			
			keys=sequfragments.keys()
			keys.sort()
			
			#print keys
			
			alnseq=''
			for key in keys:
				alnseq=alnseq+sequfragments[key].replace('-','')
				#print sequfragments[key][-20:]
			
			alnfile=open(prefix+'.mummerfna','w')
			print >> alnfile, ">"+inputfile
			print >> alnfile, alnseq
			alnfile.close()
			
			os.system('rm *.sml')
			#sys.exit()
			
			
	elif runmummer=='a':
		os.system('perl /nfs/team81/tdo/Bin/abacas.2.07b.pl -r '+ref+' -q '+inputfile+' -p promer -a 0 -b 1 -N 1 -s 7 -i 50 -v 20')
		os.system('mv '+inputfile+'_'+ref+'.NoNs.fasta '+prefix+'.mummerfna')
			
	alnfiles=alnfiles+' '+prefix+'.mummerfna'
	snpfiles.append(prefix+'.snps')


if align!='':
	print "Creating alignment with Mauve..."
	os.system('progressiveMauve --output='+ref+'.mummeraln '+ref+' '+alnfiles)
	#os.system('progressiveMauve --collinear --output='+ref+'.mummeraln '+ref+' '+alnfiles)
	alines=open(ref+'.mummeraln', 'rU').read().split('=')
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
	
	
	
	
	
	
	alnout=open(align,'w')
	
	for key in seqnames:
		print >> alnout, '>'+key
		print >> alnout, sequences[key]
		
	
	alnout.close()
	
if runmummer=='a':
	sys.exit()

	

refbases={}
bases={}#0=A 1=C 2=G 3=T
nostates=0
converter={}
convertback={}
indelposns={}
contigs={}

print "Reading MUMmer output..."

for inputfile in snpfiles:
	lines=open(inputfile,"rU").readlines()
	snps[inputfile]={}
	contigs[inputfile]={}
	indelposns[inputfile]={}
	for x, line in enumerate(lines):
		words=line.strip().split()
		contigs[inputfile][words[10]+','+words[0]]=words[11]
		if snps['all'].has_key(words[10]+','+words[0]):
			if words[0]==lines[x-1].split()[0]:
				if not indelposns[inputfile].has_key(words[10]+','+words[0]):
					indelposns[inputfile][words[10]+','+words[0]]=1
				indelposns[inputfile][words[10]+','+words[0]]=indelposns[inputfile][words[10]+','+words[0]]+1
			else:
				snps['all'][words[10]+','+words[0]]=snps['all'][words[10]+','+words[0]]+1
		else:
			snps['all'][words[10]+','+words[0]]=1
			refbases[words[10]+','+words[0]]=words[1]
			bases[words[10]+','+words[0]]=[0,0,0,0,0,0,0,0,0,0]
		#print inputfile, words[0], words[1], words[2]
		if not converter.has_key(words[2]):
			converter[words[2]]=nostates
			convertback[nostates]=words[2]
			nostates=nostates+1
		bases[words[10]+','+words[0]][converter[words[2]]] = bases[words[10]+','+words[0]][converter[words[2]]]+1
		snps[inputfile][words[10]+','+words[0]]=1

namesort=snps.keys()
namesort.sort()
subjectnames.sort()

output=open(outfile,'w')
if tabfile!='':
	tabout=open(tabfile,'w')
	print >> tabout, "ID   SNP"
lines=open(ref,'rU').readlines()
chromosomelengths=[]
curlength=0
reflength=0
for line in lines:
	if len(line)>0 and line[0]=='>':
		if curlength!=0:
			chromosomelengths.append(curlength)
			curlength=0
	else:
		curlength=curlength+len(line.strip())
		reflength=reflength+len(line.strip())


print >> output, 'Ref\tSNP_Pos',
if embl!='':
	print >> output, '\tCDS/rRNA/tRNA/Intergenic',
print >> output, '\tRef_base\tSNP_base\tTotal',
for name in subjectnames:
	if name!='all':
		print >> output, '\t'+name,
print >> output, '\tPattern\n',


snpposns={}

for key in snps['all'].keys():

	if snpposns.has_key(key.split(',')[0]):
		snpposns[key.split(',')[0]].append(key.split(',')[1])
	else:
		snpposns[key.split(',')[0]]=[key.split(',')[1]]

keys=snpposns.keys()
keys.sort()


#if embl!='':
#	print "Reading EMBL file(s)..."
#	embldata=[]
#	for i in range(0, reflength+1):
#		embldata.append('I')
#	embls=embl.split(',')
#	donelen=0
#	for em in embls:
#		lines=open(em,'rU').readlines()
#		start=-1
#		end=-1
#		chromosomelen=0
#		for line in lines:
#			words=line.strip().split()
#			#if len(words)>2 and words[1]=='source':
#			#	chromosomelen=int(words[2].split('..')[1])
#			#	for i in range(int(words[2].split('..')[1])):
#			#		embldata.append('I')
#			if len(words)>1 and words[1] in ('CDS', 'rRNA', 'tRNA'):
#			
#				joinlist=words[2].split(',')
#				for join in joinlist:
#					a=int(join.split('..')[0].replace('complement(','').replace('join(','').replace('order(','').replace('<',''))
#					b=int(join.split('..')[1].replace(')','').replace('>',''))
#					if int(a)<int(b):
#						start=a
#						end=b
#					else:
#						start=b
#						end=a
#					for x in range(start,end+1):
#						embldata[donelen+x]=words[1][0]
#		donelen=donelen+chromosomelen



if embl!='':
	print "\nReading EMBL file...",
	sys.stdout.flush()
	CDSposns={}
	CDSpools={}
	embldata=[]
	donelen=0
	CDScount=0		
	for i in range(reflen):
		embldata.append('I')

	lines=open(embl,'rU').readlines()
	
	for x, line in enumerate(lines):
		words=line.strip().split()
		if len(words)>2 and words[1] in ('CDS', 'rRNA', 'tRNA') and '..' in words[2]:
			#print words
			location=words[2].replace('<','').replace('>','')
			extralocationlines=1
			foundpool='n'
			y=x+1
			locustag=''
			while foundpool=='n' and lines[y].split()[1] not in ('CDS', 'rRNA', 'tRNA', 'gene', 'misc_feature', 'sig_peptide', 'repeat_unit', 'repeat_region', 'misc_RNA'):
				if lines[y].split()[1][0]!='/' and y==x+extralocationlines:
					location=location+lines[y].split()[1].strip()
					extralocationlines=extralocationlines+1
				if '/systematic_id' in lines[y].split()[1] or '/locus_tag' in lines[y].split()[1]:
					locustag=lines[y].split()[1].split('"')[1]
					#print locustag
					foundpool='y'
				y=y+1
			
			direction='f'
			feature=words[1]
			if feature=='CDS':
				CDScount=CDScount+1
				CDSposns[CDScount]=[]
				CDSpools[CDScount]=locustag
			embldata, CDSposns=getCDSseq(embldata, location, direction, CDScount, feature)
			

	
	snptoCDSpool={}
	for i in CDSposns.keys():
		for j in CDSposns[i]:
			snptoCDSpool[j[0]]=CDSpools[i]
		
	print "Done"
	sys.stdout.flush()





if embl!='':
	print "\nCalculating dN/dS values...\n"+ref+'...',
	sys.stdout.flush()
	
	comp={'A':'T','T':'A','G':'C','C':'G'}
	
	reflines=open(ref, "rU").readlines()
	refseq=''
	for line in reflines:
		if line.strip()[0]!='>':
			refseq=refseq+line.strip()
	refCDSseq=''
	for key in CDSposns.keys():
		for x in CDSposns[key]:
			if x[1]=='f':
				refCDSseq=refCDSseq+refseq[x[0]-1].upper()
			else:
				refCDSseq=refCDSseq+comp[refseq[x[0]-1].upper()]
	
	print "Done"
	sys.stdout.flush()
	for inputfile in snpfiles:
		print inputfile+'...',
		sys.stdout.flush()
		tempseq=''
		for key in CDSposns.keys():
			for x in CDSposns[key]:
				CDSmapped='n'
				
				if inputfile.mapped[chromosome].has_key(x[0]):
					CDSmapped=inputfile.mapped[chromosome][x[0]]#[0]
						#print inputfile.mapped[chromosome][x[0]][0]
					
				if CDSmapped!='n':
					#print CDSmapped
					if x[1]=='f':
						tempseq=tempseq+CDSmapped.upper()
					else:
						tempseq=tempseq+comp[CDSmapped.upper()]
					
				else:
					tempseq=tempseq+'-'

		
		inputfile.N, inputfile.S, inputfile.dN, inputfile.dS, inputfile.pN, inputfile.pS, inputfile.varS, inputfile.varN, inputfile.z, inputfile.goodlen, inputfile.Nd, inputfile.Ns=dnbyds(refCDSseq, tempseq)
		
		inputfile.CDSseq=tempseq
		
		if inputfile.dS!=0:
			print "%.2f" % (inputfile.dN/inputfile.dS)
		sys.stdout.flush()





print "Writing output files..."

donelength=0
for y, key in enumerate(keys):
	if y>0:
		donelength=donelength+chromosomelengths[y-1]
	snpsort=map(lambda value:int(value), snpposns[key])
	snpsort.sort()
	inindel='n'
	
	for m, j in enumerate(snpsort):
		i=key+','+str(j)
		if refbases[i]=='.':
			continue
		
		
			
		
		print >> output, '\t'.join(i.split(',')),
		if embl!='':
			print >> output, '\t'+embldata[donelength+j],
		
		print >> output, '\t'+refbases[i],
		snpbase=''
		snpbasecount=0
		for k, base in enumerate(bases[i]):
			if base>snpbasecount:
				snpbasecount=base
				snpbase=convertback[k]
		
		
		if inindel=='y' and (snpbase!='.' or snpsort[m-1]!=j-1):
			inindel='n'
		
		if tabfile!='':
			if snpbase=='.' and inindel=='n':
				inindel='a'
				countadded=1
				indelend=j
				indelstart=j
				newsnpbase=snpbase
				
				
				while newsnpbase=='.' and m+countadded<len(snpsort) and snpsort[m+countadded]==j+countadded:
					snpbasecount=0
					l=key+','+str(snpsort[m+countadded])
					for k, base in enumerate(bases[l]):
						if base>snpbasecount:
							snpbasecount=base
							newsnpbase=convertback[k]
					if newsnpbase=='.':
						indelend=snpsort[m+countadded]
					countadded=countadded+1
				
				if indelstart==indelend:
					print >> tabout, 'FT   Indel           '+str(donelength+indelstart)
				else:
					print >> tabout, 'FT   Indel           '+str(donelength+indelstart)+'..'+str(donelength+indelend)
				print >> tabout, 'FT                   /strains="', 
			
			elif inindel=='n':	
				print >> tabout, 'FT   SNP             '+str(donelength+j)
				print >> tabout, 'FT                   /colour="2"'
				print >> tabout, 'FT                   /refAllele="'+refbases[i]+'"'
				print >> tabout, 'FT                   /SNPAllele="'+snpbase+'"'
				print >> tabout, 'FT                   /SNPstrains="',
				
	
		print >> output, '\t'+snpbase+'\t'+str(snps['all'][i]),	
	
		
	
		pattern=''
		for x, name in enumerate(namesort):
			if name !='all':
				if snps[name].has_key(i):
					if tabfile!='' and inindel!='y':
						print >> tabout, subjectnames[x],
					print >> output, '\t'+str(snps[name][i]),
					pattern=pattern+'1'
				else:
					print >> output, '\t0',
					pattern=pattern+'0'
				if showcontigs=='y':
					print >> output, '\t'+contigs[name][i],
		
		if tabfile!='' and inindel!='y':
			print >> tabout, '"'
		
		if inindel=='a':
			inindel='y'
		
		print >> output, '\t'+pattern+'\n',
				
				
output.close()
if tabfile!='':
	tabout.close()
	
if dirty=='n':
	print "Cleaning up files..."
	for p in prefixes:
		cleanline='rm '+p+'.cluster '+p+'.delta '+p+'.filter '+p+'.snps '
		if align!='':
			cleanline=cleanline+p+'.tiling '+p+'.mummerfna*'
		os.system(cleanline)
	if align!='':
		cleanline=cleanline+ref+'.mummeraln* '+ref+'.sml'
		os.system(cleanline)
	
os.system('rm *.sml')
print "Done\n"		
				
