#!/usr/bin/env python

#Compares mummer snp output files.  NOTE: must be from the same reference sequence



#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------

import string, re
import os, sys, getopt


def Usage():
	print 'compare_mummer_out.py Usage:'
	print 'compare_mummer_out.py -r=[reference sequence] [input alignment(s)] > [output file name] {-h}'
	print 'Copyright Simon R Harris, Wellcome Trust Sanger Institute, UK. 2008'

	

#------------------------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------------------------


def getOptions(arg):

        try:
                opts, args = getopt.getopt(argv, "he:r:o:t:a:imcd", ["help", "dirty" "embl=" "ref=", "out=", "align=", "nomummer", "tabfile=", "indels", "circular"])
        except getopt.GetoptError:
                print "Option Error!", argv
                Usage()
                sys.exit(2)

        ref=''
        outfile=''
        inputdirs=[]
        i=0
        runmaq='y'
        tabfile=''
        indels='n'
        align=''
        circular='n'
        embl=''
        dirty='n'

        for opt, arg in opts:

                if opt in ("-h", "--help"):
                        Usage()
                        sys.exit()
                elif opt in ("-r", "--ref"):
                        ref=arg
                elif opt in ("-o", "--out"):
                        outfile=arg
                elif opt in ("-m", "--nomaq"):
                        runmaq='n'
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

        inputdirs=args


        if inputdirs==[]:
                print 'no input files selected'
                Usage()
                sys.exit()
        if outfile=='':
                outfile=ref.split('.')[0]+"Maqsnps.out"

        return ref, inputdirs, outfile, tabfile, align, embl, dirty, runmaq





class SNPanalysis:
	def __init__(self, fastq='', directory='', SNPs={}):
		self.fastq=fastq
		self.directory=directory
		self.SNPs=SNPs





if __name__ == "__main__":
        argv=sys.argv[1:]
        ref, inputdirs, outfile, tabfile, align, embl, dirty, runmaq=getOptions(argv)


snps={}

refbases={}
bases={}#0=A 1=C 2=G 3=T
nostates=0
converter={}
convertback={}
snpbases={}

snpfiles=[]
count=0
namesort=[]
for inputfile in inputdirs:
	if os.path.isfile(inputfile+'/ssaha.mapped.snp'):
		snpfiles.append(SNPanalysis())
		snpfiles[count].directory=inputfile
		namesort.append(inputfile)
		count=count+1

allsnps={}

for inputfile in snpfiles:
	lines=open(inputfile.directory+'/ssaha.mapped.snp',"rU").readlines()
	snps[inputfile.directory]={}
	for line in lines[6:-1]:
		words=line.strip().split()
		if len(words)>2 and int(words[2])>8:
			snps[inputfile.directory][words[1]+','+words[3]]=1
			if allsnps.has_key(words[1]+','+words[3]):
				allsnps[words[1]+','+words[3]]=allsnps[words[1]+','+words[3]]+1
			else:
				allsnps[words[1]+','+words[3]]=1
				refbases[words[1]+','+words[3]]=words[5]
				bases[words[1]+','+words[3]]=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
			#print inputfile, words[1], words[2], words[3]
			if not converter.has_key(words[6]):
				converter[words[6]]=nostates
				convertback[nostates]=words[6]
				nostates=nostates+1

			bases[words[1]+','+words[3]][converter[words[6]]] = bases[words[1]+','+words[3]][converter[words[6]]]+1


namesort.sort()

output=open(outfile,'w')

if embl!='':
        print >> output, '\tCDS/rRNA/tRNA/Intergenic',
print >> output, '\tRef_base\tSNP_base\tTotal',

for name in namesort:
	print >> output, '\t'+name,
print >> output, '\tPattern\n',


snpposns={}


for key in allsnps.keys():

	if snpposns.has_key(key.split(',')[0]):
		snpposns[key.split(',')[0]].append(key.split(',')[1])
	else:
		snpposns[key.split(',')[0]]=[key.split(',')[1]]



#Make pilupfile for ALL SNPs for each alignment - so we can assess the bases
pilein=open('pile_up_input.txt', 'w')
keys=snpposns.keys()
keys.sort()

for i in keys:
	for j in snpposns[i]:
		print >> pilein, i+'\t'+j
pilein.close()


avcoverage={}
covcount={}



lenstring=os.popen('grep "reference length:" '+snpfiles[0].directory+'/assemble.log').read()
reflen=int(lenstring.strip().split(":")[1])

if embl!='':
        print "Reading EMBL file(s)..."
        embldata=[]
        embls=embl.split(',')
        donelen=0
	for i in range(reflen):
        	embldata.append('I')
        for em in embls:
                lines=open(em,'rU').readlines()
                start=-1
                end=-1
                chromosomelen=0
                for line in lines:
                        words=line.strip().split()
                        #if len(words)>2 and words[1]=='source':
                         #       chromosomelen=int(words[2].split('..')[1])
                        
                        if len(words)>1 and words[1] in ('CDS', 'rRNA', 'tRNA'):

                                joinlist=words[2].split(',')
                                for join in joinlist:
                                        a=int(join.split('..')[0].replace('complement(','').replace('join(','').replace('order(',''))
                                        b=int(join.split('..')[1].replace(')',''))
                                        if int(a)<int(b):
                                                start=a
                                                end=b
                                        else:
                                                start=b
                                                end=a
                                        for x in range(start,end+1):
                                                embldata[donelen+x]=words[1][0]
                donelen=donelen+chromosomelen







#print summary file and alignment file

if align!='':
	aln={}
	alnfile=open(align, 'w')
if tabfile!='':
        tabout=open(tabfile,'w')
snpgraph=open('snpgraph.txt','w')
covgraph=open('covgraph.txt','w')
covcountgraph=open('covcountgraph.txt','w')

minquality=10
minqualityperdepth=4

sequence={}
goodhits={}

refsequence=''

#lenstring=os.popen('grep "reference length:" '+snpfiles[0].directory+'/assemble.log').read()
#reflen=int(lenstring.strip().split(":")[1])


print "Writing output files..."

if tabfile!='':
	print >> tabout, 'ID   SNP'



for key in keys:

	snpsort=map(lambda value:int(value), snpposns[key])
	snpsort.sort()
	
	for j in snpsort:
		
		outstring=''
		tabstring=''
		
		i=key+','+str(j)
		outstring=outstring+'\t'.join(i.split(','))
		if embl!='':
                        outstring=outstring+'\t'+embldata[j]
		outstring=outstring+'\t'+refbases[i]

		refsequence=refsequence+refbases[i]
		snpbase=''
		snpbasecount=0
		for k, base in enumerate(bases[i]):
			if base>snpbasecount:
				snpbasecount=base
				snpbase=convertback[k]
	
		outstring=outstring+'\t'+snpbase+'\t'+str(allsnps[i])
		
		if tabfile!='':
                	tabstring=tabstring+'FT   SNP             '+str(j)+'\n'
			tabstring=tabstring+'FT                   /refAllele="'+refbases[i]+'"\n'
			tabstring=tabstring+'FT                   /SNPAllele="'+snpbase+'"\n'
			tabstring=tabstring+'FT                   /SNPstrains="'

		

		pattern=''
		numsnpbases=1
		snpbases=[]
		for name in snpfiles:
#			if name.SNPs[j][2]>8:# and name.SNPs[j][2]>(name.SNPs[j][3]*4):
#			#print >> output, '\t'+str(name.SNPs[j][2])
				
	#		if name.SNPs[j][0]!=name.SNPs[j][1]:
	#			if tabfile!='':
	#				tabstring=tabstring+name.directory+'='+name.SNPs[j][1]+' '
	#			pattern=pattern+'1'
	#			outstring=outstring+'\t'+str(name.SNPs[j][2])
	#			if name.SNPs[j][1] not in snpbases:
	#				snpbases.append(name.SNPs[j][1])
	#				numsnpbases=numsnpbases+1
	#		else:
	#			pattern=pattern+'0'
	#			outstring=outstring+'\t0'
	#		sequence[name.directory]=sequence[name.directory]+name.SNPs[j][1]
	#		goodhits[name.directory]=goodhits[name.directory]+1
	#		else:
	#			#print >> output, '\t'+str(name.SNPs[j][2])
	#			pattern=pattern+'-'
	#			sequence[name.directory]=sequence[name.directory]+'-'
	#			outstring=outstring+'\t-'
		
			if snps[name.directory].has_key(i):
				outstring=outstring+'\t'+str(snps[name.directory][i])
				pattern=pattern+'1'
			else:
				outstring=outstring+'\t0'
				pattern=pattern+'0'
			#if snpbases[name.directory].has_key(i):
			#	sequence[name.directory]=sequence[name.directory]+snpbases[name.directory][i]
			#else:
			#	sequence[name.directory]=sequence[name.directory]+refbases[i]
		
		
#		if numsnpbases>1:
	#		print >> snpgraph, '1'
		print >> output, outstring+'\t'+pattern
		if tabfile!='':
                       	print >> tabout, tabstring+'"\nFT                   /colour='+str(numsnpbases)

		#else:
		#	print >> snpgraph, '0'
		#	for name in snpfiles:
		#		sequence[name.directory]=sequence[name.directory][:-1]
		#		goodhits[name.directory]=goodhits[name.directory]-1
		
	#	print >> covgraph, str(float(avcoverage[poscount])/len(snpfiles))
	#	print >> covcountgraph, covcount[poscount]
	#	poscount=poscount+1


#while reflen>=poscount:
#	print >> snpgraph, '0'
#	print >> covgraph, str(float(avcoverage[poscount])/len(snpfiles))
#	print >> covcountgraph, str(covcount[poscount])
#	poscount=poscount+1


if align!='':
	
	print >>alnfile, '>Reference\n'+refsequence
		
	for name in snpfiles:
		if goodhits[name.directory]>(len(sequence[name.directory])/2):
			print >> alnfile, '>'+name.directory.replace('pool','').replace('.snps','')+'\n'+sequence[name.directory]		

				
	alnfile.close()
	
if tabfile!='':
	tabout.close()			

snpgraph.close()
covgraph.close()
covcountgraph.close()
output.close()

print "Done.\n"
