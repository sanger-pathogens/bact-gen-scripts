#!/usr/bin/env python

#Removes contaminants from MIRA input files



#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------

import string, re, math
import os, sys, getopt


def Usage():
	print '\nremove_contaminants.py Usage:\n'
	print ''
	print '\nWritten by: Simon R Harris, Wellcome Trust Sanger Institute, UK. 2008\n'



#------------------------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------------------------

def getOptions(arg):

        try:
                opts, args = getopt.getopt(argv, "hf:q:x:b:i:l:o:", ["help", "dirty" "embl=" "ref=", "fasta=", "fastq=", "xml=", "blast=", "percentid=", "length=", "output="])
        except getopt.GetoptError:
		print "Option Error!", argv
                Usage()
                sys.exit(2)

        fasta=''
        fastq=''
        xml=''
	blast=''
	percentid=90
	length=80
	outfile=''
	
        for opt, arg in opts:
		
                if opt in ("-h", "--help"):
                        Usage()
                        sys.exit()
                elif opt in ("-f", "--fasta"):
                        fasta=arg
                elif opt in ("-q", "--fastq"):
                        fastq=arg
		elif opt in ("-x", "--xml"):
                        xml=arg
		elif opt in ("-b", "--blast"):
                        blast=arg
		elif opt in ("-i", "--percentid"):
                        percentid=float(arg)
		elif opt in ("-l", "--length"):
                        length=int(arg)
		elif opt in ("-o", "--output"):
                        outfile=arg
		
	infiles=args
	

        if fasta=='':
                print 'no input fasta file selected'
                Usage()
                sys.exit()
	if fastq=='':
		print 'no input fastq file selected. Setting to '+fasta+'.qual'
		fastq=fasta+'.qual'
	if xml=='':
		print 'no input xml file selected'
                Usage()
                sys.exit()
	if xml=='':
		print 'no blast file selected'
                Usage()
                sys.exit()

        return fasta, fastq, xml, blast, percentid, length, outfile




if __name__ == "__main__":
        argv=sys.argv[1:]
        fasta, fastq, xml, blast, percentid, length, outfile=getOptions(argv)
	

lines=open(blast,'rU').readlines()

f=open(fasta, 'rU').read()
q=open(fastq, 'rU').read()
x=open(xml, 'rU').read()


fout=open(fasta+'.new', 'w')
qout=open(fastq+'.new', 'w')
xout=open(xml+'.new', 'w')

cout=open('contaminants.fasta', 'w')

linecount=len(lines)
count=0
for line in lines:
	count=count+1
	words=line.strip().split()
	start=string.find(f,words[0])
	if start!=-1:
	
		seqlen=len(f[start-1:].split('\n')[1].strip())
		if ((float(words[3])/float(seqlen))*100<length) or float(words[2])<percentid:
			continue
	
		print >> fout, f[:start-1]
		print >> cout, '\n'.join(f[start-1:].split('\n')[:2])
		f='\n'.join(f[start-1:].split('\n')[2:])
	start=string.find(q,words[0])
	if start!=-1:
		print >> qout, q[:start-1]
		print >> cout, '\n'.join(q[start-1:].split('\n')[:2])
		q='\n'.join(q[start-1:].split('\n')[2:])
	start=string.find(x,words[0])
	if start!=-1:
		print >> xout, x[:start-23]
		print >> cout, '\n'.join(x[start-23:].split('\n')[:7])
		x='\n'.join(x[start-23:].split('\n')[7:])
	
	fout.flush()
	qout.flush()
	xout.flush()
	cout.flush()
	
	print str((float(count)/float(linecount))*100)+'% done...\r',
	sys.stdout.flush()

print >> fout, f
print >> qout, q
print >> xout, x

fout.close()
qout.close()
xout.close()
cout.close()

