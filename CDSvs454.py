#!/usr/bin/env python


#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------

import string, re
import os, sys, commands, getopt, gzip


def Usage():
	print 'A script to blast a set of CDSs against a set of contigs and output a tab file readable by Artemis'
        print 'CDSvs454.py Usage:'
        print 'CDSvs454.py -i [input fasta file of CDSs] -o [output file name] -d [fasta file containing contigs to blast against] -e [e-value cutoff] -p [%id cutoff] -l [log file name] -b [do blast?] {-h}'
        print 'or'
        print 'CDSvs454.py --in [input alignment] --out [output file name] --contigs [fasta file containing contigs to blast against] --evalue [e-value cutoff] --l [log file name] --blast [do blast] --help'
        print 'Copyright Simon R Harris, Wellcome Trust Sanger Institute, UK. 2008'



#------------------------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------------------------

def getOptions(arg):


        try:
                opts, args = getopt.getopt(argv, "hi:o:d:e:p:bl:", ["help", "in=", "out=", "contigs=", "evalue=", "pid=", "blast", "log="])
        except getopt.GetoptError:
                Usage()
                sys.exit(2)


        inputfile=''
        outfile=''
        db=''
        evalue='10'
        pid='0'
	runblast='n'
	logfile=''

        for opt, arg in opts:
                if opt in ("-h", "--help"):
                        Usage()
                        sys.exit()
                elif opt in ("-i", "--in"):
                        inputfile=arg
                elif opt in ("-o", "--out"):
                        outfile=arg
                elif opt in ("-d", "--contigs"):
                        db=arg
                elif opt in ("-e", "--evalue"):
                        evalue=arg
                elif opt in ("-p", "--pid"):
                        pid=arg
		elif opt in ("-b", "--blast"):
			runblast='y'
		elif opt in ("-l", "--log"):
			logfile=arg

        if inputfile=='':
                print 'No input fasta file selected!'
                Usage()
                sys.exit()
        if outfile=='':
                outfile=inputfile+".tab"
        if db=='':
		print 'No db fasta file selected!'
		Usage()
		sys.exit()

        return inputfile, outfile, db, evalue, pid, runblast, logfile


if __name__ == "__main__":
        argv=sys.argv[1:]
        inputfile, outfile, db, evalue, pid, runblast, logfile=getOptions(argv)

print "\n",

if runblast=='y':
	print "Creating blast database..."
	os.system('formatdb -p F -i '+db)
	print "Running blast..."
	os.system('blastall -p blastn -i '+inputfile+' -d '+db+' -e '+evalue+' -m 8 -o'+inputfile+'.blast')

CDSs={}
lines=open(inputfile,"rU").readlines()

for line in lines:
	if len(line)>0 and line[0]=='>':
		CDSs[line.split()[0][1:]]=int(line.split(':')[-2].split()[-1])

output=open(outfile,"w")

print >> output, "ID   CDSs\nFH   Key             Location/Qualifiers\nFH"

lines=open(inputfile+'.blast',"rU").readlines()
os.system('rm -f '+inputfile+'.blast')

lastquery=''
hitcount=[0,0,0,0,0]

print "Parsing blast output..."

for line in lines:
	words=line.strip().split()
	query=words[0]
	subject=words[1]
	if query!=lastquery:
		id=words[2]
		start=int(words[6])-1+CDSs[query]
		alnlen=words[3]
		e=words[10]
		bitscore=words[11]
		end=int(words[7])-1+CDSs[query]
		print >> output, 'FT   BLASTN_HIT      '+str(start)+'..'+str(end)
		if float(bitscore)>=200:
			colourstring='colour="255 0 0"'
			hitcount[0]=hitcount[0]+1
		elif float(bitscore)>=80:
			colourstring='colour="255 0 255"'
			hitcount[1]=hitcount[1]+1
		elif float(bitscore)>=50:
			colourstring='colour="0 255 0"'
			hitcount[2]=hitcount[2]+1
		elif float(bitscore)>=40:
			hitcount[3]=hitcount[3]+1
			colourstring='colour="0 0 255"'
		else:
			colourstring='colour="0 0 0"'
			hitcount[4]=hitcount[4]+1
		print >> output, 'FT                   /'+colourstring
		print >> output, 'FT                   /label="'+subject+'"'
		print >> output, 'FT                   /evalue="'+e+'"'
		print >> output, 'FT                   /identity="'+id+'"'
		print >> output, 'FT                   /bitscore="'+bitscore+'"'
		lastquery=query

print >> output, '//'
output.close()

nhits=len(CDSs.keys())
allhits=hitcount[0]+hitcount[1]+hitcount[2]+hitcount[3]+hitcount[4]

print '\nOf the',nhits,'CDSs Blasted:'
print '\t'+str(hitcount[0])+' (%.2f%%) had hits greater than or equal to a bitscore of 200' % (100*(float(hitcount[0])/nhits))
print '\t'+str(hitcount[1])+' (%.2f%%) had hits between 80 and 200' % (100*(float(hitcount[1])/nhits))
print '\t'+str(hitcount[2])+' (%.2f%%) had hits between 50 and 80' % (100*(float(hitcount[2])/nhits))
print '\t'+str(hitcount[3])+' (%.2f%%) had hits between 40 and 50' % (100*(float(hitcount[3])/nhits))
print '\t'+str(hitcount[4])+' (%.2f%%) had hits below 40' % (100*(float(hitcount[4])/nhits))
print '\t'+str(nhits-allhits)+' (%.2f%%) had no hits' % (100*(float(nhits-allhits)/nhits))

if logfile!='':
	output=open(logfile,"w")
	print >>output, '\nOf the',nhits,'CDSs Blasted:'
	print >>output, '\t'+str(hitcount[0])+' (%.2f%%) had hits greater than or equal to a bitscore of 200' % (100*(float(hitcount[0])/nhits))
	print >>output, '\t'+str(hitcount[1])+' (%.2f%%) had hits between 80 and 200' % (100*(float(hitcount[1])/nhits))
	print >>output, '\t'+str(hitcount[2])+' (%.2f%%) had hits between 50 and 80' % (100*(float(hitcount[2])/nhits))
	print >>output, '\t'+str(hitcount[3])+' (%.2f%%) had hits between 40 and 50' % (100*(float(hitcount[3])/nhits))
	print >>output, '\t'+str(hitcount[4])+' (%.2f%%) had hits below 40' % (100*(float(hitcount[4])/nhits))
	print >>output, '\t'+str(nhits-allhits)+' (%.2f%%) had no hits' % (100*(float(nhits-allhits)/nhits))
	output.close()

print "\nDone.\n"
