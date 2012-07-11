#!/usr/bin/env python
import string, re, gzip
import os, sys, getopt, random, math


def Usage():
	print 'plots2heatmap.py Usage:'
	print 'Transfers a set of plots to a heatmap'
	print 'plots2heatmap.py [options] <input plot files>'
	print 'Options:'
	print '-c <value>\toptional cutoff for regions to include in heatmap (makes the heatmap into presence/absence plot)'
	print '-b\t\twill change from selecting regions above the cutoff to regions below the cutoff'
	print '-o <filename>\toutput file name'
	print '-h\t\tshow this help'
	print 'Written by Simon R Harris, Wellcome Trust Sanger Institute, UK. 2009'



##############################
# Get command line arguments #
##############################


def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "ho:c:b", ["ref=", "out="])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	outfile=''
	inputdirs=[]
	abovebelow='a'
	cutoff=''

	for opt, arg in opts:
	
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-o", "--out"):
			outfile=arg
		elif opt in ("-c", "--cutoff"):
			try: cutoff=float(arg)
			except ValueError:
				print 'Non-numeric value %s found for cutoff' % argv
		elif opt in ("-b", "--below"):
			abovebelow='b'


	inputdirs=args
	
	if inputdirs==[]:
		print 'Error: No input files selected!'
		Usage()
		sys.exit()
	if outfile=='':
		print 'Error: No output file specified!'
		Usage()
		sys.exit()
		
	overwrite='n'
	
	while overwrite!='y' and os.path.isfile(outfile):
		overwrite=raw_input(outfile+' exists! Overwrite? (y/n): ')
		overwrite.lower()
		if overwrite=='n':
			outfile=raw_input('Enter a new output file name: ')
	
	
	
	return inputdirs, outfile, cutoff, abovebelow




if __name__ == "__main__":
	argv=sys.argv[1:]
	inputdirs, outfile, cutoff, abovebelow=getOptions(argv)

plotfiles=[]
for y, plotfile in enumerate(inputdirs):

	if not os.path.isfile(plotfile):
		print "Cannot find file:", plotfile

	if plotfile.split('.')[-1]=='gz':
		plotfiles.append(gzip.open(plotfile, 'r'))
	else:
		plotfiles.append(open(plotfile, 'rU'))

length=len(plotfiles[0].readlines())

output=open(outfile,"w")
count=0
total=0

for x in range(0,length):
	count=count+1
	if count>=10000:
		total=total+count
		count=0
		print str(int(100*(float(total)/length)))+'% complete\r',
		sys.stdout.flush()
	for plotfile in plotfiles:
		try: line=plotfile.next().strip()
		except:
			continue
		
		value=float(line.strip().split()[0])
		
		if cutoff=='':
			print >> output, value,
		elif abovebelow=='a' and value>cutoff:
			print >> output, "1",
		elif abovebelow=='a':
			print >> output, "0",
		elif abovebelow=='b' and value<cutoff:
			print >> output, "1",
		elif abovebelow=='b':
			print >> output, "0",
			
	print >> output, '\n',

print "100% complete"
output.close()
		