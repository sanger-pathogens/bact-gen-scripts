#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math
import time

def splitfile(readfile, run):

	tiles={}
	tilekeys=[]
	count=0
	for read in readfile:
		tile=read.split(":")[2]
		if not tiles.has_key(tile):
			tiles[tile]=[0,0,0,0,0,0,0,0,0,0,0,0,0]
			tilekeys.append(tile)
		
		tag=read.split("#")[1].split('/')[0]
		count=count+1
		tag=read.split("#")[1].split('/')[0]
		tiles[tile][int(tag)]=tiles[tile][int(tag)]+1
		for y in range(0,3):
			readfile.next()
	
	output=open(run+"_tile_counts.out" ,"w")
	print count
	for tile in tilekeys:
		print >> output, tile,
		for x in tiles[tile]:
			print >> output, "\t"+str(x)
		print >> output, "\n",
	output.close()
	



if ((len (sys.argv)!=2) and (len (sys.argv)!=3) and (len(sys.argv)!=4)) or '-h' in sys.argv[1:]:
	print "split_solexa.py <fastq file to split>"
	sys.exit()

start=time.time()
fr=sys.argv[1].split('_')[-1].split('.')[0]


if len(sys.argv)==4 and sys.argv[3].lower()=='s':
	stupid='y'
else:
	stupid='n'
	

readfile=open(sys.argv[1], "rU")

#nlines=0
#print "Counting lines..."
##for line in readfile:
##	nlines += 1
#
#nlines=int(os.popen("wc -l "+sys.argv[1]).read().strip().split()[0])
#
#
#print nlines/4, "sequences to extract..."

print sys.argv[1]+":"
run="_".join(sys.argv[1].split('/')[-1].split('_')[:2])

splitfile(readfile,run)

		
print time.time()-start

