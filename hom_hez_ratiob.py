#!/usr/bin/env python
import string, re, gzip
import os, sys


if len(sys.argv)!=3 or '-h' in sys.argv[1:]:
	print "hom_hez_ratio.py <all.pileup file (may be zipped)> <all.snp file (may be zipped)>"
	sys.exit()


translate={"A":7, "G":8, "C":9, "T":10}

#print "Name\tTot_Homozygous\tTot_Heterozygous\tTot_Het/Hom\tSNP_Homozygous\tSNP_Heterozygous\tSNP_Het/Hom\to-e/e"


if sys.argv[1].split('.')[-1]=='gz':
	lines=gzip.open(sys.argv[1], 'r').readlines()
else:
	lines=open(sys.argv[1], 'rU').readlines()

homcount=0
hezcount=0
SNPcount=0.0
othercount=0.0
ratio=0.0

for line in lines:
	if line.split()[3]=="0":
		continue
	maximum=0
	total=0
	for base in [5,6,7,8]:
		if int(line.split()[base])>maximum:
			maximum=int(line.split()[base])
		total=total+int(line.split()[base])
	if total==0:
		continue
	SNPcount=SNPcount+(float(maximum)/total)
	othercount=othercount+(float(total-maximum)/total)
	if total-maximum==0:
		homcount=homcount+1
	else:
		hezcount=hezcount+1

if homcount>0:
	hethez=str(float(hezcount)/homcount)
else:
	hethez="inf"

if SNPcount>0:
	#ratio=str(float(othercount)/SNPcount)
	ratio=str(othercount/(hezcount+homcount))
	oe=str(((float(othercount)/SNPcount)-0.01)/0.01)
else:
	ratio="inf"
	oe="-"

#print sys.argv[1]+"\t"+str(homcount)+"\t"+str(hezcount)+"\t"+hethez+"\t"+ratio+"\t"+oe


tothomcount=homcount
tothezcount=hezcount
tothethez=hethez
totratio=ratio
totoe=oe

if sys.argv[2].split('.')[-1]=='gz':
	lines=gzip.open(sys.argv[2], 'r').readlines()
else:
	lines=open(sys.argv[2], 'rU').readlines()

homcount=0
hezcount=0
SNPcount=0.0
othercount=0.0

for line in lines:
	if line.split()[0]=="SNP_hom:":
		homcount=homcount+1
		#SNPcount=SNPcount+int(line.split()[translate[line.split()[6]]])
		SNPcount=SNPcount+1
	elif line.split()[0]=="SNP_hez:":
		hezcount=hezcount+1
		maximum=0
		total=0
		for base in [7,8,9,10]:
			if int(line.split()[base])>maximum:
				maximum=int(line.split()[base])
			total=total+int(line.split()[base])
		if total==0:
			continue
		SNPcount=SNPcount+(float(maximum)/total)
		othercount=othercount+(float(total-maximum)/total)


if homcount>0:
	hethez=str(float(hezcount)/homcount)
else:
	hethez="inf"

if SNPcount>0:
	#ratio=str(float(othercount)/SNPcount)
	ratio=str(othercount/(hezcount+homcount))
	oe=str(((float(othercount)/SNPcount)-0.01)/0.01)
else:
	ratio="inf"
	oe="-"

print sys.argv[2]+"\t"+str(tothomcount)+"\t"+str(tothezcount)+"\t"+tothethez+"\t"+str(homcount)+"\t"+str(hezcount)+"\t"+hethez+"\t"+str((float(hethez)-float(tothethez))/float(tothethez))+"\t"+totratio+"\t"+ratio+"\t"+str((float(ratio)-float(totratio))/float(totratio))
		