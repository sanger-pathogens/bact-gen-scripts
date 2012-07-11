#!/usr/bin/env python

import os, sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules']))
from Si_nexus import *
from Si_general import *
from Si_SeqIO import *
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import pylab

tabfile=open(sys.argv[1], "rU")

SNPs=[]
for line in tabfile:
    words= line.strip().split()
    if len(words)==3 and words[1]=="SNP":
       SNPs.append({})
       SNPs[-1]["location"]=words[2]
    elif len(words)>1 and words[1][0]=="/":
        qualifier=words[1][1:].split("=")[0].strip()
        value=' '.join(words[1:]).replace("/"+qualifier+"=","")
        SNPs[-1][qualifier]=value.replace('"',"")



locations=set([])
homoplasylocations=set([])
nonhomlocations=set([])
SNPcount=0
homoplasycount=0

print '\t'.join(['Number of SNP bases', 'Number of non-homoplasy bases', 'number of homoplasy bases', '% of bases that are homoplasies', 'Total number of SNPs', 'Total number of non-homoplasic SNPs', 'Total number of homoplasic SNPs', '% of SNPs that are homoplasic']) 

for SNP in SNPs:

    homoplasy=False
    if "homoplasy" in SNP:
        homoplasy=True
        
    
    if not int(SNP["location"]) in locations:
        locations.add(int(SNP["location"]))
    SNPcount+=1

    if homoplasy and not int(SNP["location"]) in homoplasylocations:
        homoplasylocations.add(int(SNP["location"]))
        
    elif homoplasy:
        homoplasycount+=1
    elif not int(SNP["location"]) in nonhomlocations:
        nonhomlocations.add(int(SNP["location"]))
        


print '\t'.join(map(str,[len(locations), len(locations)-len(homoplasylocations), len(homoplasylocations), (float(len(homoplasylocations))/len(locations))*100, SNPcount, SNPcount-homoplasycount, homoplasycount, (float(homoplasycount)/SNPcount)*100]))

print 'Note: For total non-homoplasies, total homoplasies and % of SNPs that are homoplasic, for each site where homoplasy occurs, one SNP is considered non-homoplasic and the remainder are counted as homoplasic'


legend_colours=[(pylab.Rectangle((0, 0), 1, 1, fc="b")), (pylab.Rectangle((0, 0), 1, 1, fc="r"))]
n, bins, patches = plt.hist(map(list,[nonhomlocations, homoplasylocations]), bins=50, color=["b", "r"], histtype="barstacked")
plt.legend(legend_colours, ["Non homoplasy sites", "Homoplasy sites"])
plt.xlabel("Genome position")
plt.ylabel("SNP Frequency")
plt.xlim(0,1100000)
plt.show()


#print missinglist
