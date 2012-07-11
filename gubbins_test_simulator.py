#!/usr/bin/env python
import string
import os, sys
import random


genome_length=1000000

genome_snps=int(sys.argv[1])

if genome_snps>genome_length:
	print "can't have more SNPs than genome length"
	sys.exit()

high_density_length=int(sys.argv[2])

high_snp_density=float(sys.argv[3])

high_density_snps=int(((float(genome_snps)/genome_length)*high_snp_density)*high_density_length)
print high_density_snps

low_density_length=int(sys.argv[4])

low_density_snps=int(sys.argv[5])

seq=[]

for x in xrange(0,genome_length):
	seq.append(random.choice("ACGT"))

output=open("test.aln", "w")
print >> output, ">A"
print >> output, ''.join(seq)
print >> output, ">B"
print >> output, ''.join(seq)


for x in random.sample(xrange(genome_length), genome_snps):
	if (x<100000 or x>100000+high_density_length) and (x<600000 or x>600000+low_density_length):
		seq[x]=random.choice("ACGT".replace(seq[x],""))
	
for x in random.sample(xrange(100000, 100000+high_density_length), high_density_snps):
	seq[x]=random.choice("ACGT".replace(seq[x],""))
	
for x in random.sample(xrange(600000, 600000+low_density_length), low_density_snps):
	seq[x]=random.choice("ACGT".replace(seq[x],""))

print >> output, ">C"
print >> output, ''.join(seq)
print >> output, ">D"
print >> output, ''.join(seq)


output.close()