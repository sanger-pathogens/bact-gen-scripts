#!/usr/bin/env python
import os, sys
file1=open(sys.argv[1])
lines2=open(sys.argv[2]).readlines()

newscores=[]
for x, line in enumerate(file1):
	score2=float(lines2[x].strip())
	score1=float(line.strip())
	if score2==0:
		newscores.append(0.0)
	else:
		newscores.append((score1/score2)*100)

outfile=open(sys.argv[3], "w")
for score in newscores:
	print >>outfile, score
	
outfile.close()