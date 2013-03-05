#!/usr/bin/env python

import os, sys

lines=open(sys.argv[1], "rU").read().split(">")[1:]

for line in lines:
	words=line.split('\n')
	name=words[0].strip().split()[0]
	seq=''.join(words[1:])
	print name, len(seq)
