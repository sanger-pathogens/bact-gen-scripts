#!/usr/bin/env python3

import sys
for line in open(sys.argv[1], "rU"):
	words=line.strip().split(",")
	for x, word in enumerate(words):
		print(x+1, word)
	sys.exit()

		
