#!/usr/bin/env python

import string
import os, sys



for filename in sys.argv[1:]:
	try:
		lines=open(filename, "rU")
	except IOError:
		print "Could not open log file:", filename
		continue
	failed=False
	for logline in lines:
		if logline=="":
			continue
		try:
			logcommand=(' '.join(logline.strip().split()[:-1])).strip()
		except StandardError:
			logcommand=''
		try:
			logreturn=int(logline.strip().split()[-1])
		except StandardError:
			continue
		
		if logreturn!=0:
			print filename, "failed at command", logcommand, "with error code", logreturn
			failed=True
			break
	if not failed:
		print filename, "succeeded"
		