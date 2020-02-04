#!/usr/bin/env python

import os, sys
from numpy import exp, log
from scipy.stats import binom

def bounds(j, s, k=21, p=0.99):
	"""Python code to calculate the mash bounds for a particular jaccard distance (j), kmer (k) and sketch size (s) at a given probability (p)"""
	
	q2 = (1.0-p)/2.0	
	
	m2j = 1.0 / (2.0 * exp(k * j) - 1.0)
	
	x=0.0
	while x<s:
		cdfx=binom.cdf(x, s, m2j)
		if cdfx>q2:
			break
		x+=1
		
	if x>0:
		je = x/s
		j2m = -1.0 / k * log(2.0 * je / (1.0 + je))
		error = j2m - j
	else:
		error=float("Inf")
	
	return error
