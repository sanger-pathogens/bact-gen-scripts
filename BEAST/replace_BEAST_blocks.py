#!/usr/bin/env python

import os, sys
import xml.etree.ElementTree as ET
from optparse import OptionParser


def main():
        usage = "usage: %prog [options] args"
        parser = OptionParser(usage=usage)
        
        
        parser.add_option("-o", "--output", action="store", dest="outputfile", help="output xml file name [default is to overwrite input xml]", type="string", metavar="FILE", default="")
        parser.add_option("-p", "--patterns", action="store", dest="patternsfile", help="File containing constant site patterns (can be created with ~sh16/scripts/create_beast_alignment.py", type="string", metavar="FILE", default="")
        parser.add_option("-x", "--xml", action="store", dest="xmlfile", help="xml file to edit", default="", type="string", metavar="FILE")
	parser.add_option("-m", "--noMLE", action="store_false", dest="mle", help="Do not add marginal likelihood estimation for model comparison to the end of the beast block [default is to add mle block]", default=True)

	return parser.parse_args()

# in-place prettyprint formatter

def indent(elem, level=0):
	i = "\n" + level*"	"
	if len(elem):
		if not elem.text or not elem.text.strip():
			elem.text = i + "	"
		if not elem.tail or not elem.tail.strip():
			elem.tail = i
		for elem in elem:
			indent(elem, level+1)
		if not elem.tail or not elem.tail.strip():
			elem.tail = i
	else:
		if level and (not elem.tail or not elem.tail.strip()):
			elem.tail = i


if __name__=="__main__":

	(options, args)=main()
	
	xmlfile=options.xmlfile
	patternsfile=options.patternsfile
	outputfile=options.outputfile
	add_model_comparison_block=options.mle
	
	
	if xmlfile=="":
		print "No xml file selected. For help use -h"
		sys.exit()
	
	if outputfile=="":
		outputfile=xmlfile
	
	
	if patternsfile=="" and not add_model_comparison_block:
		print "No patterns file provided, and no model comparison block to be added. Nothing to do... exiting"
		sys.exit()
	
	if not os.path.isfile(xmlfile):
		print "Cannot find file", xmlfile
		sys.exit()
		
	try:
		tree=ET.parse(xmlfile)
	except StandardError:
		print "Failed to parse xml file", xmlfile
		sys.exit()
	
	try:
		doc=tree.getroot()
		mcmc=doc.find('mcmc')
		log=mcmc.find('logTree')
		filenameprefix='.'.join(log.attrib["fileName"].split('.')[:-1])
	except StandardError:
		print "Failed to extract tree log file name from xml"
		sys.exit()
	
		
	if patternsfile!="":
	
		try:
			patternsf=open(patternsfile, "rU")
		except:
			print "Failed to open patterns file", patternsfile
			sys.exit()
		
		try:
			patternstext=patternsf.readlines()[0].strip()
		except StandardError:
			print "Is your patterns file empty?"
			sys.exit()
	
		#Check the patterns text
		try:
			patternsvalues=map(int,patternstext.split())
		except StandardError:
			print "Expecting first line of patterns file to be four integers separated by whitespace"
			print "This file can be created by running ~sh16/scripts/BEAST/prepare_beast_alignment.py"
			sys.exit()
		
		if len(patternsvalues)!=4:
			print "Expecting four values in first line of patterns file"
			sys.exit()
		
		patternsnode=doc.find('patterns')
		
		if patternsnode==None:
			print "No patterns block found."
			mergeblock=doc.find('mergePatterns')
			if mergeblock==None:
				print "No mergePatterns block found either. xml file misformed."
				sys.exit()
			else:
				print "Found a mergePatterns block, though. Perhaps you already replaced the patterns block in this file? It will not be replaced again."
			patternsfile=""
			if not add_model_comparison_block:
				print "No model comparison block to be added, so nothing to do... exiting"
				sys.exit()
		else:	
			patternsnode.tag="mergePatterns"
			
			patternsattribs={}
			for attrib in patternsnode.attrib:
				if attrib=="id":
					continue
				patternsattribs[attrib]=patternsnode.attrib[attrib]
			for attrib in patternsattribs:
				del patternsnode.attrib[attrib]
		
		
			alignmentnode=patternsnode.find("alignment")
			patternsnode.remove(alignmentnode)
		
			#create new subfeatures for mergepatterns block
			ET.SubElement(patternsnode, "patterns")
			patternssubnode=patternsnode[-1]
			for attrib in patternsattribs:
				patternssubnode.attrib[attrib]=patternsattribs[attrib]
			if not "from" in patternsattribs:
				patternssubnode.attrib["from"]="1"
			if not "every" in patternsattribs:
				patternssubnode.attrib["every"]="1"
			
			ET.SubElement(patternssubnode, "alignment")
			alignnode=patternssubnode[-1]
			alignnode.attrib["idref"]="alignment"
			
			ET.SubElement(patternsnode, "constantPatterns")
			constnode=patternsnode[-1]
			ET.SubElement(constnode, "alignment")
			alignbnode=constnode[-1]
			alignbnode.attrib["idref"]="alignment"
			ET.SubElement(constnode, "counts")
			countnode=constnode[-1]
			ET.SubElement(countnode, "parameter")
			paramnode=countnode[-1]
			paramnode.attrib["value"]=patternstext
		
	mle =doc.find("marginalLikelihoodEstimator")
	if mle!=None:
		print "This xml already includes a marginal likelihood estimator block, so will not add another one"
		add_model_comparison_block=False
		if patternsfile=="":
			print "No patterns to edit, so nothing to do... exiting"
			sys.exit()
	
	if add_model_comparison_block:
		
		#create the marginalLikelihoodEstimator subelement with attributes
		ET.SubElement(doc,"marginalLikelihoodEstimator")
		node=doc[-1]
		node.attrib["chainLength"]="1000000"
		node.attrib["pathSteps"]="100"
		node.attrib["pathScheme"]="betaquantile"
		node.attrib["alpha"]="0.30"
		
		#create the samplers subelement
		ET.SubElement(node, "samplers")
		samplersnode=node[-1]
		ET.SubElement(samplersnode, "mcmc")
		mcmcnode=samplersnode[-1]
		mcmcnode.attrib["idref"]="mcmc"
		
		#create the pathLikelihood subelement
		ET.SubElement(node, "pathLikelihood")
		pathlikenode=node[-1]
		pathlikenode.attrib["id"]="pathLikelihood"
		
		ET.SubElement(pathlikenode, "source")
		sourcenode=pathlikenode[-1]
		ET.SubElement(sourcenode, "posterior")
		posteriornode=sourcenode[-1]
		posteriornode.attrib["idref"]="posterior"
		ET.SubElement(pathlikenode, "destination")
		destnode=pathlikenode[-1]
		ET.SubElement(destnode, "prior")
		priornode=destnode[-1]
		priornode.attrib["idref"]="prior"
		
		#create the log subelement
		ET.SubElement(node, "log")
		lognode=node[-1]
		lognode.attrib["id"]="MLE"
		lognode.attrib["logEvery"]="1000"
		lognode.attrib["fileName"]=filenameprefix+".MLE.log"
		ET.SubElement(lognode, "pathLikelihood")
		plnode=lognode[-1]
		plnode.attrib["idref"]="pathLikelihood"
		
	
	indent(doc)
	
	#ET.dump(tree)
		
	tree.write(outputfile)
	print "Updated xml saved in", outputfile