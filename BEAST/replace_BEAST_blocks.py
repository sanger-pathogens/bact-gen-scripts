#!/usr/bin/env python

import os, sys
import xml.etree.ElementTree as ET
from optparse import OptionParser


def main():
        usage = "usage: %prog [options]"
        parser = OptionParser(usage=usage)
        
        
        parser.add_option("-o", "--output", action="store", dest="outputfile", help="output xml file name [default is to overwrite input xml, but I take no responsibility if it destroys your lovely xml.]", type="string", metavar="FILE", default="")
        parser.add_option("-p", "--patterns", action="store", dest="patternsfile", help="File containing constant site patterns (can be created with ~sh16/scripts/create_beast_alignment.py). Whitespace separated file containing on line with four columns for the number of A, C, G and T constant sites", type="string", metavar="FILE", default="")
        parser.add_option("-P", "--precision", action="store", dest="precisionfile", help="File containing precisions. File must be whitespace delimited with two columns: 1) Name of taxon, 2) precision (in days or years after the date specified in the xml file). e.g. If your samples are dated in days since 1900, but for one sample, Bob, you only know the year of isolation is 1989, you would need to define the date for Bob to be 32509 (days since 1900 for 1/1/1989) in beauti and include the line 'Bob 365' (i.e. Bob may have been isolated any day up to 365 days after 32509) in your input precision file for this script.", type="string", metavar="FILE", default="")
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
	precisionfile=options.precisionfile
	outputfile=options.outputfile
	add_model_comparison_block=options.mle
	
	
	if xmlfile=="":
		print "No xml file selected. For help use -h"
		sys.exit()
	
	if outputfile=="":
		outputfile=xmlfile
	
	
	if patternsfile=="" and precisionfile=="" and not add_model_comparison_block:
		print "No patterns or precision files provided, and no model comparison block to be added. Nothing to do... exiting"
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
			if not add_model_comparison_block and precisionfile=="":
				print "No model comparison block to be added and no precisions file specified, so nothing to do... exiting"
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
		if patternsfile=="" and precisionfile=="":
			print "No patterns or precisions to edit, so nothing to do... exiting"
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
	
	
	
	if precisionfile!="":
	
		try:
			precisionf=open(precisionfile, "rU")
		except:
			print "Failed to open precision file", precisionfile
			sys.exit()
		
		try:
			precisiontext=precisionf.readlines()
		except StandardError:
			print "Is your precision file empty?"
			sys.exit()
	
		#Check the patterns text
		precisiondict={}
		for line in precisiontext:
			line=line.strip()
			try:
				if line.split()[0] in precisiondict:
					print "WARNING:",line.split()[0],"appears more than once in precision fileI"
				precisiondict[line.split()[0]]=float(line.split()[1])
			except StandardError:
				print "Expecting precision file to contain an isolate name followed by a precision value. Instead found:"
				print line
				sys.exit()
		
		added_precision_to=[]
		if len(precisiondict)!=0:
			
			taxanode=doc.find('taxa')
			
			if taxanode==None:
				print "No taxa block found. xml file misformed."
				sys.exit()
			else:
				for taxonblock in taxanode.findall('taxon'):
					taxon_name=""
					for attrib in taxonblock.attrib:
						if attrib=="id":
							taxon_name=taxonblock.attrib["id"]
					if taxon_name=="":
						print "WARNING: Found no ID for taxon."
						continue
					dateblock=taxonblock.find("date")
					if dateblock==None:
						print "WARNING: Missing date block!"
						continue
					if taxon_name in precisiondict:
						if not "units" in dateblock.attrib:
							print "WARNING: Date block is missing units attribute!"
						else:
							print "Adding precision of", precisiondict[taxon_name], dateblock.attrib["units"], "to", taxon_name
							dateblock.attrib["precision"]=str(precisiondict[taxon_name])
							added_precision_to.append(taxon_name)

#			You then need to add:
#
#        <leafHeight taxon="D4Brazi82">
#                <parameter id="D4Brazi82.height"/>
#        </leafHeight>
#
#to the tree model,
				
			treeModelnode=doc.find('treeModel')
			
			if treeModelnode==None:
				print "No treeModel block found. xml file misformed."
				sys.exit()
			
			for taxon in added_precision_to:
				#create new subfeatures for taxa with precision
				ET.SubElement(treeModelnode, "leafHeight")
				leafHeightnode=treeModelnode[-1]
				leafHeightnode.attrib["taxon"]=taxon
				ET.SubElement(leafHeightnode, "parameter")
				parameternode=leafHeightnode[-1]
				parameternode.attrib["id"]=taxon+".height"
				
#        <randomWalkOperator windowSize="0.1" weight="1">
#                <parameter idref="D4Brazi82.height"/>
#        </randomWalkOperator>
#
#to the operators, and if you want:
			
			operatorsnode=doc.find('operators')
			
			if operatorsnode==None:
				print "No operators block found. xml file misformed."
				sys.exit()
			
			for taxon in added_precision_to:
				#create new subfeatures for taxa with precision
				ET.SubElement(operatorsnode, "randomWalkOperator")
				randomWalkOperator=operatorsnode[-1]
				randomWalkOperator.attrib["windowSize"]="0.1"
				randomWalkOperator.attrib["weight"]="1"
				ET.SubElement(randomWalkOperator, "parameter")
				parameternode=randomWalkOperator[-1]
				parameternode.attrib["idref"]=taxon+".height"
			

#        <parameter idref="D4Brazi82.height"/>
#
#to the log file. 
			
			mcmcnode=doc.find('mcmc')
			for lognode in mcmcnode.findall('log'):
				if lognode==None:
					print "No mcmc log block found. xml file misformed."
					sys.exit()
				
				if "id" in lognode.attrib and lognode.attrib['id']=='fileLog':
					for taxon in added_precision_to:
						#create new subfeatures for taxa with precision
						ET.SubElement(lognode, "parameter")
						parameternode=lognode[-1]
						parameternode.attrib["idref"]=taxon+".height"
					
				
	
	
	
	indent(doc)
	
	#ET.dump(tree)
		
	tree.write(outputfile)
	print "Updated xml saved in", outputfile
