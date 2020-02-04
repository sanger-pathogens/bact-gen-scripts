#!/usr/bin/env python

import re, string, os, tempfile, commands, subprocess, sys

f = open('encufasta.txt', 'r')
outfile = open('ecutest.txt', 'w')

fileString = f.read()
f.close()
string = re.split('(([>].*\n.*\n){1})',fileString)

print>>outfile, string


infile = open('ecutest.txt')
outfile = open('test.txt', 'w')
n = 0

list1 = infile.readlines()
for line in list1:
    line1 = line.split("'")
    if n < 1:       
        for chars in line1:
                if ">" in chars:
                    n += 1
                    #tempFile = tempfile.TemporaryFile()
                    #tempfilename = tempfile.mkstemp()
                    #tempfile = open(tempfilename, 'w')
                    #tempfile.write(chars)
                    f=tempfile.NamedTemporaryFile()
                    f.write(chars)
                    print f.name
                    s=os.popen('python ./runmemsat '+f.name)
                    print >> outfile, s
                    #s=subprocess.Popen('python ./runmemsat (tempfile)', shell=True)
                    #output = commands.getoutput(s)
                    #print>>outfile, output
       
            #n += 1
            #filename = "test%d.txt" % n
            #outfile = open(filename, 'w')
           # outfile.write(chars)
           # outfile.close



infile.close()


infile = open('test.memsat') 
outfile = open('memsatout1.txt', 'w')

for line in infile.readlines():
    
    if "nn" in line:
        
        line=line.replace(".","%")
    
        table=line.strip().split("%")
    
        proteinID=table[0]
        ext=table[1]
    
        print>> outfile, ">"+proteinID

infile = open('test.memsat')

for line in infile.readlines()[9:]:   
    
    if ": " in line:
       
       table=line.strip().split()
    
       accession=table[0]
       length=table[1]
       value=table[2]
       
       outfile = open('memsatout1.txt','a')
           
       print >>outfile, accession+" "+length+" "+value 
       
       
       
infile = open('memsatout1.txt') 
outfile = open('memsatout2.txt', 'w')

for line in infile.readlines():   
    
    if ">" in line:
    
        print>>outfile, line
    
    elif ": " in line:
       
       table=line.strip().split()
    
       accession=table[0]
       length=table[1]
       value=table[2]
       
       outfile = open('memsatout2.txt', 'a')
           
       print >>outfile, accession+" "+length+" "+value 
       
       
       
infile = open('memsatout2.txt') 
outfile = open('memsatout.txt', 'w')

for line in infile.readlines():
    
    if "(" not in line:
    
        print>>outfile, line    

infile = open('memsatout2.txt')

numlines = infile.readlines()[1:]

outfile = open('memsatout.txt', 'a')

print >>outfile, len(numlines)




f = open('memsatout.txt', 'r')

fileString = f.read()
f.close()
string = re.sub('\n\s*\n*', '\n', fileString)

f=open('memsatout.txt', 'w')f.write(string)f.close()

infile = open('memsatout.txt')
outfile = open('memsatfinal.txt', 'w')

empty=re.compile('^$')

for char in infile.read():
    
    char=char.replace('>','')
    char=char.replace('\n', ',')
    
    outfile.write(char)
    outfile.close








