#!/usr/bin/env python
#/usr/bin/env python
import string
import os, sys
import pysam


#pysam cigar states: 0 = Mapped, 1 = Insertion, 2 = Deletion, 3 = N, 4 = S, 5 = H, 6 = P


_proc_status = '/proc/%d/status' % os.getpid()

_scale = {'kB': 1024.0, 'mB': 1024.0*1024.0,
          'KB': 1024.0, 'MB': 1024.0*1024.0}

def _VmB(VmKey):
    '''Private.
    '''
    global _proc_status, _scale
     # get pseudo file  /proc/<pid>/status
    try:
        t = open(_proc_status)
        v = t.read()
        t.close()
    except:
        return 0.0  # non-Linux?
     # get VmKey line e.g. 'VmRSS:  9999  kB\n ...'
    i = v.index(VmKey)
    v = v[i:].split(None, 3)  # whitespace
    if len(v) < 3:
        return 0.0  # invalid format?
     # convert Vm value to bytes
    return float(v[1]) * _scale[v[2]]


def memory(since=0.0):
    '''Return memory usage in bytes.
    '''
    return _VmB('VmSize:') - since


def resident(since=0.0):
    '''Return resident memory usage in bytes.
    '''
    return _VmB('VmRSS:') - since


def stacksize(since=0.0):
    '''Return stack size in bytes.
    '''
    return _VmB('VmStk:') - since
 
 
 
def get_mapped_length_from_cigar(cigarstring):
	ciglength=0
	for part in cigarstring:
		if part[0] in [0,2,3,6]:
			ciglength+=part[1]
	return ciglength

def trim_read_start(origread, newread):

	#sort out trim at start of reads
	
	#find any trim at the start of the realigned read from the cigar string
	starttrim=0
	if newread[1][0][0]==4:
		starttrim+=newread[1][0][1]
	#find any trim at the start of the original read from the cigar string
	startbtrim=0
	if origread.cigar[0][0]==4:
		startbtrim+=origread.cigar[0][1]
	
	#calculate the difference in trim length, which should be the difference in start position if nothing else has changed
	startdiff=starttrim-startbtrim
	#calculate what the start of the realigned read would have been with the original trim value
	newstart=newread[0]-startdiff
	
	
	#if the change in trim does make the start positions of the two reads is the same, we can change the read start and cigar string of the realigned read
	
	if startdiff>0 and newstart==origread.pos:
	
		#in any case, the start position can be changed
	
		newread[0]=origread.pos
		
		
		#for cases where both original and realigned reads are trimmed at the start (almost done)
		if starttrim>0 and startbtrim>0 and startdiff!=0:
			
			if newread[0]:# in startlocs:
				#if after the trimmed part, the realigned read is mapped (state 0), then we need to add the new, untrimmed part to the mapped part and replace the second value in the cigar string
				if len(newread[1])>1 and newread[1][1][0]==0:
					newread[1][0]=origread.cigar[0]
					newread[1][1]=(0,newread[1][1][1]+startdiff)
				#if after the trimmed part, the realigned read is mapped (state 0), need to do something - can this ever happen???
				elif len(newread[1])>1:
					print "if this happens, report to sh16@sanger.ac.uk. Error code: s_realign_dindel_Err1"
					print "err", newread, origread.cigar
					sys.exit()
#						else:
#							print "if this happens, report to sh16@sanger.ac.uk. Error code: realign_dindel_Err1b"
				
			
				
		#for cases where the realigned read has a trimmed start, but the original alignment does not (done)
		elif starttrim>0 and startbtrim==0:
			if len(newread[1])>1 and newread[1][1][0]==0:
				combinedlen=newread[1][0][1]+newread[1][1][1]
				newcigar=newread[1][1:]
				newread[1]=newcigar
				newread[1][0]=(0,combinedlen)
			elif len(newread[1])>1:
				newread[1][0]=(0,newread[1][0][1])


		#for cases where the original alignment has a trimmed start, but the realignment does not (not done), need to do something - does this ever happen???				
#		elif startbtrim==0 and starttrim>0:
#			print "if this happens, report to sh16@sanger.ac.uk. Error code: s_realign_dindel_Err2"
#			sys.exit()

		#for cases where neither the original or realigned reads havve a trimmed start, leave them alone
		
	#return new position and cigar for read
	return newread[0], newread[1]




def trim_read_end(origread, newread):
	#sort out trim at end of reads (not done)
					
	endtrim=0
	if newread[1][-1][0]==4:
		endtrim+=newread[1][-1][1]
	endbtrim=0
	if origread.cigar[-1][0]==4:
		endbtrim+=origread.cigar[-1][1]
	
	
	origreadend=origread.pos+get_mapped_length_from_cigar(origread.cigar)
	newreadend=newread[0]+get_mapped_length_from_cigar(newread[1])
	
	#find any trim at the start of the realigned read from the cigar string
	endtrim=0
	if newread[1][-1][0]==4:
		endtrim+=newread[1][-1][1]
	#find any trim at the start of the original read from the cigar string
	endbtrim=0
	if origread.cigar[-1][0]==4:
		endbtrim+=origread.cigar[-1][1]
	
	#calculate the difference in trim length, which should be the difference in start position if nothing else has changed
	enddiff=endtrim-endbtrim
	#calculate what the start of the realigned read would have been with the original trim value
	newend=newreadend+enddiff
	
	#print origreadend, newreadend, endtrim, endbtrim, enddiff, newend\
	
	#if the change in trim does make the end positions of the two reads the same, we can change the read end and cigar string of the realigned read
	if enddiff>0 and newend==origreadend:
		#sys.exit()
	#return newread[0], newread[1]
	
	
		
		#for cases where both original and realigned reads are trimmed at the start (almost done)
		if endtrim>0 and endbtrim>0 and enddiff!=0:
			
			#if newread[0] in startlocs:
			#if before the trimmed part, the realigned read is mapped (state 0), then we need to add the new, untrimmed part to the mapped part and replace the second to last value in the cigar string
			if len(newread[1])>1 and newread[1][-2][0]==0:
				#print "b", newread, origread.cigar
				newread[1][-1]=origread.cigar[-1]
				newread[1][-2]=(0,newread[1][-2][1]+enddiff)
				#print "a", newread
			#if after the trimmed part, the realigned read is mapped (state 0), need to do something - can this ever happen???
			elif len(newread[1])>1:
				print "if this happens, report to sh16@sanger.ac.uk. Error code: e_realign_dindel_Err1"
				print "err", newread, origread.cigar
				sys.exit()
#		else:
#			print "if this happens, report to sh16@sanger.ac.uk. Error code: realign_dindel_Err1b"
				
			
				
		#for cases where the realigned read has a trimmed end, but the original alignment does not (done)
		elif endtrim>0 and endbtrim==0:
			if len(newread[1])>1 and newread[1][-2][0]==0:
				#print "b", newread, origread.cigar
				combinedlen=newread[1][-1][1]+newread[1][-2][1]
				newcigar=newread[1][:-1]
				newread[1]=newcigar
				newread[1][-1]=(0,combinedlen)
				#print "a", newread
			elif len(newread[1])>1:
				newread[1][0]=(0,newread[1][0][1])


		#for cases where the original alignment has a trimmed start, but the realignment does not (not done), need to do something - does this ever happen???				
#		elif endtrim==0 and endbtrim>0:
#			print "if this happens, report to sh16@sanger.ac.uk. Error code: e_realign_dindel_Err2"
#			sys.exit()

		#for cases where neither the original or realigned reads havve a trimmed start, leave them alone
		
	#return new position and cigar for read
	return newread[0], newread[1]

    
    

if len(sys.argv)!=5 or "-h" in sys.argv:
	print "~sh16/scripts/realign_dindel_bams.py <dindel bam> <dindel vcf> <original bam> <output bamfile>"
	sys.exit()


if sys.argv[1].split(".")[-1]=="bam":
	samfile = pysam.Samfile( sys.argv[1], "rb" )
elif sys.argv[1].split(".")[-1]=="sam":
	samfile = pysam.Samfile( sys.argv[1], "r" )
else:
	print "Not a bam file"
	sys.exit()


	
refs=samfile.references
lengths=samfile.lengths


reads={}

count=0
for read in samfile:

	if not read.is_unmapped:
		if not refs[read.rname] in reads:
			reads[refs[read.rname]]={}
		if not read.qname in reads[refs[read.rname]]:
			reads[refs[read.rname]][read.qname]={}
		if read.is_reverse:
			reads[refs[read.rname]][read.qname]["r"]=[read.pos, read.cigar]
			count+=1
		else:
			reads[refs[read.rname]][read.qname]["f"]=[read.pos, read.cigar]
			count+=1



samfile.close()
		
for ref in reads.keys():
	print len(reads[ref])

print memory()/1000000, resident()/1000000, stacksize()


try:
	vcffile=open(sys.argv[2],"rU")
except StandardError:
	print "Could not open dindel vcf file"
	sys.exit()

locs=[]

for line in vcffile:
	if len(line)==0 or line[0]=="#":
		continue
	words=line.strip().split()
	indelstart=int(words[1])
	indelend=indelstart+len(words[3])
	locs.append([indelstart, indelend])
	



if sys.argv[3].split(".")[-1]=="bam":
	samfileb = pysam.Samfile( sys.argv[3], "rb" )
elif sys.argv[3].split(".")[-1]=="sam":
	samfileb = pysam.Samfile( sys.argv[3], "r" )
else:
	print "Not a bam file"
	sys.exit()


updatepairs={}


refsb=samfileb.references

count=0
for read in samfileb:
	
	if not read.is_unmapped and refsb[read.rname] in reads:
		
		if read.qname in reads[refsb[read.rname]]:
			if read.is_reverse and "r" in reads[refsb[read.rname]][read.qname]:
				
				#sort out trimming at the start of the read
				reads[refsb[read.rname]][read.qname]["r"][0], reads[refsb[read.rname]][read.qname]["r"][1]=trim_read_start(read, reads[refsb[read.rname]][read.qname]["r"])
				
				#sort out trimming at the end of the read
				reads[refsb[read.rname]][read.qname]["r"][0], reads[refsb[read.rname]][read.qname]["r"][1]=trim_read_end(read, reads[refsb[read.rname]][read.qname]["r"])
				
				
				#update the read positions
				if read.pos!=reads[refsb[read.rname]][read.qname]["r"][0] or read.cigar!=reads[refsb[read.rname]][read.qname]["r"][1]:
					if get_mapped_length_from_cigar(reads[refsb[read.rname]][read.qname]["r"][1])>0:
						read.cigar=reads[refsb[read.rname]][read.qname]["r"][1]
						read.pos=reads[refsb[read.rname]][read.qname]["r"][0]
					else:
						read.is_unmapped=True
						read.is_proper_pair=False
						read.isize=0
						read.rname=-1
						read.pos=-1
						read.mapq=0
				
				if not read.qname in updatepairs:
					updatepairs[read.qname]={}
				updatepairs[read.qname]["r"]=read
				count+=1
			
				
			elif not read.is_reverse and "f" in reads[refsb[read.rname]][read.qname]:
				
				#sort out trimming at the start of the read
				reads[refsb[read.rname]][read.qname]["f"][0], reads[refsb[read.rname]][read.qname]["f"][1]=trim_read_start(read, reads[refsb[read.rname]][read.qname]["f"])
				
				#sort out trimming at the end of the read
				reads[refsb[read.rname]][read.qname]["f"][0], reads[refsb[read.rname]][read.qname]["f"][1]=trim_read_end(read, reads[refsb[read.rname]][read.qname]["f"])
				
				#Update the read positions
				if read.pos!=reads[refsb[read.rname]][read.qname]["f"][0] or read.cigar!=reads[refsb[read.rname]][read.qname]["f"][1]:
					if get_mapped_length_from_cigar(reads[refsb[read.rname]][read.qname]["f"][1])>0:
						read.cigar=reads[refsb[read.rname]][read.qname]["f"][1]
						read.pos=reads[refsb[read.rname]][read.qname]["f"][0]
					else:
						read.is_unmapped=True
						read.is_proper_pair=False
						read.isize=0
						read.rname=-1
						read.pos=-1
						read.mapq=0
						
				if not read.qname in updatepairs:
					updatepairs[read.qname]={}
				updatepairs[read.qname]["f"]=read
				count+=1

print memory()/1000000, resident()/1000000, stacksize()
print count
try:
	samout = pysam.Samfile(sys.argv[4], "wb", template=samfileb)
except StandardError:
	print "Failed to open output file"
	sys.exit()

count=0
done=set()
filterforproperpair=False
samfileb.reset()
for read in samfileb:
	if read.is_proper_pair or not filterforproperpair:
		if read.qname in updatepairs:
			if len(updatepairs[read.qname])==1:
				if read.is_reverse and not "r" in updatepairs[read.qname]:
					updatepairs[read.qname]["r"]=read
				elif not read.is_reverse and not "f" in updatepairs[read.qname]:
					updatepairs[read.qname]["f"]=read
			if len(updatepairs[read.qname])==2:
			
				if updatepairs[read.qname]["f"].is_unmapped:
					updatepairs[read.qname]["r"].mate_is_unmapped
					updatepairs[read.qname]["r"].isize=0
				if updatepairs[read.qname]["r"].is_unmapped:
					updatepairs[read.qname]["f"].mate_is_unmapped
					updatepairs[read.qname]["f"].isize=0
				
				if not updatepairs[read.qname]["f"].is_unmapped and not updatepairs[read.qname]["r"].is_unmapped:
					isize=updatepairs[read.qname]["r"].aend-updatepairs[read.qname]["f"].pos
					updatepairs[read.qname]["f"].isize=isize
					updatepairs[read.qname]["r"].isize=isize*-1
			
				for fr in updatepairs[read.qname]:
	#				if updatepairs[read.qname][fr].is_unmapped:
	#					print updatepairs[read.qname][fr]
	#					print updatepairs[read.qname][fr].mapq, updatepairs[read.qname][fr].pos, updatepairs[read.qname][fr].cigar
					
					samout.write(updatepairs[read.qname][fr])
					count+=1
				
				del updatepairs[read.qname]
				done.add(read.qname)
		elif read.qname not in done:
			samout.write(read)
					
					
print count			



    	
		
