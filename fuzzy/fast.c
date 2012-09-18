/**************************************************************************
 *
 *-------------------------------------------------------------------------
 * Project Notes : various fasta and fastq manipulation routines
 *
 *
 *-------------------------------------------------------------------------
 #######################################################################
 # This software has been created by Genome Research Limited (GRL).    # 
 # GRL hereby grants permission to use, copy, modify and distribute    # 
 # this software and its documentation for non-commercial purposes     # 
 # without fee at the user's own risk on the basis set out below.      #
 # GRL neither undertakes nor accepts any duty whether contractual or  # 
 # otherwise in connection with the software, its use or the use of    # 
 # any derivative, and makes no representations or warranties, express #
 # or implied, concerning the software, its suitability, fitness for   #
 # a particular purpose or non-infringement.                           #
 # In no event shall the authors of the software or GRL be responsible # 
 # or liable for any loss or damage whatsoever arising in any way      # 
 # directly or indirectly out of the use of this software or its       # 
 # derivatives, even if advised of the possibility of such damage.     #
 # Our software can be freely distributed under the conditions set out # 
 # above, and must contain this copyright notice.                      #
 #######################################################################
 *
 *  Author : James C. Mullikin
 *
 *	Copyright (C) 1998-2000 by Genome Research Limited, All rights reserved.
 *
 **************************************************************************/


#include <math.h>
#include <values.h>
#include <stdio.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>
#include "fasta.h"

fasta *decodeFastq2 (char *fname, int *nContigs, long *tB, int qThresh, int basesonly)
{
	FILE *fil;
	char *data = NULL, *dp, *cp;
	long dataSize, totalbases=0, totalQbases=0;
	long i;
	int ret;
	int nSeg;
	fasta *seg,*segp;
	int isfasta = 0;

	if ((fil = fopen(fname,"r")) == NULL)
	{
		printf("error decodeFastq: %s not found\n",fname);
		return (NULL) ;
	}

	fseek( fil, 0, SEEK_END);
	dataSize = ftell( fil);
	rewind( fil);

	if(dataSize == 0)
	{
		printf("error decodeFastq: no data in %s\n",fname);
		fclose(fil);
		return(NULL);
	}

	if((data = (char *) malloc( (dataSize+1) * sizeof(char) )) == NULL)
	{
		printf("error decodeFastq: malloc error\n");
		fclose(fil);
		return(NULL);
	}
	
	for(i=0;i<dataSize;i+=1024*1000000){
		long size = 1024*1000000;
		if(dataSize-i < size)
			size = dataSize-i;
		//printf("1024*1000000 chunck %d\n",i);
		if(fread(data+i,sizeof(char), size, fil) != size)
		{
			printf("error decodeFastq: %s is too small\n",fname);
			fclose(fil);
			return (NULL) ;
		}
	}

	data[dataSize] = 0;
	nSeg = 0;
	dp = data;
	ret = 1;
	//printf("dataSize = %ld\n",dataSize);
	for(i=dataSize;--i>=0L;dp++)
	{
		if(*dp == '\n')
		{
			ret = 1;
			continue;
		}
		if(ret && (*dp == '@' || *dp == '>'))
		{
//			printf("%19ld %6d %70.70s\n",dataSize-i,nSeg,dp);
			nSeg++;
		}
		ret = 0;
	}
	//printf("i = %ld, dp = %ld\n",i,dp-data);

	if(nSeg == 0)
	{
		printf("error decodeFastq: no segments found\n");
		fclose(fil);
		return(NULL);
	}

	if((seg = (fasta *) calloc(nSeg, sizeof(fasta))) == NULL)
	{
		printf("error decodeFastq: calloc error\n");
		fclose(fil);
		return(NULL);
	}

	cp = dp = data;
	segp = seg;
	while(dp < data+dataSize)
	{
		if(*dp != '@' && *dp != '>')
		{
			dp++;
			continue;
		} else if ( dp > data && dp[-1] != '\n' ) {
			dp++;
			continue;
		} else {
			int tmp;
			int k;
			if(*dp == '@') isfasta = 0;
			if(*dp == '>') isfasta = 1;
			dp++;
			segp->name = (char *)(long int)(cp-data);
			while((*cp++ = tmp = *dp++) != ' ' && tmp != '\n');
			cp[-1] = 0;
			if(tmp != '\n')
			{
				segp->name2 = (char *)(long int)(cp-data);
				while((*cp++ = *dp++) != '\n');
				cp[-1] = 0;
			}
			else
				segp->name2 = NULL;
			segp->data = (char *)(long int)(cp-data);
			for(k=0;(tmp = *dp) && !((tmp == '+' || tmp == '>') && dp[-1] == '\n');dp++)
			{
				if(tmp == '\n')
					continue;
				*cp++ = tmp;
				k++;
			}
			segp->length = k;
			if(isfasta) {
				segp->finished = 1;
				segp->qual = NULL;
			} else {
			  segp->finished = 0;
			  k = 0;
			  if(tmp == '+')
			  {
				int cp_set = 0;
				int comp = 1;
				while(*dp++ != '\n');
				for(i=0;i<20;i++)
				{
					if(dp[i] == 0)
						break;
					if(dp[i] == ' ')
					{
						comp = 0;
						break;
					}
				}
				segp->qual = (char *)(long int)(cp-data);
				*cp = 0;
				if(comp) for(k=0;(tmp = *dp) && !(tmp == '@' && dp[-1] == '\n');dp++)
				{
					if(basesonly) continue;
					if(tmp == '\n')
						continue;
					*cp++ = tmp-041;
					k++;
				}
				else for(k=0;(tmp = *dp) && !(tmp == '@' && dp[-1] == '\n');dp++)
				{
					if(!basesonly) {
						if(tmp >= '0' && tmp <= '9')
						{
							if(*cp > 9) continue;
							*cp *= 10;
							*cp += tmp-'0';
							cp_set = 1;
							continue;
						}
						if(cp_set)
						{
							cp++;
							*cp = 0;
							cp_set = 0;
							k++;
						}
					}
				}
			  }
			  if(basesonly){
				  segp->finished = 1;
				  segp->qual = NULL;
			  } else if(segp->length != k) {
				  printf("Warning: %-20.20s (%d) has %d bases and %d quals, skipping entry\n",data+((long int)segp->name),strlen(data+((long int)segp->name)),segp->length,k);
				  cp = data+((long int)segp->name);
				  continue;
			  }
			}
			segp++;
			if(segp-seg > nSeg)
				break;
		}
	}
	nSeg = segp-seg;
	dataSize = cp - data;
	if(basesonly) {
		printf("decode: %ld %ld\n",data,seg);
		printf("data memloc = %ld\n",data);
		data = realloc(data,dataSize);
		printf("data memloc = %ld\n",data);
	}
	for(i=nSeg,segp=seg;--i>=0;segp++)
	{
		segp->name = data+((long int)segp->name);
		if(segp->name2)
			segp->name2 = data+((long int)segp->name2);
		segp->data = data+((long int)segp->data);
		segp->qual = data+((long int)segp->qual);
		if(segp->finished) segp->qual = NULL;
		//printf("segment %s length %d\n",segp->name,segp->length);
		totalbases+=segp->length;
		totalQbases += countQthresh(segp,qThresh);
	}
	printf("Total memory = %ld, total bases = %ld, Qbases(%d) = %ld, segments = %d\n",dataSize,totalbases,qThresh,totalQbases,nSeg);
	*nContigs = nSeg;
	*tB = totalbases;
	fclose(fil);
	return(seg);
}

fast *readFastq (char **argv, int argc, int *nContigs)
{
	FILE *fil;
	char *data = NULL, *dp, *cp;
	long dataSize;
	long i;
	int j,length;
	int ret;
	int nSeg,nseg;
	fast *seg,*segp;

	nseg = 0;
	for(j=0;j<argc;j++) 
        {
		if ((fil = fopen(argv[j],"r")) == NULL)
		{
			printf("error decodeFastq: %s not found\n",argv[j]);
			return (NULL) ;
		}

		fseek( fil, 0, SEEK_END);
		dataSize = ftell( fil);
		rewind( fil);

		if(dataSize == 0)
		{
			printf("error decodeFastq: no data in %s\n",argv[j]);
			fclose(fil);
			return(NULL);
		}

		if((data = (char *) malloc( (dataSize+1) * sizeof(char) )) == NULL)
		{
			printf("error decodeFastq: malloc error\n");
			fclose(fil);
			return(NULL);
		}
	
		for(i=0;i<dataSize;i+=1024*1000000){
			long size = 1024*1000000;
			if(dataSize-i < size)
				size = dataSize-i;
			//printf("1024*1000000 chunck %d\n",i);
			if(fread(data+i,sizeof(char), size, fil) != size)
			{
				printf("error decodeFastq: %s is too small\n",argv[j]);
				fclose(fil);
				return (NULL) ;
			}
		}
		fclose(fil);

		data[dataSize] = 0;
		nSeg = 0;
		dp = data;
		ret = 1;
		for(i=dataSize;--i>=0L;dp++)
		{
			if(*dp == '\n')
			{
				ret = 1;
				continue;
			}
			if(ret && (*dp == '@' || *dp == '>'))
			{
				nSeg++;
			}
			ret = 0;
		}

		if(nSeg == 0)
		{
			printf("error decodeFastq: no segments found\n");
			fclose(fil);
			return(NULL);
		}

		nseg += nSeg;
		if((seg = (fast *) realloc(seg, nseg*sizeof(fast))) == NULL)
		{
			printf("error decodeFastq: realloc error\n");
			fclose(fil);
			return(NULL);
		}

		cp = dp = data;
		segp = seg+nseg-nSeg;
		length = 0;
		ret = 1;
		segp->data = dp;
		for(i=dataSize;--i>=0L;dp++,length++)
		{
			if(*dp == '\n')
			{
				ret = 1;
				continue;
			}
			if(ret && (*dp == '@' || *dp == '>'))
			{
				if(length != 0)segp[-1].length = length;
				segp->data = dp;
				segp++;
				length = 0;
			}
			ret = 0;
		}
		segp[-1].length = length;
	}
	*nContigs = nseg;
	return(seg);
}

fasta *splitFastq ( fasta *iseg, int inSeg, int tB, int *nReads, int length, int step)
{
	char *data = NULL, *dp, *cp;
	long dataSize;
	int i,ret,nSeg=0,maxNseg;
	fasta *seg,*segp;
	int rc;

	maxNseg = tB/step + 1;
	if((seg = (fasta *) calloc(maxNseg, sizeof(fasta))) == NULL)
	{
		printf("error splitFastq: calloc error\n");
		return(NULL);
	}

	segp = seg;
	for(rc=0;rc<inSeg;rc++)
	{
		int size = length + strlen(iseg[rc].name)+strlen(".000000") + strlen("bases 1000000000 to 3000000000");
		if(!iseg[rc].finished) size += length;
		if(iseg[rc].length < length*2) continue;
		for(i=0;i<iseg[rc].length/step;i++)
		{
			int st, en;
			st = i*step;
			en = st+length;
			if(st >= iseg[rc].length) break;
			if(en >= iseg[rc].length)
			{
				st = iseg[rc].length-length-1;
				en = st+length;
			}
			if((segp->name = (char *) calloc(size, sizeof(char))) == NULL)
			{
				printf("error splitFastq: calloc error\n");
				return(NULL);
			}
			//sprintf(segp->name,"%s.%06d",iseg[rc].name,i);
			sprintf(segp->name,"%s",iseg[rc].name);
			segp->name2 = segp->name+strlen(segp->name)+1;
			sprintf(segp->name2,"bases %d to %d",st+1,en+1);
			segp->length = length;
			segp->data = segp->name2+strlen(segp->name2)+1;
			memcpy(segp->data,iseg[rc].data+st,length);
			segp->qual = segp->data + length;
			segp->finished = iseg[rc].finished;
			if(!iseg[rc].finished) memcpy(segp->qual,iseg[rc].qual+st,length);
			else segp->qual = NULL;
			segp++;
			nSeg++;
			if(nSeg >= maxNseg) break;
		}
		if(nSeg >= maxNseg) break;
	}
	*nReads = nSeg;
	if(nSeg == 0) return NULL;
	return(seg);
}




void fastaLC (fasta *seg, int nSeg)
{
	int i,j;
	fasta *segp;
	segp = seg;
	for(i=0;i<nSeg;i++,segp++)
	{
		char *b;
		b = (char*)segp->data;
		for(j=segp->length;--j>=0;b++)
			*b = tolower(*b);
	}
}

void fastaUC (fasta *seg, int nSeg)
{
	int i,j;
	fasta *segp;
	segp = seg;
	for(i=0;i<nSeg;i++,segp++)
	{
		char *b;
		b = (char*)segp->data;
		for(j=segp->length;--j>=0;b++)
			*b = toupper(*b);
	}
}

void error(char *out_message)
{
     printf("%s",out_message);
}

int countQthresh (fasta *seg, int qthresh)
{
	int i,j;
	fasta *segp;
	char *b;
	int n = 0;
	segp = seg;
	if(segp->finished) return(segp->length);
	b = (char*)segp->qual;
	for(j=segp->length;--j>=0;b++)
		if(*b >= qthresh)
			n++;
	return(n);
}



int reverseCompliment(fasta *seq, fasta *rseq)
{
	int i;
	int size=((seq->data != 0)+(seq->qual != 0))*seq->length;
	char *tp,*dp;

	rseq->length = 0;
	if(size <= 0 || seq->data == NULL)
		return(1);

	//printf("reverseCompliment: malloc %d",(size+10*strlen(seq->name)+100) * sizeof(char));
	//fflush(stdout);
	if((rseq->data = (char *)malloc((size+10*strlen(seq->name)+40) * sizeof(char))) == NULL)
	{
		printf("ERROR reverseCompliment: malloc\n");
		return(1);
	}
	//printf("  DONE\n");
	//fflush(stdout);

	rseq->name = rseq->data+size;
	strcpy(rseq->name,seq->name);
	rseq->path = seq->path;
	rseq->length = seq->length;
	rseq->finished = seq->finished;
	rseq->qual = NULL;
	if(seq->qual != NULL)
		rseq->qual = rseq->data+rseq->length;
	
	dp = rseq->data;
	tp = seq->data+seq->length;
	for(i=seq->length;--i>=0;)
	{
		int tmp = *--tp;
		if     (tmp == 'T') *dp++ = 'A';
		else if(tmp == 'G') *dp++ = 'C';
		else if(tmp == 'C') *dp++ = 'G';
		else if(tmp == 'A') *dp++ = 'T';
		else                *dp++ = tmp;
	}
	if(seq->qual != 0)
	{
		dp = rseq->qual;
		tp = seq->qual+seq->length;
		for(i=seq->length;--i>=0;)
			*dp++ = *--tp;
	}
	return(0);
}

void freeRead(fasta *seq)
{
	if(seq->data)
		free(seq->data);
	seq->data = NULL;
}

int duplicateRead(fasta *seq, fasta *dseq)
{
	int i;
	int size=((seq->data != 0)+(!seq->finished))*seq->length;
	char *tp,*dp;

	dseq->length = 0;
	dseq->data = NULL;
	if(size <= 0 || seq->data == NULL)
		return(1);

	//printf("duplicateRead: malloc %d",(size+10*strlen(seq->name)+100) * sizeof(char));
	//fflush(stdout);
	if((dseq->data = (char *)malloc((size + 10*strlen(seq->name) + 40 )* sizeof(char))) == NULL)
	{
		printf("ERROR reverseCompliment: malloc\n");
		return(1);
	}
	//printf("  DONE\n");
	//fflush(stdout);

	dseq->name = dseq->data+size;
	strcpy(dseq->name,seq->name);
	dseq->path = seq->path;
	dseq->SCFname = NULL;
	if(seq->SCFname != NULL)
	{
		dseq->SCFname = dseq->data+size+strlen(seq->name) + 1;
		strcpy(dseq->SCFname,seq->SCFname);
	}
	dseq->length = seq->length;
	dseq->finished = seq->finished;
	dseq->qual = NULL;
	if(!seq->finished)
		dseq->qual = dseq->data+dseq->length;

	dp = dseq->data;
	tp = seq->data;
	for(i=seq->length;--i>=0;)
		*dp++ = *tp++;

	if(!seq->finished && seq->qual != NULL)
	{
		dp = dseq->qual;
		tp = seq->qual;
		for(i=seq->length;--i>=0;)
			*dp++ = *tp++;
	}
	return(0);
}

int extractFastq(char *fname, char *pdata, long Size_pdata) {
  FILE *fil; 
  char *dp;
  long dataSize;
  long i;
  int nSeg;
  int state;

  if((fil = fopen(fname,"r")) == NULL) error("error: file not found\n");
  fseek( fil, 0, SEEK_END);
  dataSize = ftell( fil);
  rewind( fil);
  if(dataSize == 0) {
    fclose(fil);
    error("error: no data in file\n");
  }
  memset(pdata,'\0',Size_pdata);
  for(i=0;i<dataSize;i+=1024*1000000) {
    long int size = 1024*1000000;
    if(dataSize-i < size) size = dataSize-i;
    if(fread(pdata+i,sizeof(char), size, fil) != size) {
      fclose(fil);
      error("error: file is too small\n");
    }
  }
  fclose(fil);
  pdata[dataSize] = 0;
  nSeg = 0;
  dp = pdata;
  if(*dp != '>' && *dp != '@') fprintf(stderr,"Corrupt input file %s *dp=>%c< [1].\n", fname, *dp);

  /****************************************************************************
   * Use simple 4-state automaton.
   * -----------------------------
   * state 0 - character is part of sequence header.
   * state 1 - character is part of sequence bases.
   * state 2 - character is part of quality header.
   * state 3 - character is part of quality values.
   ***************************************************************************/
  state=0;
  for(i=dataSize;--i>=0L;dp++) {
    if(state == 0) {
      if(*dp == '@' || *dp == '>') {
        nSeg++;
        continue;
      }
      else if(*dp == '\n') {
        if(i == 0L) error("Corrupt input file [2].\n");
        state = 1;
        continue;
      }
      else continue;
    }
    else if(state == 1) {
      if((*dp == '\n')) {
        if(i == 0L) break;
        else if((*(dp+1) == '@') || (*(dp+1) == '>')) {
          state=0;
          continue;
        }
        else if(*(dp+1) == '+') {
          state=2;
          continue;
        }
        else continue;
      }
      else continue;
    }
    else if(state == 2) {
      if((*dp == '\n')) {
        if(i == 0L) error("Corrupt input file [3].\n");
        state = 3;
        continue;
      }
      else continue;
    }
    else if(state == 3) {
      if((*dp == '\n')) {
        if(i == 0L) break;
        else if(*(dp+1) == '@' || *(dp+1) == '>') {
          state=0;
          continue;
        }
        else continue;
      }
      else continue;

    }
    else error("Corrupt input file [4].\n");
  }
  if(nSeg == 0) error("error: no segments found\n");
  return(nSeg);
}




fasta *decodeFastq(char *fname,int *nContigs,long *tB,char* pdata,long Size_pdata,fasta *segg)
{
  FILE *fil;
  char *data = NULL, *dp, *cp;
  long dataSize, totalbases=0;
  long i;
  int ret;
  int nSeg;
  fasta *seg,*segp; 
  int state, isfastq = 0;
        
  if((fil = fopen(fname,"r")) == NULL) error("error: file not found\n");
  fseek( fil, 0, SEEK_END);
  dataSize = ftell( fil);
  rewind( fil); 
  if(dataSize == 0) {
    fclose(fil);
    error("error decodeFastq: no data in file\n");
  } 
  memset(pdata,'\0',Size_pdata);
  for(i=0;i<dataSize;i+=1024*1000000) {
    long int size = 1024*1000000;
    if(dataSize-i < size) size = dataSize-i; 
    if(fread(pdata+i,sizeof(char), size, fil) != size) {
      fclose(fil);
      printf("%s\n",fname);
      error("error: file is too small\n");
    }     
  }       
  fclose(fil);
  pdata[dataSize] = 0;
  nSeg = 0;
  dp = pdata;
    
  if     (*dp == '@') isfastq = 0;
  else if(*dp == '>') isfastq = 1;
  else                error("Corrupt input file [5].\n");
        
  /****************************************************************************
   * Use simple 4-state automaton.
   * -----------------------------
   * state 0 - character is part of sequence header.
   * state 1 - character is part of sequence bases.
   * state 2 - character is part of quality header.
   * state 3 - character is part of quality values.
   *
   * This should be used for reading data too - current code is unstable.
   *
   ***************************************************************************/
  state=0;
  for(i=dataSize;--i>=0L;dp++) {
    if(state == 0) {
      if(*dp == '@' || *dp == '>') {
        nSeg++;
        continue;
      }
      else if(*dp == '\n') {
        if(i == 0L) error("Corrupt input file [6].\n");
        state = 1;
        continue;
      }
      else continue;
    }
    else if(state == 1) {
      if((*dp == '\n')) {
        if(i == 0L) break;
        else if((*(dp+1) == '@') || (*(dp+1) == '>')) {
          state=0;
          continue;
        }
        else if(*(dp+1) == '+') {
          state=2;
          continue;
        }
        else continue;
      }
      else continue;
    }
    else if(state == 2) {
      if((*dp == '\n')) {
        if(i == 0L) error("Corrupt input file [7].\n");
        state = 3;
        continue;
      }
      else continue;
    }
    else if(state == 3) {
      if((*dp == '\n')) {
        if(i == 0L) break;
        else if(*(dp+1) == '@' || *(dp+1) == '>') {
          state=0;
          continue;
        }
        else continue;
      }
      else continue;

    }
    else error("Corrupt input file [8].\n");
  }
  if(nSeg == 0) error("no segments found\n");
  cp = dp = pdata;
  seg = segg;
  segp = seg;
  while(dp < pdata+dataSize) {
    if(*dp != '@' && *dp != '>') {
      dp++;
      continue;
    }
    else if ( dp > pdata && dp[-1] != '\n' ) {
      dp++;
      continue;
    }
    else {
      int tmp;
      int k;
      dp++;
      segp->name = (char *)(long int)(cp-pdata);
      while((*cp++ = tmp = *dp++) != ' ' && tmp != '\t' && tmp != '\n');
      cp[-1] = 0;
      if(tmp != '\n') {
        segp->name2 = (char *)(long int)(cp-pdata);
        while((*cp++ = *dp++) != '\n');
        cp[-1] = 0;
      }
      else segp->name2 = NULL;
      segp->data = (char *)(long int)(cp-pdata);
      for(k=0;(tmp = *dp) &&
              !((tmp == '+' || tmp == '>') &&
              dp[-1] == '\n');dp++) {
        if(tmp == '\n') continue;
        *cp++ = tmp;
        k++;
      }
      segp->length = k;
      segp->finished = 1;
      segp->qual = NULL;
      if(isfastq) {
        segp->finished = 1;
        segp->qual = NULL;
      }
      else {
        segp->finished = 0;
        k = 0;
        if(tmp == '+') {
          int cp_set = 0;
          int comp = 1;
          while(*dp++ != '\n');
          for(i=0;i<20;i++) {
            if(dp[i] == 0) break;
            if(dp[i] == ' ') {comp = 0; break;}
          }
          segp->qual = (char *)(long int)(cp-pdata);
          *cp = 0;
          if(comp) for(k=0;
                       (tmp = *dp) && (k < segp->length);
                       dp++) {
             if(tmp == '\n') continue;
             *cp++ = tmp-041;
             k++;
          }
          else for(k=0;
                   (tmp = *dp) && (k < segp->length);
                   dp++) {
            if(tmp >= '0' && tmp <= '9') {
              if(*cp > 9) continue;
              *cp *= 10;
              *cp += tmp-'0';
              cp_set = 1;
              continue;
            }
            if(cp_set) {
              cp++;
              *cp = 0;
              cp_set = 0;
              k++;
            }
          }
        }
        if((segp-seg + 1 < nSeg) && !(*(dp+1) == '@' && *dp == '\n')&&(segp->length!=k)) 
	{
          fprintf(stderr, "Warning: %-20.20s has %d bases and %d quals, skipping entry\n",pdata+((long int)segp->name),segp->length,k);
          cp = pdata+((long int)segp->name);
          continue;
        }
      }
      segp++;
      if(segp-seg > nSeg) break;
     }
   }
   nSeg = segp-seg;
   for(i=nSeg,segp=seg;--i>=0;segp++) {
     segp->name = pdata+((long int)segp->name);
     if(segp->name2) segp->name2 = pdata+((long int)segp->name2);
     segp->data = pdata+((long int)segp->data);
     segp->qual = pdata+((long int)segp->qual);
     if(segp->finished) segp->qual = NULL;
     totalbases+=segp->length;
   }
   *nContigs = nSeg;
   *tB = totalbases;
   return(seg);
}

