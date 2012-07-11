/**************************************************************************
 *
 *-------------------------------------------------------------------------
 * Project Notes : generates random shotgun reads from a fasta input file
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
 *	Copyright (C) 2000 by Genome Research Limited, All rights reserved.
 *
 **************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "fasta.h"
fasta *srandFastq ( fasta *iseg, int inSeg, int tB, int *nReads, int length, int step, int insert1, float percent1,
					int insert2, float percent2, int insert3, float percent3);

int main(int argc, char **argv)
{
	fasta *exp,*exptmp; 
	int i,j,ci;
	long totalBases;
	long s;
	int nReads;
	int length, step, insert1, insert2, insert3;
	float percent1, percent2, percent3, total;
	FILE *outf,*outq;
	char qual[100];
	double p[100];
	char b[4] = {'a','c','g','t'};

	p[0] = 1.;
	for(i=1;i<100;i++)
		p[i] = pow(10.,-i/10.);

	exp = decodeFastq ( argv[1], &nReads, &totalBases, 23);
	if(exp == NULL)
	{
		printf("ERROR pileup: no data found\n");
		exit(1);
	}
	outf = fopen(argv[2],"w");
	strcpy(qual,argv[2]);
	strcat(qual,".qual");
	outq = fopen(qual,"w");
	for(ci=0;ci<nReads;ci++) 
	{
		int sReads;
		int si;
		putc('>',outf);
		fprintf(outf,"%s",exp[ci].name);
		if(exp[ci].name2)
			fprintf(outf," %s", exp[ci].name2);
		putc('\n',outf);
		for(i=0;i<exp[ci].length;i++)
		{
			putc(exp[ci].data[i],outf);
			if(i%60 == 59 && i+1 != exp[ci].length)
				putc('\n',outf);
		}
		putc('\n',outf);
		putc('>',outq);
		fprintf(outq,"%s",exp[ci].name);
		if(exp[ci].name2)
			fprintf(outq," %s",exp[ci].name2);
		putc('\n',outq);
		for(i=0;i<exp[ci].length;i++)
		{
			fprintf(outq," %d",exp[ci].qual[i]);
			if(i%30 == 29 && i+1 != exp[ci].length)
				putc('\n',outq);
		}
		putc('\n',outq);
	}
}
