/**************************************************************************
 *
 *-------------------------------------------------------------------------
 * Project Notes : Assembles anchored reads onto genetic data
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
 *	Copyright (C) 1998-2001 by Genome Research Limited, All rights reserved.
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
#include "ssaha.h"

static int KLEN=21;

int main(int argc, char **argv)
{
	int i;
    int SEG_LEN=KLEN;
//	__delayed_free = 0;//so we can free memory right away

	printf("SsahaHist Version 1.02\n");
	if(argc < 2)
	{
		printf("Usage: %s <-kmer 21> <-depth 300> <-lelngth 2000> <-mismatch 4> <-reads 20000> <mates file> <input_fastq_file> <output_contig_file> \n",argv[0]);
		printf("Note: for all default values, the use of any parameters is optional.\n");
		exit(1);
	}

        for(i=1;i<argc;i++)
        {
           if(!strcmp(argv[i],"-kmer"))
           {
             sscanf(argv[++i],"%d",&SEG_LEN);
           }
        } 
	ssaha_init(argv, argc, SEG_LEN);
        return EXIT_SUCCESS;
}
