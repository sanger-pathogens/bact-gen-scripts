/***********************************************************************\
 *                                                                     * 
 *                      PROJECT       FuzzyPath                        *
 *                                                                     * 
 *---------------------------------------------------------------------*
 *                                                                     * 
 *                  Path of Euler to Extend Reads of Solexa            *
 *                                                                     *
 *                                By                                   *
 *                                                                     *
 *                            Zemin Ning                               *
 *                                                                     *
 *          Copyright (C) 2006-2008 by Genome Research Limited         *
 *                         All rights reserved                         *
 *                                                                     *
 *---------------------------------------------------------------------*
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
 *---------------------------------------------------------------------*/


/****************************************************************************/

 
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

static int *list,*readlength,*out_list,*ctg_list,*rd_group;
static int *num_pair,*mat_pair,*ins_pair,*dev_pair,*rpf_contig,*rpf_offset;
static int *kmer_nhits,*kmer_idgrp,*kmer_locus,*kmer_masks,*kmer_index,*kmer_links,*kmer_nodes;
static int *node_head,*node_list,*node_index,*kmer2node,*node_word,*node_next,num_nodes;
static int *node_rept,node_hits;
static int *repeat_word,*repeat_next; 
static int *path_kmer,*path_next,*path_idir,*path_cover;
static int *path2_kmer,*path2_next,*path2_idir,*neib_qual;
static int *n_list,*read2kmer,*reads_mask,*reads_post,*reads_mast,*reads_rept,*reads_maps,**m_qualy;
static char **m_align,**cons_reads,**exts_reads,**exts_reads2,**neib_kmer;
static char *neib1_kmer,*kmer_baseF,*kmer_baseR;
static int *kmer2_hits;
static char *kmer_cons,*rbase;
static int *kmer_qual,*m_pileup;
static signed char *patch_kmer2;
static long hmask;
static int *sm,*sm_head,*sm_list;
static int kmer_depth=1500;
static int kmer_depth2=2;
static int read_length = 3000;
static int set_offset = 4;
static int num_outs = 0;
static int rp_insert = 0;
static int num_overlap = 9;
static int reads_used = 0;
static long baseSize;
static char *path_cons,*path2_cons;
static int n_repeats = 0;
static int hit_single = 3;
static int i_contig;
static int e_contig;
static int move_dir;
static long num_akmer;
static int num_cover,num_ukmer,num_kmers;
static int kmer_dir = 0;
static int kmer2_dir = 0;
static int n_kmer2 = 2;
static int i_kmer2 = 2;
static int rcdex_pres = 0;
static int path_dir = 0;
static int num_cons;
static int i_read = 0;
static int weird_kmer = 0;
static int hit_ikk = 0;
static int sense_flag = 1;
static int a454_flag = 1;
static int abi_flag = 0;
static int i_seq = 0;
static int set_mismatch = 5;
static int rshift = 6,nshift;
static int print_flag = 0;
static int print_ctg = 0;
static int nkm_len = 5;
static int debug_flag = 0;
static int n_node = 0;
static float cfactor = 1.5;
static float set_rate = 0.59;
static int k_node = 3;
static int check_fuzzy = 0;
static int junction1=0,junction2=0;
static int qthresh=25;
static int qthresh2=15;

static void hpsortl(long n, unsigned long *ra)
{
	long l,j,ir,i;
	unsigned long rra;

	if(n < 2)
          return;
        l=(n >> 1)+1;
	ir=n;
	for (;;) {
		if (l > 1) {
			rra=ra[--l];
		} else {
			rra=ra[ir];
			ra[ir]=ra[1];
			if (--ir == 1) {
				ra[1]=rra;
				return;
			}
		}
		i=l;
		j=l << 1;
		while (j <= ir)	{
			if (j < ir && ra[j] < ra[j+1]) ++j;
			if (rra < ra[j]) {
				ra[i]=ra[j];
				j += (i=j);
			}
			else j=ir+1;
		}
		ra[i]=rra;
	}
}


/*   Subroutine to sort the DNA sequences into a matrix table   */
/* =============================================  */
void Data_Input(char **argv, int argc, int SEG_LEN)
/* =============================================  */
{
     long i,j,k,n_Sbase=SEG_LEN;
     long ps=0,ns=0,IntSeg=0,IntBase=0,lii;
     fasta *seq;
     fast  *seq2;
     FILE *fpMate;
     int nSeq;
     long totalBases;
     long mhistc = 0;
     long mhistcc = 0;
     int nsorts = 1024;
//     int tb = 0;
     int args,ac;
     int nseq = 0,seqc;
     int sshift;
     unsigned long nmask;
     int num_sect; //id_read will never exceed kmer_depth
     char line[500] = {0},outName[60]={0},outFast[60]={0},syscmd[200]={0};
     char *ptr,base[20],zero[20]={0},line2[500]={0};
     int rd,read_pair[200];
     long kmask = (1L<<(n_Sbase<<1))-1;
     double Qerr[100];
     long gtBases=0,nclip=0;
     Qerr[0] = 1.0;
     void Readpair_Stage(int brr,char **argv, int argc, int args);
     void Kmer_Process(int brr,int len,int crr,char **argv, int argc, int args);

     for(i=1;i<100;i++) Qerr[i] = pow((double)10.,(double)-.1*i);

/*   sort all the names of genes or name entries   */
     printf("Input data starts ...\n");
     args=1;
     for(i=1;i<argc;i++)
     {
       if(!strcmp(argv[i],"-kmer"))
       {
	 i++;
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-depth"))
       {
         sscanf(argv[++i],"%d",&kmer_depth);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-depth2"))
       {
         sscanf(argv[++i],"%d",&kmer_depth2);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-print"))
       {
         sscanf(argv[++i],"%d",&print_flag);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-qual"))
       {
         sscanf(argv[++i],"%d",&qthresh);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-contig"))
       {
         sscanf(argv[++i],"%d",&print_ctg);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-single"))
       {
         sscanf(argv[++i],"%d",&hit_single);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-454"))
       {
         sscanf(argv[++i],"%d",&a454_flag);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-abi"))
       {
         sscanf(argv[++i],"%d",&abi_flag);
	 a454_flag = 0;
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-mismatch"))
       {
         sscanf(argv[++i],"%d",&set_mismatch);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-shift"))
       {
         sscanf(argv[++i],"%d",&rshift);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-cover"))
       {
         sscanf(argv[++i],"%f",&cfactor);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-iread"))
       {
         sscanf(argv[++i],"%d",&i_read);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-insert"))
       {
         sscanf(argv[++i],"%d",&rp_insert);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-offset"))
       {
         sscanf(argv[++i],"%d",&set_offset);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-sense"))
       {
         sscanf(argv[++i],"%d",&sense_flag);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-node"))
       {
         sscanf(argv[++i],"%d",&k_node);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-overlap"))
       {
         sscanf(argv[++i],"%d",&num_overlap);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-rate"))
       {
         sscanf(argv[++i],"%f",&set_rate);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-fuzzy"))
       {
         sscanf(argv[++i],"%d",&check_fuzzy);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-debug"))
       {
         sscanf(argv[++i],"%d",&debug_flag);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-length"))
       {
         sscanf(argv[++i],"%d",&read_length);
         args=args+2;
       } 
       else if(!strcmp(argv[i],"-reads"))
       {
         sscanf(argv[++i],"%d",&num_outs);
         args=args+2;
       } 
     }
     for(ac=0;ac<argc;ac++)
	printf("%s ",argv[ac]);
     printf("\n");
     printf("kmer size:         %d\n",n_Sbase);
     printf("depth:             %d\n",kmer_depth);
     printf("ph2 depth:         %d\n",kmer_depth2);

/*     if((fpMate = fopen(argv[args],"r")) == NULL)
     {
       printf("Error fmate: mate input\n");
       exit(1);
     }
     nseq = 0;
     while(!feof(fpMate))
     {
       fgets(line,500,fpMate);
       if(feof(fpMate)) break;
       nseq++;
     }
     fclose(fpMate);

     Readpair_Stage(nseq,argv,argc,args);       */

     printf("mates done\n");
     Kmer_Process(nseq,SEG_LEN,kmer_depth,argv,argc,args);


     printf("Relation matrix finished\n");

     printf("All jobs finished\n");
     return;
}


int ssaha_init( char **argv, int argc, int SEG_LEN)
{

    Data_Input(argv,argc,SEG_LEN);
    return(1);   
}

/*   Subroutine to sort the DNA sequences into a matrix table   */
/* =============================================  */
void Kmer_Process(int nRead,int SEG_LEN,int n_depth,char **argv, int argc, int args)
/* =============================================  */
{
     long i,j,k,w,m,iseq,n_Sbase=SEG_LEN;
     long ps=0,ns=0,IntSeg=0,IntBase=0,lii;
     char *b,*s;
     fasta *seqp;
     fasta *seq,*segg;
     fasta  *seq2;
     int nSeq;
     FILE *fp;
     long totalBases;
     long mhistc = 0;
     long mhistcc = 0;
     int nsorts = 1024;
     int ac,max_nhit,kmer_last[3];
     int seqc,num_sect;
     unsigned long seqcharc = 0;
     int step_flag,NLINK,max_contig,sshift,kshift=2,wshift=31;
     unsigned long nmask;
     long *patch_array,*patch_head,*ray,*patch_nlink,*patch_plink;
     int *mm,*patch_index,*patch_list,*dex,*patch_ofset,*patch_2kmer;
     int **imatrix(long nrl,long nrh,long ncl,long nch);
     char **cmatrix(long nrl,long nrh,long ncl,long nch);
     unsigned int **uimatrix(long nrl,long nrh,long ncl,long nch);
     long **limatrix(long nrl,long nrh,long ncl,long nch);
     void Euler_Path(int iseq,int nRead,int n_depth,int n_Sbase,fasta *seq,int *patch_ofset,int *patch_index,int *patch_2kmer,long *patch_array,long *patch_head,long *patch_nlink,long *patch_plink);
     void Euler_Path2(int *break_kmer,int b_kmer,int joint,int nRead,int n_depth,int n_Sbase,fasta *seq,int *patch_ofset,int *patch_index,int *patch_2kmer,long *patch_array,long *patch_head,long *patch_nlink,long *patch_plink);
     void ArraySort_Mix(int n, long *arr, int *brr);
     char syscmd[200]={0};
     long kmask = (1L<<(n_Sbase<<1))-1;
     double Qerr[100];
     FILE *namef;
     long gtBases=0,nclip=0,next_link,last_link;
     Qerr[0] = 1.0;
     void ArraySort_Int2(int n, int *arr, int *brr);
     int num_seqque;
     long Size_q_pdata,mmask[100];
     char *pdata;
     void Mmask_array(long *mmask,int n_Sbase);
     void Mismatch_Kmer(long kmer1, long kmer2,int n_Sbase,long *mmask,int *num_mismatch,int *mismatch_loci);
     void Short_Path(int idp,int id_node,int kmer_len,int *sm,int *sm_list,int *short_flag,int *walk_len,int *kmer_last);

     kmer_depth=n_depth;
     sshift = (64-(n_Sbase<<1));
     nmask = (1L<<rshift)-1;

     nRead = 0;

     if((fp=fopen(argv[args+1],"rb"))==NULL) error("Cannot open file\n");
     fseek(fp, 0, SEEK_END);
     Size_q_pdata = ftell(fp) + 1;
     fclose(fp);
     if((pdata=(char*)calloc(Size_q_pdata,sizeof(char)))==NULL)
       error("calloc pdata\n");
     num_seqque = extractFastq(argv[args+1],pdata,Size_q_pdata);
     if((segg=(fasta*)calloc((num_seqque+1),sizeof(fasta)))==NULL)
       error("calloc segg\n");
     if((seq=decodeFastq(argv[args+1],&nSeq,&totalBases,pdata,Size_q_pdata,segg))==NULL)
       error("no query data found.\n");  
     if(seq == NULL)
     {
     	printf("ERROR pileup: no data found\n");
     	exit(1);
     }
     nRead += nSeq;

/*   temp mates handling        */
     if((num_pair= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - num_pair\n");
       exit(1);
     }
     if((mat_pair= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - mat_pair\n");
       exit(1);
     }
     if(rp_insert > 0)
     {
          for(i=0;i<nSeq;i++)
	  {
	    if((i%2) == 0)
	      mat_pair[i] = i+1;
	    else
	      mat_pair[i] = i-1;
	    num_pair[i] = 1;
	  }
     }
/*     if((ins_pair= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ins_pair\n");
       exit(1);
     }
     if((dev_pair= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - dev_pair\n");
       exit(1);
     }  */
     if((rpf_contig= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - rpf_contig\n");
       exit(1);
     }
     if((rpf_offset= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - rpf_offset\n");
       exit(1);
     }


     mhistc = 0;
     gtBases = 0;
     nsorts = nRead;
     for(i=0;nsorts > 10;i++)
     {
        nsorts = nsorts>>1;
     } 
     nsorts = 1L<<i;
     nshift = (n_Sbase<<1)-i;//2^i=nsorts
     hmask = (1L<<nshift)-1;
     if((patch_list= (int *)calloc(nsorts,sizeof(int))) == NULL)
     {
       printf("phusion: calloc - patch_list\n");
       exit(1);
     }
     if((patch_head= (long *)calloc(nsorts,sizeof(long))) == NULL)
     {
       printf("phusion: calloc - patch_head\n");
       exit(1);
     }


     fastaUC(seq,nSeq);
     for(iseq=0;iseq<nSeq;iseq++)
     {
        double e=0.0,eh[32];
	long IntBaseRC;
        long IntSegRC=0;
	int pos=0;
        char *q;

        IntSeg=0;
        seqp=seq+iseq;
        b = (char*)seqp->data;
        if(seqp->finished) {
          q = NULL;
        } else {
          q = NULL; //(char*)seqp->qual;
        }

        ns=(seqp->length);
        seqcharc += strlen(seqp->name);
        k=n_Sbase;
	if(k > ns) k=ns;

        for(k=n_Sbase;k>0;k--,b++,pos++)
        {
           if     (*b == 'A') IntBase=0;
           else if(*b == 'C') IntBase=1;
           else if(*b == 'G') IntBase=2;
           else if(*b == 'T') IntBase=3;
           else
           {
	     e += eh[pos] = 1.;
	     continue;
           }
           e += eh[pos] = (q != NULL) ? Qerr[q[pos]] : 0;
	   IntBaseRC=(IntBase^3)<<((n_Sbase-k)<<1);
           IntBase=IntBase<<((k-1)<<1);
           IntSeg+=IntBase;
	   IntSegRC+=IntBaseRC;
        }
        for(j=ns-n_Sbase;--j>=0;b++,pos++)
        {
	   int p=pos%n_Sbase;
	   if(e < .99) 
           {
	      if(IntSeg > IntSegRC) 
              {
                  int itt = IntSegRC>>nshift;
                  patch_list[itt]++;
	      } 
              else 
              {
                  int itt = IntSeg>>nshift;
                  patch_list[itt]++;
	      }
	      gtBases++;
	   }
	   e -= eh[p];
	   IntSeg = (IntSeg<<2)&kmask;
	   IntSegRC = (IntSegRC>>2)&kmask;
           if     (*b == 'A') IntBase=0;
           else if(*b == 'C') IntBase=1;
           else if(*b == 'G') IntBase=2;
           else if(*b == 'T') IntBase=3;
           else
           {
	     e += eh[p] = 1.;
	     continue;
           }
	   e += eh[p] = (q != NULL) ? Qerr[q[pos]] : 0; 
	   IntBaseRC=(IntBase^3)<<((n_Sbase-1)<<1);
           IntSeg+=IntBase;
	   IntSegRC+=IntBaseRC;
        }

	if(e < .99) 
        {
	  if(IntSeg > IntSegRC) 
          {
              int itt = IntSegRC>>nshift;
              patch_list[itt]++;
	  } 
          else 
          {
              int itt = IntSeg>>nshift;
              patch_list[itt]++;
	  }
	  gtBases++;
        }
     }

     if((read2kmer= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - read2kmer\n");
       exit(1);
     }
     if((reads_mask= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - reads_mask\n");
       exit(1);
     }
     if((reads_post= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - reads_post\n");
       exit(1);
     }
     if((reads_mast= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - reads_mast\n");
       exit(1);
     }
     if((reads_rept= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - reads_rept\n");
       exit(1);
     }
     if((reads_maps= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - reads_rept\n");
       exit(1);
     }
     printf("data reading finished ...\n");


     printf("histgram finished ...\n");
     printf("Input data finished one set, nsorts=%d, totalgoodBP=%ld\n",nsorts,gtBases);
     max_nhit=0;
     patch_head[0]=0;
     for(i=0;i<nsorts;i++) 
     {
        if(patch_list[i]>max_nhit)
          max_nhit = patch_list[i];
        if(i>0)
          patch_head[i] = patch_head[i-1] + patch_list[i-1];
//        patch_list[i] = 0;
     }
     memset(patch_list,0,4*nsorts);  
     mhistcc=patch_head[nsorts-1]+patch_list[nsorts-1]+20*nsorts;

     printf("setting array memory: %ld\n",mhistcc);
     if((patch_array= (long *)calloc(mhistcc,sizeof(long))) == NULL)
     {
       printf("ssaha: calloc - patch_array\n");
       exit(1);
     }
     if((patch_nlink= (long *)calloc(mhistcc,sizeof(long))) == NULL)
     {
       printf("ssaha: calloc - patch_nlink\n");
       exit(1);
     }
     if((patch_plink= (long *)calloc(mhistcc,sizeof(long))) == NULL)
     {
       printf("ssaha: calloc - patch_plink\n");
       exit(1);
     }
     if((patch_index= (int *)calloc(mhistcc,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - patch_index\n");
       exit(1);
     }
     if((patch_2kmer= (int *)calloc(mhistcc,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - patch_2kmer\n");
       exit(1);
     }
     if((kmer_index= (int *)calloc(mhistcc,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - kmer_index\n");
       exit(1);
     }
     if((patch_ofset= (int *)calloc(mhistcc,sizeof(int))) == NULL)
     {     
       printf("ssaha: calloc - patch_ofset\n");
       exit(1);
     }
     if((patch_kmer2= (signed char *)calloc(mhistcc,sizeof(signed char))) == NULL)
     {     
       printf("ssaha: calloc - patch_kmer2\n");
       exit(1);
     }

     seqc = 0;

     mhistcc=0;


     for(iseq=0;iseq<nSeq;iseq++,seqc++)
     {
        char *q;
        double e=0,eh[32],e2=0.0;
        long IntBaseRC;
        long IntSegRC=0;
//        long IntSeg=0;
        int IntBaseRC2=0,IntSegRC2=0;
        int IntBase2=0,IntSeg2=0;
        int ic,pos=0;

        IntSeg=0;
        seqp=seq+iseq;
//        if(iseq>277836)
//          printf("base: %d\n",iseq);
        b = (char*)seqp->data; 
        if(seqp->finished) {
          q = NULL;
        } else {
          q = NULL; //(char*)seqp->qual;
        } 
        ns=(seqp->length);
        last_link = -1;
        next_link = -1;

        for(k=n_Sbase;k>0;k--,b++,pos++)
        {
           if     (*b == 'A') IntBase=0;
           else if(*b == 'C') IntBase=1;
           else if(*b == 'G') IntBase=2;
           else if(*b == 'T') IntBase=3;
           else
           {
	     e += eh[pos] = 1.;
	     continue;
           }
           e += eh[pos] = (q != NULL) ? Qerr[q[pos]] : 0;
	   IntBaseRC=(IntBase^3)<<((n_Sbase-k)<<1);
           IntBase=IntBase<<((k-1)<<1);
           IntSeg+=IntBase;
	   IntSegRC+=IntBaseRC;
//                printf("kkmer: %c %d %d %ld %ld %s %ld %ld\n",*b,pos,n_Sbase-k,IntBase,IntBaseRC,seqp->name,IntSeg,IntSegRC);
        }
        for(j=ns-n_Sbase;--j>=0;b++,pos++)
        {
	   int p=pos%n_Sbase;
           long patch_pos;

/*         mapping second kmer word    */
           IntSeg2 = -1;
           IntSegRC2 = -1;
           if(pos<=(ns-n_kmer2-1))
           {
             s = (char*)seqp->data + pos+1;
             e2 = 1.0;
             IntSeg2 = 0;
             IntSegRC2 = 0;
             for(w=n_kmer2;w>0;w--)
             {
                if     (*s == 'A') IntBase2=0;
                else if(*s == 'C') IntBase2=1;
                else if(*s == 'G') IntBase2=2;
                else if(*s == 'T') IntBase2=3;

                IntBaseRC2=(IntBase2^3)<<((n_kmer2-w)<<1);
                IntBase2=IntBase2<<((w-1)<<1);
                IntSeg2=IntSeg2+IntBase2;
                IntSegRC2=IntSegRC2+IntBaseRC2;
                if((*s != 'A')&&(*s != 'C')&&(*s != 'G')&&(*s != 'T'))
                  e2 = 0.0;
                s++;
             }
           }
	   if(e < .99) 
           {
              long itt;
	      if(IntSeg >= IntSegRC) 
              {
                itt = IntSegRC>>nshift;
                patch_pos = patch_head[itt]+patch_list[itt];
                patch_array[patch_pos]=IntSegRC;
                patch_index[patch_pos]=seqc;
                patch_ofset[patch_pos]=((ns-pos)<<kshift) + 2;
//                patch_kmer2[patch_pos] = IntSegRC2;
                if(next_link>=0)
                  patch_nlink[next_link] = (itt<<nshift) + patch_list[itt];
                patch_list[itt]++;
	        mhistcc++;
	      } 
              else 
              {
                itt = IntSeg>>nshift;
                patch_pos = patch_head[itt]+patch_list[itt];
                patch_array[patch_pos]=IntSeg;
                patch_index[patch_pos]=seqc;
//                patch_ofset[patch_head[itt]+patch_list[itt]]=((pos-n_Sbase)<<kshift) + 1;
                patch_ofset[patch_pos]=((pos-n_Sbase)<<kshift) + 1;
                if(next_link>=0)
                  patch_nlink[next_link] = (itt<<nshift) + patch_list[itt];
                patch_list[itt]++;
	        mhistcc++;
	      }
              if(IntSeg==IntSegRC)
              {
                patch_array[patch_pos]=-100;
		reads_rept[seqc] = 1;
                printf("kkmer: %d %d %s %ld %ld\n",j,pos,seqp->name,IntSeg,IntSegRC);
              }
                patch_kmer2[patch_pos] = IntSeg2;
              patch_plink[patch_pos] = last_link;
              last_link = (itt<<nshift) + patch_list[itt]-1;
              next_link = patch_pos;
	      gtBases++;
	   }
	   e -= eh[p];
	   IntSeg = (IntSeg<<2)&kmask;
	   IntSegRC = (IntSegRC>>2)&kmask;
           if     (*b == 'A') IntBase=0;
           else if(*b == 'C') IntBase=1;
           else if(*b == 'G') IntBase=2;
           else if(*b == 'T') IntBase=3;
           else
           {
	     e += eh[p] = 1.;
             continue;
           }
	   e += eh[p] = 0;
	   IntBaseRC=(IntBase^3)<<((n_Sbase-1)<<1);
           IntSeg+=IntBase;
	   IntSegRC+=IntBaseRC;
        }


        IntSeg2 = -1;
        IntSegRC2 = -1;
        if(e < .99) 
        {
          long patch_pos,itt;
	  if(IntSeg >= IntSegRC) 
          {
            itt = IntSegRC>>nshift;
            patch_pos = patch_head[itt]+patch_list[itt];
            patch_array[patch_pos]=IntSegRC;
            patch_index[patch_pos]=seqc;
            patch_ofset[patch_pos]=((ns - pos)<<kshift) + 2;
//            patch_nlink[next_link] = (itt<<nshift) + patch_list[itt];
            patch_nlink[next_link] = (itt<<nshift) + patch_list[itt];
            patch_list[itt]++;
	    mhistcc++;
	  } 
          else 
          {
            itt = IntSeg>>nshift;
            patch_pos = patch_head[itt]+patch_list[itt];
            patch_array[patch_pos]=IntSeg;
            patch_index[patch_pos]=seqc;
//            patch_ofset[patch_head[itt]+patch_list[itt]]=((pos-n_Sbase)<<kshift) + 1;
            patch_ofset[patch_pos]=((pos-n_Sbase)<<kshift) + 1;
            patch_nlink[next_link] = (itt<<nshift) + patch_list[itt];
            patch_list[itt]++;
	    mhistcc++;
	  }
          if(IntSeg==IntSegRC)
          {
            patch_array[patch_pos]=-100;
                printf("kkmer: %d %d %s %ld %ld\n",j,pos,seqp->name,IntSeg,IntSegRC);
	    reads_rept[seqc] = 1;
          }
          patch_plink[patch_pos] = last_link;
          last_link = (itt<<nshift) + patch_list[itt]-1;
          next_link = patch_pos;
        }
        patch_nlink[next_link] = -1;
     }

     printf("matrix memory setting ...\n");
     m_align=cmatrix(0,kmer_depth,0,n_depth);
     neib_kmer=cmatrix(0,n_depth,0,nkm_len+1);
     m_qualy=imatrix(0,kmer_depth,0,n_depth);
     if((neib1_kmer= (char *)calloc(n_depth,sizeof(char))) == NULL)
     {
       printf("ssaha: calloc - neib1_kmer\n");
       exit(1);
     }
     if((neib_qual= (int *)calloc(n_depth,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - neib_qual\n");
       exit(1);
     }
     if((kmer2_hits= (int *)calloc(n_depth,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - kmer2_hits\n");
       exit(1);
     }
     if((n_list= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - n_list\n");
       exit(1);
     }

     memset(read2kmer,-1,nRead*4);
     ray = patch_array;
     dex = kmer_index;
     max_nhit=0;
     num_kmers = 0;
     for(i=0;i<nsorts;i++) 
     {
        if(patch_list[i] > 1)
        {
          for(k=0;k<patch_list[i];k++)
             kmer_index[patch_head[i]+k]=k;
          for(k=0;k<patch_list[i];k++)
          {
             j = k+1;
             num_sect = 1;
             while((patch_index[patch_head[i]+k]==patch_index[patch_head[i]+j])&&(patch_array[patch_head[i]+j]==patch_array[patch_head[i]+k])&&(patch_index[patch_head[i]+j]==patch_index[patch_head[i]+k]))
             {
               int idk = patch_index[patch_head[i]+k];
               int idk2 = patch_index[patch_head[i]+j];
/*               int ref_offset = (int)(patch_ofset[patch_head[i]+k]>>kshift);
               int ref_offset2 = (int)(patch_ofset[patch_head[i]+j]>>kshift);
               int ref_rcdex = (int)(patch_ofset[patch_head[i]+k]&003);
               int ref_rcdex2 = (int)(patch_ofset[patch_head[i]+j]&003);
               printf("kmer0: %d %ld %ld %d %d %d %d %d %s %s\n",i,patch_array[patch_head[i]+k],patch_array[patch_head[i]+j],patch_list[i],ref_rcdex,ref_rcdex2,ref_offset,ref_offset2,(seq+idk)->name,(seq+idk2)->name);
                                                   */            
               patch_array[patch_head[i]+k] = -100;
               patch_array[patch_head[i]+j] = -100;
//		 reads_rept[idk] = 1;
//		 reads_rept[idk2] = 1;
//		 printf("wwww: %s %s\n",(seq+idk)->name,(seq+idk2)->name);
                 num_sect++;
                 j++;
             }
             k = j - 1;
          }
          ArraySort_Mix(patch_list[i],ray+patch_head[i],dex+patch_head[i]);
          for(k=0;k<patch_list[i];k++)
          {
             j = k+1;
             num_sect = 1;
             while((k<patch_list[i])&&(patch_array[patch_head[i]+j]==patch_array[patch_head[i]+k]))
             {
/*               if(patch_array[patch_head[i]+k] == -100)
               {
                 int idk = patch_index[patch_head[i]+k];
                 int idk2 = patch_index[patch_head[i]+j];
                 int ref_offset = (int)(patch_ofset[patch_head[i]+k]>>kshift);
                 int ref_offset2 = (int)(patch_ofset[patch_head[i]+j]>>kshift);
                 int ref_rcdex = (int)(patch_ofset[patch_head[i]+k]&003);
                 int ref_rcdex2 = (int)(patch_ofset[patch_head[i]+j]&003);
                 printf("kmer1: %d %ld %ld %d %d %d %d %d %s %s\n",i,patch_array[patch_head[i]+k],patch_array[patch_head[i]+j],patch_list[i],ref_rcdex,ref_rcdex2,ref_offset,ref_offset2,(seq+idk)->name,(seq+idk2)->name);
               }       */
               num_sect++;
               j++;
             }
             num_kmers++;
             k = j - 1;
          }
        }
     }
     printf("ssaha: %ld %ld %ld\n",mhistcc,num_kmers,nsorts);
     if((kmer_cons= (char *)calloc(200,sizeof(char))) == NULL)
     {
       printf("ssaha: calloc - kmer_cons\n");
       exit(1);
     }
     if((rbase= (char *)calloc(200,sizeof(char))) == NULL)
     {
       printf("ssaha: calloc - kmer_cons\n");
       exit(1);
     }
     if((mm= (int *)calloc(20000,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - mm\n");
       exit(1);
     }
     if((kmer_qual= (int *)calloc(200,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - kmer_qual\n");
       exit(1);
     }
     if((m_pileup= (int *)calloc(200,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - m_pileup\n");
       exit(1);
     }
     if((num_kmers==0)||(mhistcc==0))
     {
       printf("No valid reads in the file!\n");
       exit(1);
     }
     if((kmer_nhits= (int *)calloc(num_kmers,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - kmer_nhits\n");
       exit(1);
     }
     if((kmer_idgrp= (int *)calloc(num_kmers,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - kmer_idgrp\n");
       exit(1);
     }
     if((kmer_locus= (int *)calloc(num_kmers,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - kmer_locus\n");
       exit(1);
     }
     if((kmer_masks= (int *)calloc(num_kmers,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - kmer_masks\n");
       exit(1);
     }
     if((kmer_links= (int *)calloc(num_kmers,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - kmer_links\n");
       exit(1);
     }
     if((kmer_nodes= (int *)calloc(num_kmers,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - kmer_nodes\n");
       exit(1);
     }
     if((sm_head= (int *)calloc(num_kmers+2000,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - sm_head\n");
       exit(1);
     }
     if((sm_list= (int *)calloc(num_kmers+2000,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - sm_list\n");
       exit(1);
     }
     if((kmer_baseF= (char *)calloc(num_kmers+2000,sizeof(char))) == NULL)
     {
       printf("ssaha: calloc - kmer_baseF\n");
       exit(1);
     }
     if((kmer_baseR= (char *)calloc(num_kmers+2000,sizeof(char))) == NULL)
     {
       printf("ssaha: calloc - kmer_baseR\n");
       exit(1);
     }
     if((repeat_word= (int *)calloc(num_kmers,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - repeat_word\n");
       exit(1);
     }
     if((repeat_next= (int *)calloc(num_kmers,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - repeat_next\n");
       exit(1);
     }
     if((path_kmer= (int *)calloc(read_length,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - path_kmer\n");
       exit(1);
     }
     if((path_cover= (int *)calloc(read_length,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - path_cover\n");
       exit(1);
     }
     if((path_next= (int *)calloc(read_length,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - path_next\n");
       exit(1);
     }
     if((path_idir= (int *)calloc(read_length,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - path_idir\n");
       exit(1);
     }
     if((path2_kmer= (int *)calloc(read_length,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - path2_kmer\n");
       exit(1);
     }
     if((path2_next= (int *)calloc(read_length,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - path2_next\n");
       exit(1);
     }
     if((path2_idir= (int *)calloc(read_length,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - path2_idir\n");
       exit(1);
     }
     if((path_cons= (char *)calloc(read_length,sizeof(char))) == NULL)
     {
       printf("ssaha: calloc - path_cons\n");
       exit(1);
     }
     if((path2_cons= (char *)calloc(read_length,sizeof(char))) == NULL)
     {
       printf("ssaha: calloc - path2_cons\n");
       exit(1);
     }

     memset(kmer_baseF,'n',num_kmers);
     memset(kmer_baseR,'n',num_kmers);
     num_kmers = 0;
     num_ukmer = 0;
     num_akmer = 0;
     for(i=0;i<nsorts;i++) 
     {
        if(patch_list[i] > 1)
        {
          for(k=0;k<patch_list[i];k++)
          {
             j = k+1;
             num_sect = 1;
             while((k<patch_list[i])&&(patch_array[patch_head[i]+j]==patch_array[patch_head[i]+k]))
             {
               num_sect++;
               j++;
             }
             if(num_sect<1800)
               kmer_nhits[num_kmers] = num_sect;
             else
               kmer_nhits[num_kmers] = 1800;
             kmer_idgrp[num_kmers] = i;
             kmer_locus[num_kmers] = k;
             kmer_masks[num_kmers] = 0;
             for(w=k;w<j;w++)
             {
                int idd = kmer_index[patch_head[i]+w];
                patch_2kmer[patch_head[i]+idd] = num_kmers;
             }
//             if((num_sect>=2)&&(num_sect<750))
             if((num_sect>=3)&&(num_sect<250))
             {
               for(w=k;w<j;w++)
               {
                  int idd = kmer_index[patch_head[i]+w];
                  int idk = patch_index[patch_head[i]+idd];
                  int ref_offset = (int)(patch_ofset[patch_head[i]+idd]>>kshift);
                  int ref_rcdex = (int)(patch_ofset[patch_head[i]+idd]&003);
//                  patch_2kmer[patch_head[i]+idd] = num_kmers;
                  if(read2kmer[idk] < 0)
                  {
                    read2kmer[idk] = num_kmers;
                  }
               }
               num_akmer = num_akmer + num_sect;
               num_ukmer++;
             }
             num_kmers++;
             k = j - 1;
          }
        }
     }

     printf("number kmers %d %d %d\n",num_kmers,num_ukmer,num_cover);
     num_cover = num_akmer/num_ukmer;
     printf("number kmers %d %d %d\n",num_kmers,num_ukmer,num_cover);
     if(num_outs == 0)
       num_outs = num_ukmer*25/3000;
     cons_reads=cmatrix(0,num_outs,0,read_length+1);
     printf("number of extended reads: %d\n",num_outs);

     num_kmers = 0;
     for(i=0;i<nsorts;i++)
     {
        if(patch_list[i] > 1)
        {
          for(k=0;k<patch_list[i];k++)
          {

             j = k+1;
             num_sect = 1;
             while((k<patch_list[i])&&(patch_array[patch_head[i]+j]==patch_array[patch_head[i]+k]))
             {
               num_sect++;
               j++;
             }
             for(w=k;w<j;w++)
             {
                long next_link,last_link;
                int idd = kmer_index[patch_head[i]+w];
                int idk = patch_index[patch_head[i]+idd];
                int ref_offset = (int)(patch_ofset[patch_head[i]+idd]>>kshift);
                int ref_rcdex = (int)(patch_ofset[patch_head[i]+idd]&003);
                int i_node = num_kmers;

                next_link = patch_nlink[patch_head[i]+idd];
                last_link = patch_plink[patch_head[i]+idd];
                if(next_link>=0)
                {
                  int uph_index = (int)(patch_nlink[patch_head[i]+idd]>>nshift);
                  int low_index = (int)(patch_nlink[patch_head[i]+idd]&hmask);
                  int i_kmer = patch_2kmer[patch_head[uph_index]+low_index];
                  int cidk = sm_list[i_node];
                  int match = 0,jj;

                  for(jj=0;jj<cidk;jj++)
                  {
                     if(mm[jj] == i_kmer)
                     {
                       match = 1;
                       break;
                     }
                  }
                  if(match==0)
                  {
                    mm[sm_list[i_node]] = i_kmer;
                    sm_list[i_node]++;
                  }
                }
                if(last_link>=0)
                {
                  int uph_index = (int)(patch_plink[patch_head[i]+idd]>>nshift);
                  int low_index = (int)(patch_plink[patch_head[i]+idd]&hmask);
                  int i_kmer = patch_2kmer[patch_head[uph_index]+low_index];
                  int cidk = sm_list[i_node];
                  int match = 0,jj;

                  for(jj=0;jj<cidk;jj++)
                  {
                     if(mm[jj] == i_kmer)
                     {
                       match = 1;
                       break;
                     }
                  }
                  if(match==0)
                  {
                    mm[sm_list[i_node]] = i_kmer;
                    sm_list[i_node]++;
                  }
                }
             }
             num_kmers++;
             k = j - 1;
          }
        }
     }

     sm_head[0] = 0;
     for(i=1;i<=num_kmers;i++)
        sm_head[i] = sm_head[i-1] + sm_list[i-1]; 
     num_akmer = sm_head[num_kmers] + 2000;  
     if((sm= (int *)calloc(num_akmer,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - sm\n");
       exit(1);
     }
     for(i=0;i<num_kmers;i++)
        sm_list[i] = 0;

     num_kmers = 0;
     for(i=0;i<nsorts;i++) 
     {
        if(patch_list[i] > 1)
        {
          for(k=0;k<patch_list[i];k++)
          {
             j = k+1;
             num_sect = 1;
             while((k<patch_list[i])&&(patch_array[patch_head[i]+j]==patch_array[patch_head[i]+k]))
             {
               num_sect++;
               j++;
             }
             for(w=k;w<j;w++)
             {
	        long next_link,last_link;
                int idd = kmer_index[patch_head[i]+w];
                int idk = patch_index[patch_head[i]+idd];
                int ref_offset = (int)(patch_ofset[patch_head[i]+idd]>>kshift);
                int ref_rcdex = (int)(patch_ofset[patch_head[i]+idd]&003);
                int i_node = num_kmers;

		next_link = patch_nlink[patch_head[i]+idd];
		last_link = patch_plink[patch_head[i]+idd];
		if(next_link>=0)
		{
		  int uph_index = (int)(patch_nlink[patch_head[i]+idd]>>nshift);
		  int low_index = (int)(patch_nlink[patch_head[i]+idd]&hmask);
		  int i_kmer = patch_2kmer[patch_head[uph_index]+low_index];
		  int cidk = sm_list[i_node];
		  int offset = sm_head[i_node];
		  int match = 0,jj;

                  for(jj=0;jj<cidk;jj++)
                  {
		     if(sm[offset+jj] == i_kmer)
		     {
		       match = 1;
		       break;
		     }
		  }
		  if(match==0)
		  {
		    sm[offset+sm_list[i_node]] = i_kmer;
		    sm_list[i_node]++;
		  }
		}
		if(last_link>=0)
		{
		  int uph_index = (int)(patch_plink[patch_head[i]+idd]>>nshift);
		  int low_index = (int)(patch_plink[patch_head[i]+idd]&hmask);
		  int i_kmer = patch_2kmer[patch_head[uph_index]+low_index];
		  int cidk = sm_list[i_node];
		  int offset = sm_head[i_node];
		  int match = 0,jj;

                  for(jj=0;jj<cidk;jj++)
                  {
		     if(sm[offset+jj] == i_kmer)
		     {
		       match = 1;
		       break;
		     }
		  }
		  if(match==0)
		  {
		    sm[offset+sm_list[i_node]] = i_kmer;
		    sm_list[i_node]++;
		  }
		}
             }
             num_kmers++;
             k = j - 1;
          }
        }
     }

     for(i=0;i<num_kmers;i++)
     {
        long kmer_long;
        int num_rd = 0;
        int rd_array[sm_list[i]];

        j = kmer_idgrp[i];
        k = kmer_locus[i];
        kmer_long = patch_array[patch_head[j]+k];
        for(j=0;j<sm_list[i];j++)
        {
           int num_jkb = kmer_nhits[sm[sm_head[i]+j]];
           int short_flag,id_node,walk_len;
           if(num_jkb==1)
           {
             int r_off;
             char qbase;
             int jj = kmer_idgrp[sm[sm_head[i]+j]];
             int kk = kmer_locus[sm[sm_head[i]+j]];
             int idd = kmer_index[patch_head[jj]+kk];
             int idp = patch_index[patch_head[jj]+idd];
             int ref_offset = (int)(patch_ofset[patch_head[jj]+idd]>>kshift);
             int ref_rcdex = (int)(patch_ofset[patch_head[jj]+idd]&003);
             if(ref_rcdex == 1)
               r_off = ref_offset - 1 + n_Sbase;
             else
               r_off = (seq+idp)->length - 1 - ref_offset;
             if((seq+idp)->qual[r_off]>qthresh)
             {
               qbase = (seq+idp)->data[r_off];
               id_node = sm[sm_head[i]+j];
               if(sm_list[id_node]>0) 
                 Short_Path(i,id_node,n_Sbase,sm,sm_list,&short_flag,&walk_len,kmer_last);
               else
                 short_flag = 1;
               if(short_flag == 0)
               {
                 rd_array[num_rd] = sm[sm_head[i]+j];
                 num_rd++;
               }
             }
             else
             {
               qbase = tolower((seq+idp)->data[r_off]);
             }
           }
           else if(num_jkb==2)
           {
             int idp = 0;
             id_node = sm[sm_head[i]+j];
             if(sm_list[id_node]>0) 
               Short_Path(i,id_node,n_Sbase,sm,sm_list,&short_flag,&walk_len,kmer_last);
             else
               short_flag = 1;
             if(short_flag == 0)
             {
               rd_array[num_rd] = sm[sm_head[i]+j];
               num_rd++;
             }
           }
           else
           {
             rd_array[num_rd] = sm[sm_head[i]+j];
             num_rd++;
           }
        }
        if(num_rd<sm_list[i])
        {
          for(j=0;j<num_rd;j++)
            sm[sm_head[i]+j] = rd_array[j];
          sm_list[i] = num_rd;
        }
     }

     for(i=0;i<num_kmers;i++)
     {
        long kmer_long;
        int num_rd = 0;
        int rd_array[sm_list[i]];

        j = kmer_idgrp[i];
        k = kmer_locus[i];
        kmer_long = patch_array[patch_head[j]+k];
        for(j=0;j<sm_list[i];j++)
        {
           int num_jkb = kmer_nhits[sm[sm_head[i]+j]];
           int short_flag,id_node,walk_len;
           if(num_jkb==1)
           {
             int r_off;
             char qbase;
             int jj = kmer_idgrp[sm[sm_head[i]+j]];
             int kk = kmer_locus[sm[sm_head[i]+j]];
             int idd = kmer_index[patch_head[jj]+kk];
             int idp = patch_index[patch_head[jj]+idd];
             int ref_offset = (int)(patch_ofset[patch_head[jj]+idd]>>kshift);
             int ref_rcdex = (int)(patch_ofset[patch_head[jj]+idd]&003);
             if(ref_rcdex == 1)
               r_off = ref_offset - 1 + n_Sbase;
             else
               r_off = (seq+idp)->length - 1 - ref_offset;
             if((seq+idp)->qual[r_off]>qthresh)
             {
               qbase = (seq+idp)->data[r_off];
               id_node = sm[sm_head[i]+j];
               if(sm_list[id_node]>0) 
                 Short_Path(i,id_node,n_Sbase,sm,sm_list,&short_flag,&walk_len,kmer_last);
               else
                 short_flag = 1;
               if(short_flag == 0)
               {
                 rd_array[num_rd] = sm[sm_head[i]+j];
                 num_rd++;
               }
             }
             else
             {
               qbase = tolower((seq+idp)->data[r_off]);
             }
           }
           else if(num_jkb==2)
           {
             int idp = 0;
             id_node = sm[sm_head[i]+j];
             if(sm_list[id_node]>0) 
               Short_Path(i,id_node,n_Sbase,sm,sm_list,&short_flag,&walk_len,kmer_last);
             else
               short_flag = 1;
             if(short_flag == 0)
             {
               rd_array[num_rd] = sm[sm_head[i]+j];
               num_rd++;
             }
           }
           else
           {
             rd_array[num_rd] = sm[sm_head[i]+j];
             num_rd++;
           }
        }
        if(num_rd<sm_list[i])
        {
          for(j=0;j<num_rd;j++)
            sm[sm_head[i]+j] = rd_array[j];
          sm_list[i] = num_rd;
        }
     }

     Mmask_array(mmask,n_Sbase);
/*   process kmer to find gaps in kmers   **/
/*     for(j=0;j<num_kmers;j++)
     {
        long kmer1 = 5887022159168,kmer_long;
	int num_mismatch;

        i = kmer_idgrp[j];
        k = kmer_locus[j];
        kmer_long = patch_array[patch_head[i]+k];
        Mismatch_Kmer(kmer1,kmer_long,n_Sbase,mmask,&num_mismatch);
        num_sect = kmer_nhits[j];
     }     */


     num_cons = 0;
     if(num_kmers<(3*num_outs))
       exit(1);
     for(iseq=0;iseq<num_outs;iseq++)
     {
        int locus,i_kmer;
	int seeds;
        i = (int)(drand48()*nRead);
        j = (int)(drand48()*num_outs);
	seeds = j%10;
	if((i+seeds) < nRead)
	  i = i+seeds;
        i_kmer = read2kmer[i];
        while((read2kmer[i]<=0)||(kmer_masks[i_kmer]==1))
        {
          i = (int)(drand48()*nRead);
          j = (int)(drand48()*num_outs);
	  seeds = j%10;
  	  if((i+seeds) < nRead)
	    i = i+seeds;
          i_kmer = read2kmer[i];
        }
        reads_mask[i] = 1;
        kmer_masks[i_kmer] = 1;
/*        for(j=locus;j<(locus+num_sect);j++)
        {
           int idd = kmer_index[patch_head[kmer_idgrp[i_kmer]]+j];
           int idk = patch_index[patch_head[kmer_idgrp[i_kmer]]+idd];
        }    */
     }
     w = read2kmer[0];
     memset(kmer_masks,0,4*num_kmers);

     if((node_head= (int *)calloc(num_kmers,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - node_head\n");
       exit(1);
     }
     if((node_list= (int *)calloc(num_kmers,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - node_list\n");
       exit(1);
     }
     if((node_index= (int *)calloc(num_kmers,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - node_index\n");
       exit(1);
     }
     if((kmer2node= (int *)calloc(num_kmers,sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - kmer2node\n");
       exit(1);
     }
     if((node_rept= (int *)calloc((2*num_kmers),sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - node_rept\n");
       exit(1);
     }
     if((node_word= (int *)calloc((n_depth),sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - node_word\n");
       exit(1);
     }
     if((node_next= (int *)calloc((n_depth),sizeof(int))) == NULL)
     {
       printf("ssaha: calloc - node_next\n");
       exit(1);
     }
     num_nodes = 0;

     i_contig  = -1;
     reads_used = 0;
     for(iseq=i_read;iseq<nRead;iseq++)
     {
       if(reads_mask[iseq]==1)
       {
         node_hits = 0;
	 if((iseq%5000) == 0)
	 {
	   memset(rpf_contig,0,4*nSeq);
	   memset(rpf_offset,0,4*nSeq);
	   memset(reads_post,0,4*nSeq);
	   memset(reads_mast,0,4*nSeq);
	   reads_used = 0;
	 }
	 for(j=0;j<reads_used;j++)
	 {
	    rpf_contig[reads_maps[j]] = 0;
	    rpf_offset[reads_maps[j]] = 0;
	    reads_post[reads_maps[j]] = 0;
	    reads_mast[reads_maps[j]] = 0;
	 }
         Euler_Path(iseq,nRead,n_depth,n_Sbase,seq,patch_ofset,patch_index,patch_2kmer,patch_array,patch_head,patch_nlink,patch_plink);
         if(strlen(cons_reads[i_contig]) < read_length)
         {
           int ikk = kmer2node[node_word[0]];
//           int ikj = node_rept[node_word[0]];
//           node_head[i_contig] = kmer2node[node_word[0]];
//           node_list[i_contig] = node_rept[node_word[0]];
           node_head[node_word[0]] = kmer2node[node_word[0]];
           node_list[node_word[0]] = strlen(cons_reads[i_contig]);
           node_index[node_word[0]] = i_contig;
         }
       }
     }

     e_contig  = i_contig;
     printf("Phase one finished! %d %d\n",i_contig,n_repeats);
     exts_reads=cmatrix(0,n_repeats,0,read_length*2);
//     exts_reads2=cmatrix(0,n_repeats,0,read_length*2);
     num_cons = 0;
     print_flag = 1;
     for(k=0;k<n_repeats;k++)
     {
        long kk1,kk2;
        j = repeat_next[k];
        if(kmer_nhits[j]>=2)
        {
          int idd,idk,break_kmer,move_len;
          i = kmer_idgrp[repeat_word[k]];
          w = kmer_locus[repeat_word[k]];
          kk1 = patch_array[patch_head[i]+w];
          i = kmer_idgrp[j];
          w = kmer_locus[j];
          kk2 = patch_array[patch_head[i]+w];
          idd = kmer_index[patch_head[i]+w]; 
          idk = patch_index[patch_head[i]+idd];
          iseq = idk;
          path_dir = 1;
          move_len = 700;
          printf("repeat_kmer: %d %d %ld %ld %s %ld %d\n",k,kmer_nhits[j],repeat_word[k],j,(seq+idk)->name,kk1,kmer_links[repeat_word[k]]);
          junction1 = repeat_word[k];
          junction2 = j;
          memset(path_kmer,0,read_length*4);
          memset(path_next,0,read_length*4);
          memset(path_idir,0,read_length*4);
          memset(path_cons,'\0',read_length);
          memset(path2_cons,'\0',read_length);
          n_node = k;
          i_contig  = k;
          Euler_Path2(&break_kmer,j,k,move_len,n_depth,n_Sbase,seq,patch_ofset,patch_index,patch_2kmer,patch_array,patch_head,patch_nlink,patch_plink);
          printf("case0: %d %d %d %d\n",k,j,break_kmer,num_cons);
          n_node = -1;
          if((num_cons<2)&&(break_kmer==repeat_word[k]))
          {
            int num_tkmer = 0,num_lcover = 0;
            int num2_cons;
	    int ext_len = 0;

            path_dir = 0;
            Euler_Path2(&break_kmer,j,k,move_len,n_depth,n_Sbase,seq,patch_ofset,patch_index,patch_2kmer,patch_array,patch_head,patch_nlink,patch_plink);
            num2_cons = num_cons;
            for(m=0;m<num2_cons;m++)
               num_tkmer = num_tkmer + path_cover[m];
            strcpy(path2_cons,path_cons);
            memset(path_cons,'\0',read_length);
          printf("case01: %d %d %d %d\n",break_kmer,num_cons,node_head[junction1],node_list[junction1]);
            path_dir = 1;
            Euler_Path2(&break_kmer,j,k,move_len,n_depth,n_Sbase,seq,patch_ofset,patch_index,patch_2kmer,patch_array,patch_head,patch_nlink,patch_plink);
          printf("case11: %d %d %d %d %d %d %d\n",k,num_cons,num2_cons,node_head[junction1],node_list[junction1],junction1,junction2);
            for(m=0;m<num2_cons;m++)
            {
               if     (path2_cons[num2_cons-m-1] == 'A')
                 exts_reads[k][m] = 'T';
               else if(path2_cons[num2_cons-m-1] == 'C')
                 exts_reads[k][m] = 'G';
               else if(path2_cons[num2_cons-m-1] == 'G')
                 exts_reads[k][m] = 'C';
               else if(path2_cons[num2_cons-m-1] == 'T')
                 exts_reads[k][m] = 'A';
               else
                 exts_reads[k][m] = path2_cons[num2_cons-m-1];
            }
            if(num2_cons>=n_Sbase)
            {
              for(m=n_Sbase;m<num_cons;m++)
              {
                 exts_reads[k][num2_cons+m-n_Sbase] = path_cons[m];
                 num_tkmer = num_tkmer + path_cover[m];
              }
            }
            else
            {
              for(m=0;m<num_cons;m++)
              {
                 exts_reads[k][m] = path_cons[m];
                 num_tkmer = num_tkmer + path_cover[m];
              }
            }
	    ext_len = strlen(exts_reads[k]);
	    if(ext_len>0)
	       num_lcover = num_tkmer/ext_len;
	    else
	       num_lcover = 0;
            printf("cover:1 %d %d %d %d %d\n",k,num_cons,num2_cons,num_cover,num_lcover);
            printf("read: %d %d %d %d www1\n",k,ext_len,num_cons,num2_cons);
//            if(num2_cons>3)
//              strncpy(exts_reads2[k],exts_reads[k],strlen(exts_reads[k])-2); 
          }
          else
          {
            int num_tkmer = 0,num_lcover = 0;
            int num2_cons = num_cons;
	    int ext_len = 0;
            path_dir = 0;
            for(m=0;m<num2_cons;m++)
               num_tkmer = num_tkmer + path_cover[m];
            strncpy(path2_cons,path_cons,num2_cons);
            memset(path_cons,'\0',read_length);
            Euler_Path2(&break_kmer,j,k,move_len,n_depth,n_Sbase,seq,patch_ofset,patch_index,patch_2kmer,patch_array,patch_head,patch_nlink,patch_plink);
            for(m=0;m<num2_cons;m++)
            {
               if     (path2_cons[num2_cons-m-1] == 'A')
                 exts_reads[k][m] = 'T';
               else if(path2_cons[num2_cons-m-1] == 'C')
                 exts_reads[k][m] = 'G';
               else if(path2_cons[num2_cons-m-1] == 'G')
                 exts_reads[k][m] = 'C';
               else if(path2_cons[num2_cons-m-1] == 'T')
                 exts_reads[k][m] = 'A';
               else
                 exts_reads[k][m] = path2_cons[num2_cons-m-1];
            }
          printf("case02: %d %d %s\n",break_kmer,num2_cons,exts_reads[k]);
            if(num2_cons>=n_Sbase)
            {
              for(m=n_Sbase;m<num_cons;m++)
              {
                 exts_reads[k][num2_cons+m-n_Sbase] = path_cons[m];
                 num_tkmer = num_tkmer + path_cover[m];
              }
            }
            else
            {
              for(m=0;m<num_cons;m++)
              {
                 exts_reads[k][m] = path_cons[m];
                 num_tkmer = num_tkmer + path_cover[m];
              }
            }
	    ext_len = strlen(exts_reads[k]);
	    if(ext_len>0)
	       num_lcover = num_tkmer/ext_len;
	    else
	       num_lcover = 0;
            printf("cover:2 %d %d %d %d %d %s\n",k,num_cons,num2_cons,num_cover,num_lcover,exts_reads[k]);
            printf("read: %d %d %d %d www2\n",k,ext_len,num_cons,num2_cons);
//            if(num2_cons>3)
//              strncpy(exts_reads2[k],exts_reads[k],strlen(exts_reads[k])-2); 
          }
          
        }
     }

     printf("Phase two finished! %d %d\n",i_contig,n_repeats);
     if((namef = fopen(argv[args+2],"w")) == NULL)
     {          
       printf("namef - args name\n");
       exit(1);
     }
     for(j=0;j<e_contig;j++)
     {
        int rc,st,ed;
        ed = strlen(cons_reads[j]);
        if(ed>50)
        {
          fprintf(namef,"@extend_seq_%08d\n",j);
          for(rc=0;rc<ed;rc++)
             fprintf(namef,"%c",cons_reads[j][rc]);
          fprintf(namef,"\n");
          fprintf(namef,"+extend_seq_%08d\n",j);
          putc(041,namef);
          for(rc=1;rc<20;rc++)
             putc(5+041,namef);
          for(rc=20;rc<(ed-20);rc++)
             putc(40+041,namef);
          for(rc=(ed-20);rc<ed;rc++)
             putc(5+041,namef);
          fprintf(namef,"\n");
        }
     }
     for(j=0;j<n_repeats;j++)
     {
        int rc,st,ed;
        ed = strlen(exts_reads[j]);
        if(ed>50)
        {
          fprintf(namef,"@extra_seq_%08d\n",j);
          for(rc=0;rc<ed;rc++)
             fprintf(namef,"%c",exts_reads[j][rc]);
          fprintf(namef,"\n");
          fprintf(namef,"+extra_seq_%08d\n",j);
          putc(041,namef);
          for(rc=1;rc<20;rc++)
             putc(5+041,namef);
          for(rc=20;rc<(ed-20);rc++)
             putc(40+041,namef);
          for(rc=(ed-20);rc<ed;rc++)
             putc(5+041,namef);
          fprintf(namef,"\n");
        }
     }

/*     for(j=0;j<n_repeats;j++)
     {
        int rc,st,ed;
        ed = strlen(exts_reads2[j]);
        if(ed>50)
        {
          fprintf(namef,"@extra2_seq_%08d\n",j);
          for(rc=0;rc<ed;rc++)
             fprintf(namef,"%c",exts_reads2[j][rc]);
          fprintf(namef,"\n");
          fprintf(namef,"+extra2_seq_%08d\n",j);
          putc(041,namef);
          for(rc=1;rc<ed;rc++)
          {
             putc(40+041,namef);
          }
          fprintf(namef,"\n");
        }
     } */
     fclose(namef);
     
     free(patch_list);
     free(patch_head);

     free(patch_array);
     free(patch_index);
}

/*   Subroutine to walk euler path to extend read length   */
/* =============================================  */
void Euler_Path(int iseq,int nRead,int n_depth,int n_Sbase,fasta *seq,int *patch_ofset,int *patch_index,int *patch_2kmer,long *patch_array,long *patch_head,long *patch_nlink,long *patch_plink)
/* =============================================  */
{
     int i,j,k,w,i_node;
     void ArraySort_Int2(int n, int *arr, int *brr);
     void ArraySort_Mix(int n, long *arr, int *brr);
     int Break_Point,B_joint,hit_depth,kmer_last[3];
     void Check_454maps(int num,int nRead,int bulb_flag,int B_joint,int *map_454post,int *map_454next,int *map_454path,int *map_454qual,int max_kmerhits,int max_kmernext,int *Break_454,int *map_kmernext);
     void Check_Break(int nRead,int *id_post2, int n_depth, int *Break_Point,int *B_joint);
     void Check_Break_c3(int nRead,int *id_post2, int n_depth, int *Break_Point);
     void Check_Break_pile(int nRead,int *id_post2, int n_depth, int *Break_Point);
     int num_sect,num_sect2;
     fasta *seqp;
     long kmer_long;
     int kk,next_kmer,km,i_kmer,last_hit=0,num_454reads;
     int ref_offset = 0,ref_rcdex,stop_flag;
     int next_move1,next_move2,max_start;
     float path_rate = 1.0;
     int kshift = 2,single_hit;
     int r_depth = n_depth*2.5;
     long next_link,last_link,id_nkmer[n_depth];
     char cons_base,base_pile[n_depth],kmerch[n_Sbase],kmerrc[n_Sbase];
     int reads_last[r_depth],reads_kmer[r_depth],reads_idex[r_depth],reads_rcdx[r_depth],reads_rclt[r_depth];
     int reads_loci[r_depth],reads_lolt[r_depth],reads_dupy[r_depth],mismatch_loci[100];
     int id_read[n_depth],id_post[n_depth],id_post2[n_depth],id_rddex[n_depth],id_offset[n_depth]; //id_read will never exceed kmer_depth
     int num_lasthit,ref_rcdex2,next_kmer_p,next_kmer_p2;
     int word_hit,next_flag = 0,prev_kmer,num_maxs;
     void Short_Path(int idp,int id_node,int n_Sbase,int *sm,int *sm_list,int *short_flag,int *walk_len,int *kmer_last);
     void Maximum_Occurr(int nHit,int *hit_array, int *hit_index,int *hit_list);
     void Patch_Occurr(int nHit,int *hit_array, int *hit_index,int *hit_next,int *hit_list);
     void Read_Align(fasta *seq, int *kmer_index, int *patch_index, char **m_align,int **m_qualy, int *patch_ofset,long *patch_head,int *id_read, int *id_post, int *id_rcdx,int num_sect,int kshift,int max_read,int n_Sbase,int n_depth,int i,int k);
     void Read_AlignFuzzy(fasta *seq, int *kmer_index, int *patch_index, char **m_align,char *fuzzy_ray,int **m_qualy, int *patch_ofset,long *patch_head,int *id_read, int *id_post, int *id_rcdx,int num_sect,int kshift,int max_read,int n_Sbase,int n_depth,int i,int k,int num_sum,int num_mismatch,int *mismatch_loci,long kmer1);
     void Read_Pileup(fasta *seq, int *kmer_index, int *patch_index, char **m_align,int **m_qualy, int *patch_ofset,int *patch_2kmer,long *patch_array,long *patch_head,long *patch_nlink,long *patch_plink,long *id_nkmer,int *id_read, int *id_post, int *id_rcdx,int *id_post2,int *id_rddex,int *id_offset,int *reads_kmer,int *reads_last,int *reads_idex,int *reads_rcdx,int *reads_rclt,int *reads_loci,int *reads_lolt,int *reads_dupy,int num_sect,int kshift,int max_len,int max_read,int n_Sbase,int n_depth,int num_lasthit,long kmer_long, int i, int k,int i_node,int *last_hit);
     void Call_Mismatchkmers(long *patch_array,long *patch_head,long kmer1,int n_Sbase);
     int pair_next[n_depth],pair_index[n_depth],pair_patch[n_depth],hit_list[4];
     int pair_pnext[n_depth],pair_pgaps[n_depth],pair_quals[n_depth],pair_pdist[n_depth];
     long mmask[100];
     void Mmask_array(long *mmask,int n_Sbase);
     void Mismatch_Kmer(long kmer1, long kmer2,int n_Sbase,long *mmask,int *num_mismatch,int *mismatch_loci);
     void Check_Break_Fuzzy(int ikk,int num_sectsum,int *id_post2,int n_depth,char *fuzzy_ray,int *Break_Point,int *got_next);
     void Get_Kmerbase(long IntSeg,int n_Sbase,char *kmerch,char *kmerrc);

        i_contig++;
	hit_depth = 0;
        single_hit = 0;
//        if(num_cons>=100)
        if(num_cons>=2)
        {
          int num_tkmer = 0,num_lcover = 0;
          for(j=0;j<num_cons;j++)
             num_tkmer = num_tkmer + path_cover[j];
          num_lcover = num_tkmer/num_cons;
          printf("read: %d %d %d %d %d %s\n",i_seq,num_cons,i_contig,num_cover,num_lcover,(seq+n_list[i_contig-1])->name);
        }
        w = read2kmer[iseq];
        printf("contig: %d %d %s %d %d\n",i_contig,iseq,(seq+iseq)->name,w,reads_used);
	reads_used = 0;
	i_node = w;
        next_kmer = w;
        next_kmer_p = -1;
        next_kmer_p2 = -1;
        prev_kmer = 0;
	weird_kmer = 0;
        if(print_ctg)
        {
          if(abs(i_contig-print_ctg)<2)
            print_flag = 1;
          else
            print_flag = 0;
        }
        n_list[i_contig] = iseq;
        num_cons = 0;
        stop_flag=0;
        memset(reads_last,0,4*r_depth);
        memset(reads_kmer,0,4*r_depth);
        memset(reads_idex,0,4*r_depth);
        memset(reads_rcdx,0,4*r_depth);
        memset(reads_rclt,0,4*r_depth);
        memset(reads_loci,0,4*r_depth);
        memset(reads_lolt,0,4*r_depth);
        num_lasthit = 0;
        rcdex_pres = 0;

	move_dir = 1;
        while((read2kmer[iseq]>=0)&&(stop_flag == 0))
        {
           long i_base;
           i = kmer_idgrp[next_kmer];
           k = kmer_locus[next_kmer];
           kmer_long = patch_array[patch_head[i]+k];
           num_sect = kmer_nhits[next_kmer];

	   i_node = next_kmer;
           num_cons++;
           if(num_cons>=read_length)
             stop_flag  = 1;
           if((num_sect <= kmer_depth)&&(num_sect>=1))
           {
             if((num_sect>=1)&&(kmer_long!= -100))
             {
               int dupst_flag = 0;
               int exten_len = 0;
               int duped_flag = 0,id_rcdx[kmer_depth];
               int idd, idk,idd2,idk2,ref_offset2;
               int rc,st,ed,max_len,out_flag = 1,max_read,bulb_flag;
               float rate = 0.0;
               int max_move1,max_move2,max_node2,jj,ikk,rcdex1,rcdex2,rcdex_nows,num_hits;

               max_read = 0;
               for(kk=0;kk<num_sect;kk++)
               {
                  idd = kmer_index[patch_head[i]+kk+k];
                  idk = patch_index[patch_head[i]+idd];
                  seqp = seq + idk;
                  if(seqp->length>max_read)
		  {
                    max_read = seqp->length;
		    if(seqp->length>=50)
		      max_read = 50;
		  }
               }

               memset(base_pile,'\0',n_depth);
               Read_Align(seq,kmer_index,patch_index,m_align,m_qualy,patch_ofset,patch_head,id_read,id_post,id_rcdx,num_sect,kshift,max_read,n_Sbase,n_depth,i,k);

               ArraySort_Int2(num_sect,id_post,id_read);
               idd = kmer_index[patch_head[i]+k+id_read[0]];
               idk = patch_index[patch_head[i]+idd];
               ref_offset = (int)(patch_ofset[patch_head[i]+idd]>>kshift);
               idd2 = kmer_index[patch_head[i]+k+id_read[num_sect-1]];
               ref_offset2 = (int)(patch_ofset[patch_head[i]+idd2]>>kshift);
               max_len = max_read + max_read - n_Sbase - ref_offset; 
	       if(max_len < 20)
	         max_len = 50;
               exten_len = ref_offset2-ref_offset;
               if(exten_len>=0)
               {
	           hit_ikk = 0;
                   Read_Pileup(seq,kmer_index,patch_index,m_align,m_qualy,patch_ofset,patch_2kmer,patch_array,patch_head,patch_nlink,patch_plink,id_nkmer,id_read,id_post,id_rcdx,id_post2,id_rddex,id_offset,reads_kmer,reads_last,reads_idex,reads_rcdx,reads_rclt,reads_loci,reads_lolt,reads_dupy,num_sect,kshift,max_len,max_read,n_Sbase,n_depth,num_lasthit,kmer_long,i,k,i_node,&last_hit);
                   ikk = hit_ikk;
                   if(ikk<=1)
                     single_hit++;
                   else
                     single_hit = 0;
                   if(num_cons>1)
                     kmer_links[prev_kmer] = 1;
                   prev_kmer = next_kmer;
                   path_cover[num_cons-1] = ikk;
                   ArraySort_Int2(ikk,id_post,id_read);
                   max_move1 = 0;
                   max_move2 = 0;
                   max_node2 = 0;
                   max_start = 0;
		   num_maxs = 0;
                   for(kk=0;kk<ikk;kk++)
                   {
                      jj = kk+1;
		      num_454reads = 1;
                      while((jj<ikk)&&(id_post[kk]==id_post[jj]))
                      {
		        if((seq+id_rddex[id_read[jj]])->length>30)
			  num_454reads++;
                        jj++;
                      }
                      if((jj-kk)>max_move1)
                      {
                        if(max_move1 > 0)
			{
                          max_move2 = max_move1;
			  max_node2 = id_post[kk];
			}
                        max_move1 = jj-kk;
                        max_start = kk;
                      }
                      if(((jj-kk)>max_move2)&&(max_start!=kk))
		      {
                        max_move2 = jj-kk;
			max_node2 = id_post[max_start];
		      }
                 if(print_flag)
                      printf("max: %d %d %d %d %d %d %d\n",kk,id_post[kk],i_contig,max_move1,max_move2,max_start,num_454reads);
                      num_maxs++;
                      kk = jj - 1;
                   }

                   if((max_move1<ikk)&&(max_move2==0))
                     max_move2 = ikk - max_move1;
                   num_lasthit = last_hit;

                   for(kk=0;kk<ikk;kk++)
                   {
                      id_read[kk] = kk;
                      base_pile[kk] = neib_kmer[kk][0];
                   }
                   ArraySort_Mix(ikk,id_nkmer,id_read);
		   /* Pay attention to this  */
		   bulb_flag = 0;
                   if((max_move2>=2)&&(sm_list[i_node]>=3))
                   {
		     int short_flag=0,n_long = 0,walk_len,kj;
                     int bulb_list[sm_list[i_node]];
                     int bulb_list1[sm_list[i_node]];
                     int bulb_list2[sm_list[i_node]];
                    
                    printf("check: %d %d %d %d\n",max_move2,i_node,max_node2,sm_list[i_node]);
                     if(sm_list[i_node]==2)
                       short_flag = 1;
                     else
                     {
                       n_long = 0;
                       for(kk=0;kk<sm_list[i_node];kk++)
                       {
                          max_node2 = sm[sm_head[i_node]+kk];
                          Short_Path(i_node,max_node2,n_Sbase,sm,sm_list,&short_flag,&walk_len,kmer_last);
                 if(print_flag)
		          printf("node: %d %d %d %d %d %d %d\n",kk,max_node2,short_flag,walk_len,kmer_last[0],kmer_last[1],kmer_last[2]);
			  bulb_list[kk] = kmer_last[0];
			  bulb_list1[kk] = kmer_last[1];
			  bulb_list2[kk] = kmer_last[2];
                          if(short_flag == 0)
                            n_long++;
                       }

		       bulb_flag = 0;
		       if(sm_list[i_node]==3)
		       {
                         for(kk=0;kk<sm_list[i_node];kk++)
		         {
		            for(kj=0;kj<sm_list[i_node];kj++)
		  	    {
			       if((bulb_list[kj]!=0)&&(kj!=kk)&&(bulb_list[kj]==bulb_list[kk]))
			       {
			         if((bulb_list1[kj]!=bulb_list1[kk])&&(bulb_list2[kj]!=bulb_list2[kk]))
			         {
			           bulb_flag = 1;
			         }
			       }
			    }
		         }
		       }
                       if((sm_list[i_node]==3)&&(n_long==2))
                       {
                         short_flag = 1;
                         max_move2 = 0;
                       }
                       else
                         short_flag=0;      
                    printf("checkwww: %d %d %d %d %d\n",max_move2,n_long,bulb_flag,i_node,sm_list[i_node]);
                     }
		     if(bulb_flag)
                       printf("crush bulb: %d %d %d %d\n",max_move2,bulb_flag,i_node,sm_list[i_node]);
		     if(short_flag==0)
		     {
                       Check_Break_pile(ikk,id_post2,n_depth,&Break_Point); 
		       if((Break_Point==1)&&(rp_insert==0))
		         stop_flag = 1;
		     }
		     else
		       max_move2 = 0;
                   }
                   next_move1 = 0;
                   next_move2 = 0;
                   for(kk=0;kk<ikk;kk++)
                   {
                      jj = kk+1;
                      while((jj<ikk)&&(id_nkmer[kk]==id_nkmer[jj]))
                      {
                        jj++;
                      }
                      if((jj-kk)>next_move1)
                      {
                        if(next_move1>0)
                          next_move2 = next_move1;
                        next_move1 = jj-kk;
                      }
                      kk = jj - 1;
                   }
		   if((ikk==0)||(max_move1==0))
		   {
		     ikk = 1;
		     max_move1 = 1;
		     stop_flag = 1;
		   }
                   path_rate = next_move1;
                   path_rate = path_rate/ikk;
                   rate  = max_move2;
                   rate  = rate/max_move1;
                   i_seq = iseq;
		   if((ikk==3)&&(max_move1==2))
		   {
                     Check_Break_c3(ikk,id_post2,n_depth,&Break_Point);
		     if(Break_Point==1)
		       stop_flag = 1;
		   }
                   num_454reads = 0;
		   for(kk=0;kk<ikk;kk++)
		   {
		      if((seq+id_rddex[id_read[kk]])->length>50)
		        num_454reads++;
		   }
		   hit_depth = ikk;
                   if((path_rate < 0.65)||(rate > 0.52)||(max_move2>=set_mismatch))
                   {
                     int m,n_pairs = 0;
                     int pair_flag  = 0;
                     int NUM_pairs=2;
                     int n_rept = n_repeats;
                     int num_tkmer = 0,num_lcover = 0;
                     int ctgkmer_next[ikk],i_rept = 0;
                     int map_post[ikk],map_kmer[ikk];
	             int num_ppatch = 0;
		     int max_offset,max_kmer;
		     int fuzzy_flag = 0;
		     int num_454hits, num_454maps,num_454pats,num_454kmer,num_454allhits,Break_454;
		     int map_454post[ikk],map_454next[ikk],map_454path[ikk],map_454qual[ikk];
                     int num_454negs,num_454cuts,num_454kmer2;
                     int max_kmerhits,max_kmernext;
                     float rate2;

		     max_kmerhits = 0;
		     max_kmernext = 0;
                     num_454maps = 0;
                     num_454pats = 0;
                     num_454negs = 0;
                     num_454allhits = 0;
                     printf("repeat: %d %f %f %d %d\n",ikk,path_rate,rate,max_move2,num_cons);
                     Check_Break(ikk,id_post2,n_depth,&Break_Point,&B_joint);
                     for(kk=0;kk<num_cons;kk++)
                        num_tkmer = num_tkmer + path_cover[kk];
                     if(num_cons>0)
                       num_lcover = num_tkmer/num_cons;
                     else
                       num_lcover = 0;
		     if((Break_Point==1)&&(max_move1<=2)&&(check_fuzzy))
		     {
		       int num_sectsum = 0;
		       int num_sectsum1 = 0;
		       int got_next;
		  
                       memset(base_pile,'\0',n_depth);
		       Mmask_array(mmask,n_Sbase);
                       printf("gap: %d %f %ld %d %d\n",ikk,path_rate,kmer_long,max_move1,max_move2);
		       for(kk=0;kk<num_kmers;kk++)
	               {
     	                  int num_mismatch,num_sect2;
			  long kmer1;
                          int ix = kmer_idgrp[kk];
			  int iy = kmer_locus[kk];

			  kmer1 = patch_array[patch_head[ix]+iy];
   	                  Mismatch_Kmer(kmer_long,kmer1,n_Sbase,mmask,&num_mismatch,mismatch_loci);
			  num_sect2 = kmer_nhits[kk];
//			    printf("mkmer: %ld %d %d\n",kmer1,num_mismatch,num_sect2);
			  if((num_mismatch <=2)&&(num_sect2<50))
			  {
			    if(num_sectsum>20)
			    {
			      if((num_mismatch <= 1)&&(num_sectsum<(n_depth-2)))
			      {
                                Read_AlignFuzzy(seq,kmer_index,patch_index,m_align,base_pile,m_qualy,patch_ofset,patch_head,id_read,id_post,id_rcdx,num_sect2,kshift,max_read,n_Sbase,n_depth,ix,iy,num_sectsum,num_mismatch,mismatch_loci,kmer1);
			        num_sectsum = num_sectsum + num_sect2;
			      }
			    }
			    else
			    {
                                Read_AlignFuzzy(seq,kmer_index,patch_index,m_align,base_pile,m_qualy,patch_ofset,patch_head,id_read,id_post,id_rcdx,num_sect2,kshift,max_read,n_Sbase,n_depth,ix,iy,num_sectsum,num_mismatch,mismatch_loci,kmer1);
			        num_sectsum = num_sectsum + num_sect2;
			    }
			  }
		       }
                       Check_Break_Fuzzy(ikk,num_sectsum,id_post2,n_depth,base_pile,&Break_Point,&got_next);
		       if(Break_Point==1)
		         stop_flag = 1;
		       else
		       {
                         printf("got: %d %d\n",ikk,got_next);
			 single_hit = 0;
			 fuzzy_flag = 1;
                         next_kmer = got_next;
		       }
		     }
                     else if(Break_Point==1)
                     {
                       for(kk=0;kk<ikk;kk++)
                       {
                          jj = kk+1;
                          while((jj<ikk)&&(id_nkmer[kk]==id_nkmer[jj]))
                          {
                            jj++;
                          }
			  max_offset = 0;
			  max_kmer = 0;
			  num_454hits = 0;
			  num_454cuts = 0;
                          for(m=kk;m<jj;m++)
                          {
                             int pdk = id_read[m];
                             int rd = mat_pair[id_rddex[pdk]];
                             int plen = abs(num_cons-rpf_offset[rd]+n_Sbase);
                             int gap = abs(rp_insert - abs(num_cons-rpf_offset[rd])-n_Sbase);
			     int map_dis;

                             if((reads_rept[id_rddex[pdk]]==0)&&(neib_qual[pdk]>=2))
			     {
			       int bridge = 0;
			       if(reads_mast[id_rddex[pdk]]>0)
			       {
			         if((reads_post[id_rddex[pdk]]==id_offset[pdk])&&(reads_post[id_rddex[pdk]]>=num_overlap))
				   bridge = 1;
			       }
                               if(reads_post[id_rddex[pdk]]>max_kmerhits)
			       {
			         max_kmerhits = reads_post[id_rddex[pdk]];
				 max_kmernext = id_post2[id_read[m]];
			       }
			       if(rp_insert > 0) bridge = 0;
                               if(bridge)
			       {
                                 printf("454_bridge: %-2d %-25s %d %d %d %d %d %d\n",m,(seq+id_rddex[pdk])->name,id_post2[id_read[m]],num_cover,num_lcover,reads_post[id_rddex[pdk]],(seq+id_rddex[pdk])->length,abs(id_offset[pdk]-reads_post[id_rddex[pdk]]));
                                 map_454post[num_454allhits] = reads_post[id_rddex[pdk]];
                                 map_454next[num_454allhits] = id_post2[id_read[m]];
                                 map_454qual[num_454allhits] = neib_qual[pdk];
                                 map_454path[num_454allhits] = num_454maps;
				 num_454allhits++;
			         num_454hits++;
			       }
			       if((reads_post[id_rddex[pdk]] == num_cons)&&(num_cons>=20))
			       {
			         printf("454_cut: %d %d\n",reads_post[id_rddex[pdk]],num_cons);
				 num_454cuts++;
			       }
			       if((reads_post[id_rddex[pdk]]-id_offset[pdk])>n_Sbase)
			         num_454negs++;
			     }
                 if(print_flag)
                             printf("pass: %-2d %-25s %-42s %ld %ld %d %d %d %d %d %d %d %d %d\n",m,(seq+id_rddex[pdk])->name,m_align[pdk],kmer_long,i_node,id_post2[id_read[m]],num_cover,num_lcover,reads_post[id_rddex[pdk]],(seq+id_rddex[pdk])->length,id_offset[pdk],reads_mast[id_rddex[pdk]],reads_rept[id_rddex[pdk]],neib_qual[pdk]);
                             if((num_pair[rd]==1)&&(rpf_contig[rd] == i_contig)&&(gap<100)&&(rpf_offset[rd]>0))
//                             if((rpf_contig[rd] == i_contig)&&(num_pair[rd]==1)&&(rpf_offset[rd]>0))
                             {
                               printf("pair: %s %s %d %d %d %d %d\n",(seq+id_rddex[pdk])->name,(seq+rd)->name,rpf_contig[rd],rpf_offset[rd],gap,plen,num_cons);
                               pair_next[n_pairs] = id_post2[id_read[m]];
                               pair_index[n_pairs] = n_pairs;
                               pair_pnext[n_pairs] = id_post2[id_read[m]];
                               pair_pgaps[n_pairs] = gap;
                               pair_pdist[n_pairs] = plen;
			       pair_quals[n_pairs] = neib_qual[pdk];
                               pair_patch[n_pairs] = num_ppatch;
                               n_pairs++;
                             }
			     if(reads_post[id_rddex[pdk]]>max_offset)
			       max_offset = reads_post[id_rddex[pdk]];
			     max_kmer = id_post2[id_read[m]];
                          }
                          if(((jj-kk)>=2)&&(kmer_masks[next_kmer]==0)&&(num_454allhits==0))
                          {
                             int pdk = id_read[kk];
                             int rd = mat_pair[id_rddex[pdk]];
                             int idgrp = kmer_idgrp[id_post2[id_read[kk]]];
                             int locus = kmer_locus[id_post2[id_read[kk]]];
                             int iidd =  kmer_index[patch_head[idgrp]+locus];
                             int iidk = patch_index[patch_head[idgrp]+iidd];

/*                             if((next_kmer==next_kmer_p)&&(id_post2[id_read[kk]]==next_kmer_p2))
                             {
                               max_move2 = 10;
                               Break_Point = 1;
                               kk = ikk;
                               printf("next2break: %d %d %d %d %d\n",ikk,i_contig,next_kmer,next_kmer_p,next_kmer_p2);
                             }
 
                             if(next_flag==0)
                             {
                               next_kmer_p = next_kmer;
                               next_kmer_p2 = id_post2[id_read[kk]];
                             next_flag = 1;
                             }       */
                             if(Break_Point!=0)
                             {
                               path_rate = max_move1;
                               if(max_move2==0)
                                 path_rate = 2.0;
                               else
                                 path_rate = path_rate/max_move2;
                               if((path_rate<1.4)||(max_move2>=set_mismatch))
			       {
                                 ctgkmer_next[i_rept] = id_post2[id_read[kk]];
                                 kmer2node[next_kmer] = n_rept; 
//                                 node_rept[next_kmer] = n_repeats-n_rept+1; 
                                 i_rept++;
                               }
                             }
                          }
                          if((jj-kk)>=2)
                          {
                            node_word[node_hits] = next_kmer;
                            node_next[node_hits] = id_post2[id_read[kk]];
			    map_post[num_ppatch] = max_offset;
			    map_kmer[num_ppatch] = max_kmer;
			    if(num_454hits>0)
			    {
			      num_454kmer = id_post2[id_read[kk]];
			      num_454maps++;
			    }
			    if(num_454cuts>0)
			    {
			      num_454kmer2 = id_post2[id_read[kk]];
			      num_454pats++;
			    }
			    printf("down: %d %d %d %d %d %d\n",ikk,i_contig,num_cons,num_ppatch,map_post[num_ppatch],num_454hits); 
                            node_hits++;
			    num_ppatch++;
                          }
                          word_hit = 0;
                          for(m=0;m<node_hits;m++)
                          {
                             if(next_kmer==node_word[m])
                               word_hit++;
                          }
                          if(word_hit>3)
                          {
                             max_move2 = 10;
                             Break_Point = 1;
                             kk = ikk;
                             printf("next2break: %d %d\n",ikk,i_contig);
                          }
                          if(n_pairs>=NUM_pairs)
                          {
                            next_kmer = id_post2[id_read[kk]];
                          }
                          kk = jj - 1;
                       }
                     }
                     printf("pairs: %d %d %d %d %d\n",i_contig,n_pairs,i_rept,next_kmer,num_454maps);
                     if(n_pairs == 0)
		     {
		       if(num_454maps == 1)
		         Break_454 = 0;
		       else if(num_454maps > 1 )
		         Check_454maps(num_454reads,num_454allhits,bulb_flag,B_joint,map_454post,map_454next,map_454path,map_454qual,max_kmerhits,max_kmernext,&Break_454,&num_454kmer); 
		       else
		         Break_454=1;
		       if(num_454pats == 1)
		       {
		         Break_454 = 0;
		         num_454kmer = num_454kmer2;
		       }
		       if(num_ppatch == 2)
		       {
		         if(abs(map_post[1] - map_post[0]) <= 3 )
		         {
                           stop_flag  = 1; 
		           Break_454=1;
			   printf("repeat stop\n");
		         }
		       }
		     }
		     else
		     {
		       Break_454 = 1;
		       Break_Point = 1;
		     }
		     if((Break_454==0)&&(num_cons<read_length))
		     {
		       next_kmer = num_454kmer;
		       Break_Point = 0;
                       printf("454hit: %d %d %d\n",i_contig,num_454maps,next_kmer);
                       stop_flag  = 0;
		     }
                     else if((Break_Point==0)&&(fuzzy_flag==0))
		     {
                       next_kmer = B_joint;
		     }
		     else
		     {
                       if((n_pairs==0)&&(fuzzy_flag==0))
                       {
		         int flag_map = 0;
   		         if(num_ppatch==2)
		         {
		           if((map_post[0]<set_offset)&&(map_post[1]>=(set_offset+1)))
		           {
                             next_kmer = map_kmer[1];
			     flag_map = 1;
 			   }
			   else if((map_post[1]<4)&&(map_post[0]>=5))
			   {
                             next_kmer = map_kmer[0];
			     flag_map = 1;
			   }
		           if(abs(map_post[1] - map_post[0]) <= 3 )
                             stop_flag  = 1; 
		         }
		         if(flag_map==0)
		         {
                           path_rate = max_move1;
                           if(max_move2==0)
                             path_rate = 2.0;
                           else
                             path_rate = path_rate/max_move2;
                             if(((path_rate>1.4)&&(max_move2<set_mismatch))||(Break_Point==0))
                             {
                             printf("rate: %d %d %f %d %d\n",next_move1,max_move2,path_rate,Break_Point,id_post[max_start]);
                             next_kmer = id_post[max_start];
                           }
                           else
                           {
                         
                             for(kk=0;kk<i_rept;kk++)
                             {
                                if(kmer_masks[next_kmer]==0)
                                {
                                  printf("next: %-2d %-25s %ld %d %d %d %d\n",kk,(seq+iseq)->name,kmer_long,next_kmer,ctgkmer_next[kk],num_cons,ikk);
                                  repeat_word[n_repeats] = next_kmer;
                                  repeat_next[n_repeats] = ctgkmer_next[kk];
                                  n_repeats++;
                                }
                             }
                             kmer_masks[next_kmer]= 1;
                             stop_flag  = 1;
                           }
                         }
                       }
                       else if(fuzzy_flag == 0)
                       {
                         Maximum_Occurr(n_pairs,pair_next,pair_index,hit_list);
                         ArraySort_Int2(n_pairs,pair_next,pair_index);
                         if(hit_list[1]==0)
                         {
                           next_kmer = pair_next[hit_list[2]];
			   /*  pay attention to this   */
			   if(n_pairs>0)
			   {
			     if((n_pairs == 1)&&(pair_pdist[0]<10))
			       stop_flag = 1;
			     else
			       stop_flag = 0;
			   }
                           printf("pair-next: %d %d %d\n",next_kmer,n_pairs,stop_flag);
                         }
                         else if(hit_list[1]>=2)
		         {
                           Patch_Occurr(n_pairs,pair_pgaps,pair_patch,pair_pnext,hit_list);
			   if(hit_list[0]>=0)
		  	     next_kmer = hit_list[0];
		  	   else
			     stop_flag  = 1;
		         }
			 else if((hit_list[1]==1)&&(n_pairs>=3))
			 {
			   if(pair_pnext[0] == pair_pnext[1])
			   {
			     if(pair_quals[n_pairs-1] < 20)
			     {
			       stop_flag = 0;
			       next_kmer = pair_pnext[0];
			     }
			     else
			       stop_flag  = 1;
			   }
                           else
                           {
                             if(pair_quals[0] < 20)
			     {
			       stop_flag = 0;
			       next_kmer = pair_pnext[1];
			     }
			     else
			       stop_flag  = 1;
			   }
                         }
		         else
                           stop_flag  = 1; 
                       }
		     }
		     if(a454_flag&&num_454negs>=2)
		     {
		       stop_flag  = 1;
		     }
                   }
                   else
                   {
                     next_kmer = id_post[max_start];
                   }
                   if(single_hit >=hit_single)
                     stop_flag  = 1;
               }
               else
               {
                 stop_flag  = 1;
                 printf("out-extend: %d\n",exten_len);
               }
             }
	     else
	       stop_flag  = 1;
             if(kmer_baseR[i_node]=='n')
             {
               kmer_baseR[i_node] = kmer_baseF[i_node];
               kmer_nodes[i_node] = next_kmer;
             }
//	   printf("kkkk: %d %d\n",i_node,next_kmer);
             if(num_cons >= read_length)
	       stop_flag  = 1;
           }
           else
           {
             stop_flag  = 1;
             printf("out-num_sect: %d %d %d %ld %s\n",num_sect,kmer_depth,kmer_depth2,kmer_long,(seq+iseq)->name);
           }
        }
	Get_Kmerbase(kmer_long,n_Sbase,kmerch,kmerrc);
//        printf("kmer: %d %d %d %ld %s %s %d\n",num_cons,i_contig,move_dir,kmer_long,kmerch,kmerrc,hit_depth);
        if(hit_depth > 1)
	{
          for(k=1;k<n_Sbase;k++)
	  {
            if((num_cons+k)<read_length)
	    {
	       if(move_dir == 1)
                 cons_reads[i_contig][num_cons-1+k] = kmerch[k];
	       else
                 cons_reads[i_contig][num_cons-1+k] = kmerrc[k];
	    }
	  }
	}
}

/*   Subroutine to walk euler path to extend read length   */
/* =============================================  */
void Euler_Path2(int *break_kmer,int b_kmer,int i_joint,int move_length,int n_depth,int n_Sbase,fasta *seq,int *patch_ofset,int *patch_index,int *patch_2kmer,long *patch_array,long *patch_head,long *patch_nlink,long *patch_plink)
/* =============================================  */
{
     int i,j,k,w,iseq=0;
     void ArraySort_Int2(int n, int *arr, int *brr);
     void ArraySort_Mix(int n, long *arr, int *brr);
     void Check_Break_pile(int nRead,int *id_post2, int n_depth, int *Break_Point);
     int Break_Point;
     int num_sect;
     fasta *seqp;
     long kmer_long;
     int kk,next_kmer,km,i_kmer,last_hit=0;
     int ref_offset = 0,ref_rcdex,stop_flag;
     int next_move1,next_move2,max_start;
     float path_rate = 1.0;
     int kshift = 2;
     int r_depth = n_depth*2.5;
     long next_link,last_link,id_nkmer[n_depth];
     char cons_base;
     int reads_last[r_depth],reads_kmer[r_depth],reads_idex[r_depth],reads_rcdx[r_depth],reads_rclt[r_depth],reads_loci[r_depth];
     int id_read[n_depth],id_post[n_depth],id_post2[n_depth],id_rddex[n_depth]; //id_read will never exceed kmer_depth
     void Maximum_Occurr(int nHit,int *hit_array, int *hit_index,int *hit_list);
     int num_lasthit;
     int seq_flag;

     if(move_length==read_length)
       seq_flag = 1;
     else
       seq_flag = 0;
//        if(num_cons>=100)
        if(num_cons>=2)
        {
          printf("read2: %d %d %d %s\n",i_joint,num_cons,i_contig,(seq+n_list[i_contig-1])->name);
        }
        printf("contig: %d %d %d %s\n",i_contig,repeat_word[i_joint],repeat_next[i_joint],(seq+iseq)->name);
        w = b_kmer;
        next_kmer = w;
        print_flag = 0;
//        n_list[i_contig] = iseq;
        num_cons = 0;
        stop_flag=0;
        memset(reads_last,0,4*r_depth);
        memset(reads_kmer,0,4*r_depth);
        memset(reads_idex,0,4*r_depth);
        memset(reads_rcdx,0,4*r_depth);
        memset(reads_rclt,0,4*r_depth);
        num_lasthit = 0;
        rcdex_pres = path_dir;

        move_dir = 1;
        while((next_kmer>=0)&&(stop_flag == 0))
        {
             long i_base;
             i = kmer_idgrp[next_kmer];
             k = kmer_locus[next_kmer];
             kmer_long = patch_array[patch_head[i]+k];
             num_sect = kmer_nhits[next_kmer];
             num_cons++;
             if(num_cons>=move_length)
               stop_flag  = 1;
             if((num_sect <= kmer_depth)&&(num_sect>=kmer_depth2))
             {
               int dupst_flag = 0;
               int exten_len = 0;
               int duped_flag = 0,id_rcdx[kmer_depth];
               int idd, idk,idd2,idk2,ref_offset2;
               int rc,st,ed,max_len,out_flag = 1,max_read;
               float rate = 0.0;
               int max_move1,max_move2,jj,ikk,rcdex1,rcdex2,rcdex_nows,num_hits;

               max_read = 0;
               for(kk=0;kk<num_sect;kk++)
               {
                  idd = kmer_index[patch_head[i]+kk+k];
                  idk = patch_index[patch_head[i]+idd];
                  seqp = seq + idk;
                  if(seqp->length>max_read)
                    max_read = seqp->length;
               }
               for(kk=0;kk<num_sect;kk++)
               {
                  idd = kmer_index[patch_head[i]+kk+k];
                  idk = patch_index[patch_head[i]+idd];
                  seqp = seq + idk;
                  memset(m_align[kk],'\0',100);
//                  memset(m_qualy[kk],0,100*4);
                  ref_offset = (int)(patch_ofset[patch_head[i]+idd]>>kshift);
                  ref_rcdex = (int)(patch_ofset[patch_head[i]+idd]&003);
                  id_read[kk] = kk;
                  id_post[kk] = ref_offset;
                  id_rcdx[kk] = ref_rcdex;
                  if(ref_rcdex == 1)
                  {
                    st = max_read - n_Sbase - ref_offset; 
                    for(rc=0;rc<st;rc++)
                       m_align[kk][rc] = ' ';
                    ed = seqp->length; 
                    for(rc=0;rc<ed;rc++)
                    {
                       m_align[kk][st+rc] = seqp->data[rc];
//                       m_qualy[kk][st+rc] = seqp->qual[rc];
/*                quality screen    */
                       if(seqp->qual[rc]<qthresh)
                         m_align[kk][st+rc] = tolower(m_align[kk][st+rc]); 
                    }
                  }
                  else
                  {
                    st = max_read - n_Sbase - ref_offset; 
                    for(rc=0;rc<st;rc++)
                       m_align[kk][rc] = ' ';
                    ed = seqp->length-1; 
                    for(rc=0;rc<seqp->length;rc++)
                    {
                       if(seqp->data[ed-rc] == 'A')
                         m_align[kk][st+rc] = 'T'; 
                       else if(seqp->data[ed-rc] == 'C')
                         m_align[kk][st+rc] = 'G'; 
                       else if(seqp->data[ed-rc] == 'G')
                         m_align[kk][st+rc] = 'C'; 
                       else if(seqp->data[ed-rc] == 'T')
                         m_align[kk][st+rc] = 'A'; 
                       else
                         m_align[kk][st+rc] = seqp->data[rc]; 
/*                quality screen    */
                       if(seqp->qual[ed-rc]<qthresh)
                         m_align[kk][st+rc] = tolower(m_align[kk][st+rc]); 
//                       m_qualy[kk][st+rc] = seqp->qual[ed-rc];
                    }
                  }
               }
               ArraySort_Int2(num_sect,id_post,id_read);
               idd = kmer_index[patch_head[i]+k+id_read[0]];
               idk = patch_index[patch_head[i]+idd];
               ref_offset = (int)(patch_ofset[patch_head[i]+idd]>>kshift);
               idd2 = kmer_index[patch_head[i]+k+id_read[num_sect-1]];
               ref_offset2 = (int)(patch_ofset[patch_head[i]+idd2]>>kshift);
//               max_len = max_read + (seq+idk)->length - n_Sbase - ref_offset; 
               max_len = max_read + max_read - n_Sbase - ref_offset; 
               exten_len = ref_offset2-ref_offset;
//               if(exten_len>=1)
               if(exten_len>=0)
               {
                 int m;

                 memset(kmer_cons,'\0',max_len+1);
                 memset(kmer_qual,0,(max_len+1)*4);
                 memset(m_pileup,0,(max_len+1)*4);
//                   printf("=================================================================================\n");
                 for(m=0;m<max_len;m++)
                 {
                   int num_hit = 0, max_hit = 0, hit[4]= {0},hit_idx;
                   hit[0] = 0;
                   hit[1] = 0;
                   hit[2] = 0;
                   hit[3] = 0;
                   hit_idx = -1;
                   memset(rbase,'\0',max_len); 
                   for(kk=0;kk<num_sect;kk++)
                   {
                      if((m_align[id_read[kk]][m]!= ' ')&&(m_align[id_read[kk]][m]!= '\0'))
                      {
                        if((m_align[id_read[kk]][m] == 'A')||(m_align[id_read[kk]][m] == 'a'))
                          hit[0]++;
                        else if((m_align[id_read[kk]][m] == 'C')||(m_align[id_read[kk]][m] == 'c'))
                          hit[1]++;
                        else if((m_align[id_read[kk]][m] == 'G')||(m_align[id_read[kk]][m] == 'g'))
                          hit[2]++;
                        else if((m_align[id_read[kk]][m] == 'T')||(m_align[id_read[kk]][m] == 't'))
                          hit[3]++;
                        rbase[num_hit] = m_align[id_read[kk]][m];
                        num_hit++;
                      }
                   }
                   for(kk=0;kk<4;kk++)
                   {
                      if(hit[kk]>max_hit)
                      {
                        max_hit = hit[kk];
                        hit_idx = kk;
                      }
                   }
                   m_pileup[m] = num_hit;
/*                   rate = max_hit;
                   rate = rate/num_hit;
                   if(rate < 0.2)
                   {
                     out_flag  = 0;
                     printf("out: %d %d %d %f %ld %d\n",max_hit,num_hit,num_sect,rate,kmer_long,m);
                     m = max_len;
                   }       */
                   if(hit_idx == 0)
                     kmer_cons[m] = 'A';
                   else if(hit_idx == 1)
                     kmer_cons[m] = 'C';
                   else if(hit_idx == 2)
                     kmer_cons[m] = 'G';
                   else if(hit_idx == 3)
                     kmer_cons[m] = 'T';
                   else
                     kmer_cons[m] = 'N';
                   kmer_qual[m] = 40; 
                 }
                   for(kk=0;kk<num_lasthit;kk++)
                   {
                      reads_kmer[kk] = reads_last[kk];
                      reads_rcdx[kk] = reads_rclt[kk];
                      reads_loci[kk] = reads_loci[kk];
                      reads_idex[kk] = kk;
                   }
                   ikk = 0;
                   for(kk=0;kk<num_sect;kk++)
                   {
                      int low_index,uph_index;
                      long kmer_next;

                      idd = kmer_index[patch_head[i]+k+id_read[kk]];
                      idk = patch_index[patch_head[i]+idd];
                      ref_rcdex = (int)(patch_ofset[patch_head[i]+idd]&003);
                      next_link = patch_nlink[patch_head[i]+idd];
                      last_link = patch_plink[patch_head[i]+idd];
                      if(ref_rcdex==1)
                      {
                          reads_kmer[num_lasthit+ikk] = idk;
                          reads_rcdx[num_lasthit+ikk] = ref_rcdex;
                          reads_idex[num_lasthit+ikk] = num_lasthit+ikk;
                          reads_last[ikk] = idk;
                          reads_rclt[ikk] = ref_rcdex;
                          ikk++;
                      }
                      else
                      {
                          reads_kmer[num_lasthit+ikk] = idk;
                          reads_rcdx[num_lasthit+ikk] = ref_rcdex;
                          reads_idex[num_lasthit+ikk] = num_lasthit+ikk;
                          reads_last[ikk] = idk;
                          reads_rclt[ikk] = ref_rcdex;
                          ikk++;
                      }
                   }
                   last_hit = ikk;
                 if(num_lasthit>0)
                 {
                   num_hits = num_lasthit+ikk;
                   ArraySort_Int2(num_lasthit+ikk,reads_kmer,reads_idex);
                   rcdex1 = 0;
                   rcdex2 = 0;
                   for(kk=0;kk<num_hits;kk++)
                   {
                      jj = kk+1;
                      while((jj<num_hits)&&(reads_kmer[kk]==reads_kmer[jj]))
                      {
                        jj++;
                      }
                      if((jj-kk)==2)
                      {
                        if(reads_rcdx[reads_idex[jj-1]]==reads_rcdx[reads_idex[kk]])
                          rcdex1++;
                        else
                          rcdex2++;
                      }
                      kk = jj - 1;
                   }
                 }
                 else
                 {
                   rcdex1 = 10;
                   rcdex2 = 0;
                 }

                 ikk = 0;
                 if(print_flag)
                 printf("============================= %d %d\n",num_cons,num_lasthit);
                 for(kk=0;kk<num_sect;kk++)
                 {
                    int low_index,uph_index,j2,jk;
                    long kmer_next;

                    idd = kmer_index[patch_head[i]+k+id_read[kk]];
                    idk = patch_index[patch_head[i]+idd];
                    ref_offset = (int)(patch_ofset[patch_head[i]+idd]>>kshift);
                    ref_rcdex = (int)(patch_ofset[patch_head[i]+idd]&003);
                    next_link = patch_nlink[patch_head[i]+idd];
                    last_link = patch_plink[patch_head[i]+idd];
                    reads_kmer[num_lasthit+ikk] = idk;
                    reads_rcdx[num_lasthit+ikk] = ref_rcdex;
                    reads_idex[num_lasthit+ikk] = num_lasthit+ikk;
                    if(print_flag)
                      printf("%d %d %-3d %-30s %-55s ",num_sect,ref_rcdex,ref_offset,(seq+idk)->name,m_align[id_read[kk]]);
/*                    if((num_pair[idk]==1)&&(num_cons>100))
                    {
                      int rd = mat_pair[idk];
                      if((rpf_contig[rd] == i_contig)&&((num_cons-rpf_offset[rd])>100)&&((num_cons-rpf_offset[rd])<1400))
                      {
                        rpf_contig[idk] = i_contig;
                        rpf_offset[idk] = num_cons;
                      }        
                    }
                    else
                    {
                      rpf_contig[idk] = i_contig;
                      rpf_offset[idk] = num_cons;
                    }           */
                    rpf_contig[idk] = i_contig;
                    rpf_offset[idk] = num_cons;
                    rcdex_nows = 0;
                    if(rcdex1>=rcdex2)
                    {
                      if(rcdex_pres==0)
                        rcdex_nows = 0;
                      else
                        rcdex_nows = 1;
                    }
                    else
                    {
                      if(rcdex_pres==0)
                        rcdex_nows = 1;
                      else
                        rcdex_nows = 0;
                    }
                    kmer_dir = rcdex_pres;
                    kmer2_dir = rcdex_nows;
                    if(rcdex_nows == 0)
                    {
                      if(ref_rcdex==1)
                      {
                        if(next_link>=0)
                        {
                          char *ss;
                          int *sq;
                          ss = m_align[id_read[kk]]+max_read;
                          sq = m_qualy[id_read[kk]]+max_read;
                          uph_index = (int)(patch_nlink[patch_head[i]+idd]>>nshift);
                          low_index = (int)(patch_nlink[patch_head[i]+idd]&hmask);
                          i_kmer = patch_2kmer[patch_head[uph_index]+low_index];
			  j2= kmer_idgrp[i_kmer];
			  jk = kmer_locus[i_kmer];
			  kmer_next = patch_array[patch_head[j2]+jk];
                         
                          neib_kmer[ikk][nkm_len] = '\0'; 
                          for(m=0;m<nkm_len;m++)
                          {
                             if((*ss!='\0')&&(*sq>=qthresh))
                               neib_kmer[ikk][m] = *(ss);
                             else
                               neib_kmer[ikk][m] = 'N';
                             ss++;
                             sq++;
                          }
                          if(print_flag)
                            printf("w1: %ld %ld %d %d %d\n",kmer_long,kmer_next,i_kmer,rcdex1,rcdex2);
                          id_nkmer[ikk] = kmer_next;
                          id_read[ikk] = ikk;
                          id_rddex[ikk] = idk;
                          id_post[ikk] = i_kmer;
                          id_post2[ikk] = i_kmer;
                          move_dir = 1;
                          ikk++;
                        }
                        else
                          if(print_flag) printf("\n");
                      }
                      else
                      {
                        if(last_link>=0)
                        {
                          char *ss;
                          int *sq;
                          ss = m_align[id_read[kk]]+max_read;
                          sq = m_qualy[id_read[kk]]+max_read;
                          uph_index = (int)(patch_plink[patch_head[i]+idd]>>nshift);
                          low_index = (int)(patch_plink[patch_head[i]+idd]&hmask);
                          i_kmer = patch_2kmer[patch_head[uph_index]+low_index];
			  j2= kmer_idgrp[i_kmer];
			  jk = kmer_locus[i_kmer];
			  kmer_next = patch_array[patch_head[j2]+jk];

                          neib_kmer[ikk][nkm_len] = '\0'; 
                          for(m=0;m<nkm_len;m++)
                          {
                             if((*ss!='\0')&&(*sq>=qthresh))
                               neib_kmer[ikk][m] = *(ss);
                             else
                               neib_kmer[ikk][m] = 'N';
                             ss++;
                             sq++;
                          }
                          if(print_flag)
                            printf("w2: %ld %ld %d %d %d\n",kmer_long,kmer_next,i_kmer,rcdex1,rcdex2);
                          id_read[ikk] = ikk;
                          id_post[ikk] = i_kmer;
                          id_post2[ikk] = i_kmer;
                          id_rddex[ikk] = idk;
                          id_nkmer[ikk] = kmer_next;
                          move_dir = 1;
                          ikk++;
                        }
                        else
                          if(print_flag) printf("\n");
                      }
                    }
                    else
                    {
                      if(ref_rcdex==2)
                      {
                        if(next_link>=0)
                        {
                          char *ss;
                          int *sq;
                          ss = m_align[id_read[kk]]+max_read;
                          sq = m_qualy[id_read[kk]]+max_read;
                          uph_index = (int)(patch_nlink[patch_head[i]+idd]>>nshift);
                          low_index = (int)(patch_nlink[patch_head[i]+idd]&hmask);
                          i_kmer = patch_2kmer[patch_head[uph_index]+low_index];
			  j2= kmer_idgrp[i_kmer];
			  jk = kmer_locus[i_kmer];
			  kmer_next = patch_array[patch_head[j2]+jk];
                         
                          neib_kmer[ikk][nkm_len] = '\0'; 
                          for(m=0;m<nkm_len;m++)
                          {
                             if((*ss!='\0')&&(*sq>=qthresh))
                               neib_kmer[ikk][m] = *(ss);
                             else
                               neib_kmer[ikk][m] = 'N';
                             ss++;
                             sq++;
                          }
                          if(print_flag)
                            printf("w3: %ld %ld %d %d %d\n",kmer_long,kmer_next,i_kmer,rcdex1,rcdex2);
                          id_read[ikk] = ikk;
                          id_post[ikk] = i_kmer;
                          id_post2[ikk] = i_kmer;
                          id_rddex[ikk] = idk;
                          id_nkmer[ikk] = kmer_next;
                          move_dir = 2;
                          ikk++;
                        }
                        else
                          if(print_flag) printf("\n");
                      }
                      else
                      {
                        if(last_link>=0)
                        {
                          char *ss;
                          int *sq;
                          ss = m_align[id_read[kk]]+max_read;
                          sq = m_qualy[id_read[kk]]+max_read;
                          uph_index = (int)(patch_plink[patch_head[i]+idd]>>nshift);
                          low_index = (int)(patch_plink[patch_head[i]+idd]&hmask);
                          i_kmer = patch_2kmer[patch_head[uph_index]+low_index];
			  j2= kmer_idgrp[i_kmer];
			  jk = kmer_locus[i_kmer];
			  kmer_next = patch_array[patch_head[j2]+jk];
                         
                          neib_kmer[ikk][nkm_len] = '\0'; 
                          for(m=0;m<nkm_len;m++)
                          {
                             if((*ss!='\0')&&(*sq>=qthresh))
                               neib_kmer[ikk][m] = *(ss);
                             else
                               neib_kmer[ikk][m] = 'N';
                             ss++;
                             sq++;
                          }
                          if(print_flag)
                            printf("w4: %ld %ld %d %d %d\n",kmer_long,kmer_next,i_kmer,rcdex1,rcdex2);
                          id_read[ikk] = ikk;
                          id_post[ikk] = i_kmer;
                          id_post2[ikk] = i_kmer;
                          id_rddex[ikk] = idk;
                          id_nkmer[ikk] = kmer_next;
                          move_dir = 2;
                          ikk++;
                        }
                        else
                          if(print_flag) printf("\n");
                      }
                    }
                 }
                 if(print_flag)
                   printf("---------------------------------------------------------------------------------\n");
                   if(move_dir==1)
                   {
                     i_base = kmer_long >> (2*(n_Sbase-1));
                     if(i_base==0)
                       cons_base = 'A';
                     else if(i_base==1)
                       cons_base = 'C';
                     else if(i_base==2)
                       cons_base = 'G';
                     else if(i_base==3)
                       cons_base = 'T';
                     else
                       cons_base = 'C';
                   }
                   else
                   {
                     i_base = kmer_long&003;
                     if(i_base==0)
                       cons_base = 'T';
                     else if(i_base==1)
                       cons_base = 'G';
                     else if(i_base==2)
                       cons_base = 'C';
                     else if(i_base==3)
                       cons_base = 'A';
                     else
                       cons_base = 'C';
                   }
                 if(print_flag)
                   printf("%-7d %d %s  %-49s %ld %c\n",num_cons,i_contig,(seq+n_list[i_contig])->name,kmer_cons,ikk,cons_base);

                   path_kmer[num_cons-1] = next_kmer;
                   path_idir[num_cons-1] = rcdex_pres;
                   path_cons[num_cons-1] = cons_base;
                   path_cover[num_cons-1] = ikk;
                   rcdex_pres = rcdex_nows;
                   ArraySort_Int2(ikk,id_post,id_read);
                   max_move1 = 0;
                   max_move2 = 0;
                   max_start = 0;
                   for(kk=0;kk<ikk;kk++)
                   {
                      jj = kk+1;
                      while((jj<num_sect)&&(id_post[kk]==id_post[jj]))
                      {
                        jj++;
                      }
                      if((jj-kk)>max_move1)
                      {
                        if(max_move1 > 0)
                          max_move2 = max_move1;
                        max_move1 = jj-kk;
                        max_start = kk;
                      }
                      if(((jj-kk)>max_move2)&&(max_start!=kk))
                        max_move2 = jj-kk;
                 if(print_flag)
                      printf("max: %d %d %d %d %d %d\n",kk,id_post[kk],i_contig,max_move1,max_move2,max_start);
                      kk = jj - 1;
                   }
                   if((max_move1<ikk)&&(max_move2==0))
                     max_move2 = ikk - max_move1;
                   num_lasthit = last_hit;

                   for(kk=0;kk<ikk;kk++)
                      id_read[kk] = kk;
                   ArraySort_Mix(ikk,id_nkmer,id_read);
                   if(max_move2>=2)
                   {
                     Check_Break_pile(ikk,id_post2,n_depth,&Break_Point);
                     if(Break_Point==1)
                       stop_flag = 1;
                   }

                   next_move1 = 0;
                   next_move2 = 0;
                   for(kk=0;kk<ikk;kk++)
                   {
                      jj = kk+1;
                      while((jj<ikk)&&(id_nkmer[kk]==id_nkmer[jj]))
                      {
                        jj++;
                      }
                      if((jj-kk)>next_move1)
                      {
                        if(next_move1>0)
                          next_move2 = next_move1;
                        next_move1 = jj-kk;
                      }
                      kk = jj - 1;
                   }
		   if((ikk==0)||(max_move1==0))
		   {
		     ikk = 1;
		     max_move1 = 1;
		     stop_flag = 1;
		   }
                   path_rate = next_move1;
                   path_rate = path_rate/ikk;
                   rate  = max_move2;
                   rate  = rate/max_move1;
                   if((path_rate < 0.65)||(rate > 0.52)||(max_move2>=set_mismatch))
                   {
                     int m,n_pairs = 0;
                     int pair_flag  = 0;
                     int NUM_pairs=2;
                       printf("repeat: %d %f %d %s\n",ikk,path_rate,iseq,(seq+iseq)->name);
                     for(kk=0;kk<ikk;kk++)
                     {
                        n_pairs = 0;
                        jj = kk+1;
                        while((jj<ikk)&&(id_nkmer[kk]==id_nkmer[jj]))
                        {
                          jj++;
                        }
                        for(m=kk;m<jj;m++)
                        {
                           int pdk = id_read[m];
                           int rd = mat_pair[id_rddex[pdk]];
//                           printf("pass: %-2d %-25s %-42s %ld %d %d\n",m,(seq+id_rddex[pdk])->name,m_align[pdk],kmer_long,id_post2[id_read[m]],max_move2);
                           if((rpf_contig[rd] == i_contig)&&((num_cons-rpf_offset[rd])>50)&&((num_cons-rpf_offset[rd])<1400))
                             n_pairs++; 
                        }
                        if((jj-kk)>=2)
                        {
                           int pdk = id_read[kk];
                           int rd = mat_pair[id_rddex[pdk]];
                           int idgrp = kmer_idgrp[id_post2[id_read[kk]]];
                           int locus = kmer_locus[id_post2[id_read[kk]]];
                           int iidd =  kmer_index[patch_head[idgrp]+locus];
                           int iidk = patch_index[patch_head[idgrp]+iidd];
//                           if(num_cons<100)
//                             printf("small: %s %d\n",path_cons,max_move2);
                        }
                        if(n_pairs>=NUM_pairs)
                          next_kmer = id_post2[id_read[kk]];
                        kk = jj - 1;
                     }
//                     if(n_pairs==0)
                     {
                       path_rate = max_move1;
                       if(max_move2==0)
                         path_rate = 2.0;
                       else
                         path_rate = path_rate/max_move2;
                       if((path_rate>1.4)&&(max_move2<set_mismatch))
                       {
                         printf("rate: %d %d %f %d\n",max_move1,max_move2,path_rate,ikk);
                       }
                       else
                       {
                         kmer_masks[next_kmer]= 1;
                         stop_flag  = 1;
                       }
                     }
//                     else
//                       printf("rate: %d %d %f %d\n",next_move1,NUM_pairs,path_rate,ikk);
              
                   }
                   else
                   {
                   }
                   *break_kmer = next_kmer;
                   next_kmer = id_post[max_start];
                 if((print_flag)&&(num_cons<2))
                      printf("yyy: %d %d %d\n",next_kmer,n_node,repeat_word[n_node]);
                   if((num_cons<4)&&(next_kmer == repeat_word[n_node]))
                   {
                     stop_flag  = 1;
                     *break_kmer = repeat_word[n_node];
                   }
//                     printf("rate: %d %d %f %d\n",next_move1,next_move2,path_rate,ikk);
               }
               else
               {
                 stop_flag  = 1;
                 printf("out-extend: %d\n",exten_len);
               }
               
             }
             else
             {
               stop_flag  = 1;
               printf("out-num_sect: %d %d %d %ld %s\n",num_sect,kmer_depth,kmer_depth2,kmer_long,(seq+iseq)->name);
             }
        }
}


/*   Subroutine to the DNA sequences close to the break point   */
/* =============================================  */
void Check_Break(int nRead,int *id_post2, int n_depth, int *Break_Point, int *B_joint)
/* =============================================  */
{
     int i,j,k,nseq;
     int seqc, nb_array[n_depth],nb_index[n_depth];
     int hitBase[4],hitIndex[4];
     int next_move1,next_move2,stop_flag,n_hits,n_diffs;
     float path_rate = 0.0,path_rate2;
     void ArraySort_Int2(int n, int *arr, int *brr);
     void ArraySort2_Int2(int n, int *arr, int *brr);
     char neib_seq[n_depth];
     int i_idpost = 0;
     int i_hits;

     *Break_Point = 1;
     *B_joint = 0;
     n_diffs = 0;
     for(i=0;i<1;i++)
     {
        memset(neib_seq,'\0',n_depth);
        memset(hitBase,0,16);
        for(j=0;j<nRead;j++)
        {
           nb_index[j] = j;
           nb_array[j] = 0;
           neib_seq[j] = neib_kmer[j][i];
           if(neib_kmer[j][i] == 'A')
	   {
             nb_array[j] = 1;
	     hitIndex[0] = j;
	     hitBase[0]++;
	   }
           else if(neib_kmer[j][i] == 'C')
	   {
             nb_array[j] = 2;
	     hitIndex[1] = j;
	     hitBase[1]++;
	   }
           else if(neib_kmer[j][i] == 'G')
	   {
             nb_array[j] = 3;
	     hitIndex[2] = j;
	     hitBase[2]++;
	   }
           else if(neib_kmer[j][i] == 'T')
	   {
             nb_array[j] = 4;
	     hitIndex[3] = j;
	     hitBase[3]++;
	   }
	   if((i==0)&&(kmer2_hits[j]<=1))
             nb_array[j] = 0;
        }
	n_hits = 0;
        next_move1 = 0;
        next_move2 = 0;
	i_idpost = 0;
	i_hits = 0;
        for(j=0;j<4;j++)
        {
	   if(hitBase[j]>0)
	   {
	     i_idpost = hitIndex[j];
	     i_hits = j;
             n_hits++;
	   }
        }
	if(n_hits==1)
	{
	  next_move1 = hitBase[i_hits];
	  next_move2 = 0;
	}
        else
	{
	  int hit_index[4];
	  for(j=0;j<4;j++)
	     hit_index[j] = j;
          ArraySort2_Int2(4,hitBase,hit_index);
	  next_move1 = hitBase[0];
	  next_move2 = hitBase[1];
	  i_idpost = hitIndex[hit_index[0]]; 
	}
        if(sense_flag==1)
        {
          if((i==0)&&(next_move1>0)&&(next_move2==0))
          {
            *Break_Point = 0;
            i = nkm_len;
            printf("Go: %d %d\n",next_move1,next_move2);
          }
        }
        else
        {
          float rate = next_move2;
          rate = rate/next_move1;
          if(((next_move1>0)&&(rate<0.1))||((next_move1>3)&&(next_move2==1)))
          {
            *Break_Point = 0;
            i = nkm_len;
            printf("Go-1: %d %d %f\n",next_move1,next_move2,rate);
          }
        }
        if((i==0)&&(n_hits==1))
        {
          *Break_Point = 0;
          i = nkm_len;
        }
        path_rate = next_move1;
        path_rate2 = next_move2;
        if(next_move2==0)
          path_rate = 2.0;
        else
          path_rate = path_rate/nRead;
        path_rate2 = path_rate2/nRead;
	if(n_hits==0)
	  i_idpost = 0;
	if(i_idpost < 0) 
	  i_idpost = 0;
        *B_joint = id_post2[i_idpost];
        printf("NKM: %d %d %f %f %s %d %d\n",i,nRead,path_rate,path_rate2,neib_seq,next_move1,next_move2);
        if((path_rate>1.4)&&(next_move2<3))
          printf("rate-c: %d %d %d %f %d\n",i,next_move1,next_move2,path_rate,nRead);
        else
        {
          stop_flag  = 1;
        }
     }
     if(*Break_Point == 0)
     {
       int equal_flag = 0;
       if(next_move1==0)
       {
          n_hits = 0;
          for(k=0;k<nRead;k++)
	  {
             if(neib_seq[k] != 'N')
	     {
	       nb_index[n_hits] = k; 
	       n_hits++;
	     }
	  }
	  if(n_hits == 1)
	    i_idpost = nb_index[0];
	  else
	  {
	    if(kmer2_hits[nb_index[0]]==kmer2_hits[nb_index[1]])
	      equal_flag = 1;
	    else if(kmer2_hits[nb_index[0]]>kmer2_hits[nb_index[1]])
              i_idpost = nb_index[0];
	    else
              i_idpost = nb_index[1];
	  }
       }
       if(equal_flag==1)
         *Break_Point = 1;
       else
       {
         *B_joint = id_post2[i_idpost];
         printf("Go: %d %d %d %d\n",next_move1,next_move2,i_idpost,id_post2[i_idpost]);
       }
     }
}

/*   Subroutine to traverse short path to remove base errors   */
/* =============================================  */
void Short_Path(int p_node,int id_node,int kmer_len,int *sm,int *sm_list,int *short_flag,int *path_len,int *kmer_last)
/* =============================================  */
{
     int i,j,k,i_node;
     int walk_len = 0,break_flag;

     *path_len = 0;
     *short_flag = 0;
     kmer_last[0] = 0;
     kmer_last[1] = 0;
     kmer_last[2] = 0;
     i_node = id_node;
     break_flag = 0;
     while((walk_len<40)&&(sm_list[i_node]>0)&&(break_flag==0))
     {
       if(debug_flag)
         printf("loop: %d %d %d %d\n",walk_len,i_node,sm_list[i_node],id_node);
       if(walk_len == (kmer_len-1))
         kmer_last[1] = i_node;
       if(walk_len == (kmer_len-2))
         kmer_last[2] = i_node;
       if(walk_len == kmer_len)
         kmer_last[0] = i_node;
       if(sm_list[i_node]==0)
       {
         *short_flag = 1;
	 break_flag = 1;
       }
       else if(sm_list[i_node]==1)
       {
         if(p_node==sm[sm_head[i_node]])
         {
           *short_flag = 1;
	   break_flag = 1;
         }
         else
         {
           p_node = i_node;
           i_node = sm[sm_head[i_node]];
         }
       }
       else if(sm_list[i_node]<=2)
       {
         for(j=0;j<sm_list[i_node];j++)
         {
            int idd = sm[sm_head[i_node]+j];
            if(idd!=p_node)
            {
              p_node = i_node;
              i_node = idd;
              j = sm_list[i_node];
            }
         }
       }
       else
       {
         walk_len = 40;
       }
       walk_len++; 
     }
     if((walk_len>=15)&&(walk_len<60))
       *short_flag = 0;
     else
       *short_flag = 1;
     *path_len = walk_len;
}

/*   Subroutine to the DNA sequences close to the break point   */
/* =============================================  */
void Check_Break_pile(int nRead,int *id_post2, int n_depth, int *Break_Point)
/* =============================================  */
{
     int i,j,k,nseq;
     int seqc, nb_array[n_depth],nb_index[n_depth];
     int next_move1,next_move2,stop_flag,idex,idex1,break_flag;
     void ArraySort_Int2(int n, int *arr, int *brr);
     char neib_seq[n_depth];

        *Break_Point = 0;
        i = 0;
	idex = 0;
	idex1 = 0;
        memset(neib_seq,'\0',n_depth);
        for(j=0;j<nRead;j++)
        {
           nb_index[j] = j;
           nb_array[j] = 0;
           neib_seq[j] = neib_kmer[j][i];
           if(neib_kmer[j][i] == 'A')
             nb_array[j] = 1;
           else if(neib_kmer[j][i] == 'C')
             nb_array[j] = 2;
           else if(neib_kmer[j][i] == 'G')
             nb_array[j] = 3;
           else if(neib_kmer[j][i] == 'T')
             nb_array[j] = 4;
        }
        ArraySort_Int2(nRead,nb_array,nb_index);
        next_move1 = 0;
        next_move2 = 0;
        for(j=0;j<nRead;j++)
        {
           k = j+1;
           while((k<nRead)&&(nb_array[j]>0)&&(nb_array[k]==nb_array[j]))
           {
              k++;
           }
           if((k-j)>=2)
           {
             if((k-j)>=next_move1)
             {
               if(next_move1>0)
	       {
                 next_move2 = next_move1;
		 idex = idex1;
	       }
               next_move1 = k-j;
	       idex1 = nb_index[j];
             }
             else
             {
               if((k-j)>=next_move2)
	       {
                 next_move2 = k-j;
		 idex = nb_index[j];
	       }
             }
           }
           j = k - 1;
        }
	break_flag = 0;
        if((next_move2==1)||(next_move2 == 2))
	{
	  if(kmer2_hits[idex] >= 4 )
	    break_flag = 1;
	}
        if((next_move2 >= 3 )|| (break_flag))
        {
          printf("pile: %d %s %d %d %d\n",nRead,neib_seq,next_move1,next_move2,kmer2_hits[idex]);
          *Break_Point = 1;
        }
/*	if(next_move2>0)
	{
          for(j=0;j<nRead;j++)
	     printf("next-hit: %d %d %c\n",j,kmer2_hits[j],neib_seq[j]);
	}   */
}

/*   Subroutine to the DNA sequences close to the break point   */
/* =============================================  */
void Check_Break_Fuzzy(int nRead,int nFuzzy,int *id_post2, int n_depth, char *fuzzy_ray,int *Break_Point,int *got_next)
/* =============================================  */
{
     int i,j,k,nseq,p_hit,i_hit;
     int seqc, nb_array[n_depth],nb_index[n_depth];
     int next_move1,next_move2,stop_flag,n_hits,n_diffs;
     float path_rate = 0.0,path_rate2;
     void ArraySort_Int2(int n, int *arr, int *brr);
     char neib_seq[n_depth],neib_seq2[n_depth];

     *Break_Point = 0;
     n_diffs = 0;
     i = 0;
     memset(neib_seq,'\0',n_depth);
     for(j=0;j<nRead;j++)
        neib_seq2[j] = neib1_kmer[j];
     for(j=0;j<nFuzzy;j++)
     {
        nb_index[j] = j;
        nb_array[j] = 0;
        neib_seq[j] = fuzzy_ray[j];
        if(fuzzy_ray[j] == 'A')
          nb_array[j] = 1;
        else if(fuzzy_ray[j] == 'C')
          nb_array[j] = 2;
        else if(fuzzy_ray[j] == 'G')
	{
          nb_array[j] = 3;
	}
        else if(fuzzy_ray[j] == 'T')
          nb_array[j] = 4;
     }
     ArraySort_Int2(nFuzzy,nb_array,nb_index);
     next_move1 = 0;
     next_move2 = 0;
     n_hits = 0;
     for(j=0;j<nFuzzy;j++)
     {
        if(nb_array[j]>0)
	{
	  i_hit = nb_array[j];
          n_hits++;
	}
     }
     for(j=0;j<nFuzzy;j++)
     {
        k = j+1;
        while((k<nFuzzy)&&(nb_array[j]>0)&&(nb_array[k]==nb_array[j]))
        {
           k++;
        }
        if((k-j)>=2)
        {
          if((k-j)>=next_move1)
          {
            if(next_move1>0)
              next_move2 = next_move1;
            next_move1 = k-j;
	    i_hit = nb_array[j];
          }
          else
          {
            if((k-j)>=next_move2)
              next_move2 = k-j;
          }
        }
        j = k - 1;
     }
     if(n_hits == 1)
     {
       next_move1 = 1;
       next_move2 = 0;
     }
     if((next_move2 > 0)||(n_hits == 0 ))
     {
       *Break_Point = 1;
     }
     else
     {
       printf("fpile1: %d %d %s\n",n_hits,i_hit,neib_seq2);
       for(j=0;j<nRead;j++)
       {
       printf("fpile2: %d %d %c\n",n_hits,i_hit,neib1_kmer[j]);
          if((neib1_kmer[j] == 'A')||(neib1_kmer[j] == 'a'))
	    p_hit = 1;
          else if((neib1_kmer[j] == 'C')||(neib1_kmer[j] == 'c'))
	    p_hit = 2;
          else if((neib1_kmer[j] == 'G')||(neib1_kmer[j] == 'g'))
	    p_hit = 3;
          else if((neib1_kmer[j] == 'T')||(neib1_kmer[j] == 't'))
	    p_hit = 4;
	  else
	    p_hit = 0;
          if(p_hit == i_hit)
	  {
            *Break_Point = 0;
	    *got_next = id_post2[j];
            printf("fpile3: %d %c %d %d %s\n",j,neib1_kmer[j],id_post2[j],i_hit,neib_seq);
	    j = nRead;
	  }
       }
     }
}

/*   Subroutine to select a path using 454 reads  */
/* =============================================  */
void Check_454maps(int num_454reads,int nRead,int bulb_flag,int B_joint,int *map_454post,int *map_454next,int *map_454path, int *map_454qual, int max_kmerhits,int max_kmernext,int *Break_454,int *map_kmernext)
/* =============================================  */
{
     int i,k,j,num_hits,map_maxhits[nRead],map_next[nRead],map_number[nRead];
     int max_post,max_number,i_numhit,i_poshit,ave_qual;
     int map_index1[nRead],map_index2[nRead],map_quality[nRead];
     void ArraySort2_Int2(int n, int *arr, int *brr);
     float rate;

     i = 0;
     *Break_454 = 0;
     num_hits = 0;
//     for(k=0;k<nRead;k++)
//        printf("www: %d %d %d %d %d\n",k,map_454post[k],map_454next[k],map_454path[k],map_454qual[k]); 
     for(k=0;k<nRead;k++)
     {
        j = k+1;
        while((j<nRead)&&(map_454next[k]==map_454next[j]))
        {
           j++;
        }
	if((j-k)>=2)
	{
	  max_post = 0;
	  ave_qual = 0;
	  for(i=k;i<j;i++)
	  {
	     if(map_454post[i]>max_post)
	       max_post = map_454post[i];
	     ave_qual = ave_qual + map_454qual[i];
	  }
	  map_maxhits[num_hits] = max_post;
	  map_next[num_hits] = map_454next[k];
	  map_number[num_hits] = j-k;
	  map_quality[num_hits] = ave_qual/(j-k);
	  map_index1[num_hits] = num_hits;
	  map_index2[num_hits] = num_hits;
	}
	else
	{
	  map_maxhits[num_hits] = map_454post[k];
	  map_next[num_hits] = map_454next[k];
	  map_number[num_hits] = 1;
	  map_quality[num_hits] = map_454qual[k];
	  map_index1[num_hits] = num_hits;
	  map_index2[num_hits] = num_hits;
	}
        num_hits++;
        k = j - 1;
     }

     ArraySort2_Int2(num_hits,map_number,map_index1);
     for(i=0;i<num_hits;i++)
        printf("map: %d %d %d %d\n",i,map_number[i],map_maxhits[map_index1[i]],map_quality[map_index1[i]]);

     if((num_hits==2)&&(map_number[1]==1)&&(map_number[0]>=3)&&(map_quality[map_index1[1]]<10))
     {
       *map_kmernext = map_next[map_index1[0]];
       *Break_454 = 0;
     }
     else
     {
       ArraySort2_Int2(num_hits,map_maxhits,map_index2);
       i_numhit = map_index1[0];
       i_poshit = map_index2[0];
       rate = map_maxhits[1];
       rate = rate/map_maxhits[0];
          printf("abi: %d %d %f\n",map_maxhits[1],map_maxhits[0],rate);
       if(rate<set_rate)
       {
         *map_kmernext = map_next[i_poshit];
         *Break_454 = 0;
//         if(*map_kmernext==max_kmernext)
          printf("abi: %d %d\n",*Break_454,*map_kmernext);
       }
       else if((rate>0.8)&&(bulb_flag))
       {
         *map_kmernext = B_joint;
         *Break_454 = 0;
       }
       else
         *Break_454 = 1;
     }
     if((*Break_454==0)&&(*map_kmernext!=max_kmernext))
       *Break_454 = 0;
}

/*   Subroutine to the DNA sequences close to the break point   */
/* =============================================  */
void Check_Break_c3(int nRead,int *id_post2, int n_depth, int *Break_Point)
/* =============================================  */
{
     int i,j,num_low;

     i = 0;
     *Break_Point = 0;
     num_low = 0;
     for(j=0;j<nRead;j++)
     {
        if(neib_kmer[j][i] == 'N')
	  num_low++;
     }
     if(num_low==0)
       *Break_Point = 1;
}

/*   Subroutine to the number of mismatches within kmer words  */
/* =========================================*/
void Mismatch_Kmer(long kmer1, long kmer2,int n_Sbase,long *mmask,int *num_mismatch,int *mismatch_loci)
/* =========================================*/
{
     int hit_num,k;
     long DifIntbase,mask2;

     DifIntbase = kmer1^kmer2;
     hit_num = 0;
     for(k=0;k<(2*n_Sbase);k++)
     {
        if((mmask[k]&DifIntbase) > 0)
        {
          if(mask2 > 0)
          {
            if(k%2 == 0)
	    {
	      mismatch_loci[hit_num] = (k/2) + 1; 
              hit_num++;
	    }
          }
          else
	  {
	    mismatch_loci[hit_num] = (k/2) + 1; 
            hit_num++;
	  }
        }
        mask2 = mmask[k]&DifIntbase;
     }
     *num_mismatch = hit_num;
}

/*   Subroutine to calculate mmask array  */
/* =========================================*/
void Mmask_array(long *mmask,int n_Sbase)
/* =========================================*/
{
     int k;
     long long_kmer;

     long_kmer = 2L;
     mmask[0] = 1L;
     for(k=1;k<(2*n_Sbase);k++)
        mmask[k] = long_kmer<<(k-1);  
}


/*   Subroutine to sort the maximum and second maximum occurrences   */
/* =============================================  */
void Maximum_Occurr(int nHit,int *hit_array, int *hit_index,int *hit_list)
/* =============================================  */
{
     int i,j,k;
     void ArraySort_Int2(int n, int *arr, int *brr);
     int max_move1,max_move2,max_start,num_hits;

     ArraySort_Int2(nHit,hit_array,hit_index);
     max_move1 = 0;
     max_move2 = 0;
     max_start = 0;
     num_hits = 0;
     for(k=0;k<nHit;k++)
     {
        j = k+1;
        while((j<nHit)&&(hit_array[k]==hit_array[j]))
        {
           j++;
        }
        if((j-k)>=2)
          num_hits++;
        if((j-k)>max_move1)
        {
          if(max_move1 > 0)
            max_move2 = max_move1;
          max_move1 = j-k;
          max_start = k;
        }
        if(((j-k)>max_move2)&&(max_start!=k))
          max_move2 = j-k;
        k = j - 1;
     }
     hit_list[0] = max_move1;
     hit_list[1] = max_move2;
     hit_list[2] = max_start;
     hit_list[3] = num_hits;
}

/*   Subroutine to calulate insert gap size distributions   */
/* =============================================  */
void Patch_Occurr(int nHit,int *hit_array, int *hit_index,int *hit_next,int *hit_list)
/* ===========================================  */
{
     int i,j,k;
     void ArraySort_Int2(int n, int *arr, int *brr);
     int max_move1,max_move2,max_start,num_hits,sum_score;
     int m_score[nHit],m_index[nHit],m_pairs[nHit];

     max_move1 = 5000000;
     max_move2 = 5000000;
     max_start = 0;
     num_hits = 0;
     hit_list[0] = -1;
     for(k=0;k<nHit;k++)
     {
        j = k+1;
        sum_score = hit_array[k];
        while((j<nHit)&&(hit_index[k]==hit_index[j])&&(hit_next[k]==hit_next[j]))
        {
	   sum_score = sum_score + hit_array[j];
     printf("score: %d %d %d %d\n",k,j,sum_score,hit_array[j]);
           j++;
        }
        if((j-k)>=2)
	{
	  m_score[num_hits] = sum_score/(j-k);
	  m_index[num_hits] = k;
	  m_pairs[num_hits] = j-k;
     printf("sum0: %d %d %d\n",num_hits,sum_score,m_score[num_hits]);
          num_hits++;
	}
        k = j - 1;
     }
     ArraySort_Int2(num_hits,m_score,m_index);
     if(num_hits == 1)
       hit_list[0] = hit_next[0];
     else if(num_hits > 1)
     {
       if(m_score[0] < m_score[1])
         hit_list[0] = hit_next[m_index[0]];
       if(m_pairs[m_index[1]] >= 2 )
         hit_list[0] = -1;
     }
     printf("sum: %d %d %d %d\n",m_score[0],m_score[1],m_index[0],hit_list[0]);
}


/* =============================================  */
void Read_Align(fasta *seq, int *kmer_index, int *patch_index, char **m_align,int **m_qualy, int *patch_ofset,long *patch_head,int *id_read, int *id_post, int *id_rcdx,int num_sect,int kshift,int max_read,int n_Sbase,int n_depth,int i,int k)
/* =============================================  */
{
     int j,kk;
     fasta *seqp;
 
     for(kk=0;kk<num_sect;kk++)
     {
        int ref_offset,ref_rcdex,st,ed,rc;
        int idd = kmer_index[patch_head[i]+kk+k];
        int idk = patch_index[patch_head[i]+idd];
        seqp = seq + idk;
        memset(m_align[kk],'\0',100);
        memset(m_qualy[kk],0,400);
        ref_offset = (int)(patch_ofset[patch_head[i]+idd]>>kshift);
        ref_rcdex = (int)(patch_ofset[patch_head[i]+idd]&003);
        id_read[kk] = kk;
        id_post[kk] = ref_offset;
        id_rcdx[kk] = ref_rcdex;
        if(ref_rcdex == 1)
        {
         st = max_read - n_Sbase - ref_offset;
	 if(st>=0)
	 {
          for(rc=0;rc<st;rc++)
             m_align[kk][rc] = ' ';
          ed = seqp->length;
	  if(ed>=50)
	    ed = 54-st;
          for(rc=0;rc<ed;rc++)
          {
             m_align[kk][st+rc] = seqp->data[rc];
             m_qualy[kk][st+rc] = seqp->qual[rc];
/*           quality screen    */
             if(seqp->qual[rc]<qthresh)
                m_align[kk][st+rc] = tolower(m_align[kk][st+rc]);
          }
	 }
	 else
	 {
	  st = ref_offset - max_read + n_Sbase;
	  ed =  max_read + 1;
          for(rc=0;rc<ed;rc++)
          {
	     if((st+rc) < seqp->length)
	     {
               m_align[kk][rc] = seqp->data[st+rc];
               m_qualy[kk][rc] = seqp->qual[st+rc];
/*             quality screen    */
               if(seqp->qual[st+rc]<qthresh)
                 m_align[kk][rc] = tolower(m_align[kk][rc]);
	     }
          }
	 }
        }
        else
        {
	 int st2 = 0;
         st = max_read - n_Sbase - ref_offset;
	 if(st>=0)
	 {
          for(rc=0;rc<st;rc++)
             m_align[kk][rc] = ' ';
          ed = seqp->length-1;
	  if(ed>=50)
	  {
	    ed = max_read-st;
	    st2 = seqp->length - ref_offset - n_Sbase -1;
	  }
          for(rc=0;rc<=ed;rc++)
          {
             if(seqp->data[ed+st2-rc] == 'A')
               m_align[kk][st+rc] = 'T';
             else if(seqp->data[ed+st2-rc] == 'C')
               m_align[kk][st+rc] = 'G';
             else if(seqp->data[ed+st2-rc] == 'G')
               m_align[kk][st+rc] = 'C';
             else if(seqp->data[ed+st2-rc] == 'T')
               m_align[kk][st+rc] = 'A';
             else
               m_align[kk][st+rc] = seqp->data[rc];
/*           quality screen    */
             if(seqp->qual[ed+st2-rc]<qthresh)
               m_align[kk][st+rc] = tolower(m_align[kk][st+rc]);
             m_qualy[kk][st+rc] = seqp->qual[ed+st2-rc];
          }
	 }
	 else
	 {
	  st = seqp->length - ref_offset + max_read - n_Sbase -1;
	  ed =  max_read+3;
	  if((st-ed) < -1)
	    ed = st+1;
          for(rc=0;rc<ed;rc++)
          {
             if(seqp->data[st-rc] == 'A')
               m_align[kk][rc] = 'T';
             else if(seqp->data[st-rc] == 'C')
               m_align[kk][rc] = 'G';
             else if(seqp->data[st-rc] == 'G')
               m_align[kk][rc] = 'C';
             else if(seqp->data[st-rc] == 'T')
               m_align[kk][rc] = 'A';
             else
               m_align[kk][rc] = seqp->data[st+rc];
/*           quality screen    */
             if(seqp->qual[st-rc]<qthresh)
               m_align[kk][rc] = tolower(m_align[kk][rc]);
             m_qualy[kk][rc] = seqp->qual[st-rc];
	  }
	 }
        }
     }
}

/*   Subroutine to buildup multiple read alignment   */
/* =============================================  */
void Read_AlignFuzzy(fasta *seq, int *kmer_index, int *patch_index, char **m_align,char *fuzzy_ray,int **m_qualy, int *patch_ofset,long *patch_head,int *id_read, int *id_post, int *id_rcdx,int num_sect,int kshift,int max_read,int n_Sbase,int n_depth,int i,int k,int num_sum,int num_mismatch,int *mismatch_loci,long kmer_long)
/* =============================================  */
{
     int j,kk,pos;
     fasta *seqp;
     char NC;

     for(kk=0;kk<num_sect;kk++)
     {
        int ref_offset,ref_rcdex,st,ed,rc;
        int idd = kmer_index[patch_head[i]+kk+k];
        int idk = patch_index[patch_head[i]+idd];
        seqp = seq + idk;
        memset(m_align[kk+num_sum],'\0',n_depth);
        memset(m_qualy[kk+num_sum],0,n_depth*4);
        ref_offset = (int)(patch_ofset[patch_head[i]+idd]>>kshift);
        ref_rcdex = (int)(patch_ofset[patch_head[i]+idd]&003);
        id_read[kk+num_sum] = kk+num_sum;
        id_post[kk+num_sum] = ref_offset;
        id_rcdx[kk+num_sum] = ref_rcdex;
        if(ref_rcdex == 1)
        {
          st = max_read - n_Sbase - ref_offset;
          for(rc=0;rc<st;rc++)
             m_align[kk+num_sum][rc] = ' ';
          ed = seqp->length;
          for(rc=0;rc<ed;rc++)
          {
             m_align[kk+num_sum][st+rc] = seqp->data[rc];
             m_qualy[kk+num_sum][st+rc] = seqp->qual[rc];
/*           quality screen    */
             if(seqp->qual[rc]<qthresh)
                m_align[kk+num_sum][st+rc] = tolower(m_align[kk+num_sum][st+rc]);
          }
        }
        else
        {
          st = max_read - n_Sbase - ref_offset;
          for(rc=0;rc<st;rc++)
             m_align[kk+num_sum][rc] = ' ';
          ed = seqp->length-1;
          for(rc=0;rc<seqp->length;rc++)
          {
             if(seqp->data[ed-rc] == 'A')
               m_align[kk+num_sum][st+rc] = 'T';
             else if(seqp->data[ed-rc] == 'C')
               m_align[kk+num_sum][st+rc] = 'G';
             else if(seqp->data[ed-rc] == 'G')
               m_align[kk+num_sum][st+rc] = 'C';
             else if(seqp->data[ed-rc] == 'T')
               m_align[kk+num_sum][st+rc] = 'A';
             else
               m_align[kk+num_sum][st+rc] = seqp->data[rc];
/*           quality screen    */
             if(seqp->qual[ed-rc]<qthresh)
               m_align[kk+num_sum][st+rc] = tolower(m_align[kk+num_sum][st+rc]);
             m_qualy[kk+num_sum][st+rc] = seqp->qual[ed-rc];
          }
        }
	if(move_dir == 1)
	  pos = st+ref_offset+n_Sbase;
	else
	  pos = st+ref_offset-1;
      
	if((m_align[kk+num_sum][pos]== ' ')||(m_align[kk+num_sum][pos]== '\0'))
	  NC = 'N';
	else
	  NC = m_align[kk+num_sum][pos];
	fuzzy_ray[kk+num_sum] = NC;
	printf("fuzzy: %-3d %-3d %-3d %-3d %c %-3d %-25ld %-30s %s\n",kk+num_sum,num_mismatch,num_sect,i_contig,NC,move_dir,kmer_long,seqp->name,m_align[kk+num_sum]);
     }
}

/*   Subroutine to buildup multiple read alignment   */
/* =============================================  */
void Read_Pileup(fasta *seq, int *kmer_index, int *patch_index, char **m_align,int **m_qualy, int *patch_ofset,int *patch_2kmer,long *patch_array,long *patch_head,long *patch_nlink,long *patch_plink,long *id_nkmer,int *id_read, int *id_post, int *id_rcdx,int *id_post2,int *id_rddex,int *id_offset,int *reads_kmer,int *reads_last,int *reads_idex,int *reads_rcdx,int *reads_rclt,int *reads_loci,int *reads_lolt,int *reads_dupy,int num_sect,int kshift,int max_len,int max_read,int n_Sbase,int n_depth,int num_lasthit,long kmer_long, int i, int k,int i_node,int *last_hit)
/* =============================================  */
{
     int j,m,kk,ikk;
     fasta *seqp;
     char cons_base;
     int idd,idk,ref_rcdex,ref_offset,num_hits,rcdex1,rcdex2;
     int jj,long_hits,hit2reads[num_sect],hit2index[num_sect],hitnouniq[num_sect];
     int rcdex_nows = 0;
     void ArraySort_Int2(int n, int *arr, int *brr);
     void Four_Cases(fasta *seq,char **m_align, int **m_qualy, long *patch_nlink,long *patch_plink,long *patch_array,long *patch_head,int *patch_2kmer,int *id_read,int *id_post,int *id_post2,int *id_rddex,long *id_nkmer,int i,int idd,int idk,int kk,int max_read,int i_case,int n_Sbase,int *id_offset,int offset,long kmer_long);
     long next_link,last_link;

     memset(kmer_cons,'\0',max_len+1);
     memset(kmer_qual,0,(max_len+1)*4);
     memset(m_pileup,0,(max_len+1)*4);

      ikk = 0;
      long_hits = 0;
      for(m=0;m<max_len;m++)
      {
         int num_hit = 0, max_hit = 0, hit[4]= {0},hit_idx;
         hit[0] = 0;
         hit[1] = 0;
         hit[2] = 0;
         hit[3] = 0;
         hit_idx = -1;
         memset(rbase,'\0',max_len);
         for(kk=0;kk<num_sect;kk++)
         {
            if((m_align[id_read[kk]][m]!= ' ')&&(m_align[id_read[kk]][m]!= '\0'))
            {
              if((m_align[id_read[kk]][m] == 'A')||(m_align[id_read[kk]][m] == 'a'))
                hit[0]++;
              else if((m_align[id_read[kk]][m] == 'C')||(m_align[id_read[kk]][m] == 'c'))
                hit[1]++;
              else if((m_align[id_read[kk]][m] == 'G')||(m_align[id_read[kk]][m] == 'g'))
                hit[2]++;
              else if((m_align[id_read[kk]][m] == 'T')||(m_align[id_read[kk]][m] == 't'))
                hit[3]++;
              rbase[num_hit] = m_align[id_read[kk]][m];
              num_hit++;
            }
        }
        for(kk=0;kk<4;kk++)
        {
           if(hit[kk]>max_hit)
           {
             max_hit = hit[kk];
             hit_idx = kk;
           }
        }
        m_pileup[m] = num_hit;
        if(hit_idx == 0)
          kmer_cons[m] = 'A';
        else if(hit_idx == 1)
          kmer_cons[m] = 'C';
        else if(hit_idx == 2)
          kmer_cons[m] = 'G';
        else if(hit_idx == 3)
          kmer_cons[m] = 'T';
        else
          kmer_cons[m] = 'N';
        kmer_qual[m] = 40;
     }
     for(kk=0;kk<num_lasthit;kk++)
     {
        reads_kmer[kk] = reads_last[kk];
        reads_rcdx[kk] = reads_rclt[kk];
        reads_loci[kk] = reads_lolt[kk];
        reads_idex[kk] = kk;
     }
     ikk = 0;
     for(kk=0;kk<num_sect;kk++)
     {
        int low_index,uph_index;
        long kmer_next;

        idd = kmer_index[patch_head[i]+k+id_read[kk]];
        idk = patch_index[patch_head[i]+idd];
        ref_rcdex = (int)(patch_ofset[patch_head[i]+idd]&003);
        ref_offset = (int)(patch_ofset[patch_head[i]+idd]>>kshift);
        next_link = patch_nlink[patch_head[i]+idd];
        last_link = patch_plink[patch_head[i]+idd];
        reads_kmer[num_lasthit+ikk] = idk;
        reads_rcdx[num_lasthit+ikk] = ref_rcdex;
        reads_idex[num_lasthit+ikk] = num_lasthit+ikk;
        reads_loci[num_lasthit+ikk] = ref_offset;
        reads_last[ikk] = idk;
	hit2reads[kk] = idk;
	hit2index[kk] = kk;
	hitnouniq[kk] = 0;
        reads_rclt[ikk] = ref_rcdex;
        reads_lolt[ikk] = ref_offset;
        ikk++;
	if((seq+idk)->length>100)
	  long_hits++;
     }
     *last_hit = ikk;
     if(num_lasthit>0)
     {
       num_hits = num_lasthit+ikk;
       rcdex1 = 0;
       rcdex2 = 0;
       if(num_hits != 2)
       {
         ArraySort_Int2(num_sect,hit2reads,hit2index);
         for(kk=0;kk<num_sect;kk++)
	 {
            jj = kk+1;
            while((jj<num_sect)&&(hit2reads[kk]==hit2reads[jj]))
            {
              jj++;
            }
            if((jj-kk)>=2)
            {
	      for(m=(kk+1);m<jj;m++)
	         hitnouniq[hit2index[m]] = 1;
	    }
            kk = jj - 1;
	 }
         ArraySort_Int2(num_lasthit+ikk,reads_kmer,reads_idex);
         for(kk=0;kk<num_hits;kk++)
         {
            jj = kk+1;
            while((jj<num_hits)&&(reads_kmer[kk]==reads_kmer[jj]))
            {
              jj++;
            }
            if((jj-kk)>=2)
            {
              int ij1,ij2,jk1,jk2,jk,jtmp;
              ij1 = reads_kmer[kk];
              ij2 = reads_kmer[jj-1];
              jk1 = reads_idex[kk];
              jk2 = reads_idex[jj-1];
              if(jk1 > jk2)
              {
                jk = jk1 - num_lasthit;
                jtmp = jk1;
                jk1 = jk2;
                jk2 = jtmp;
                ij1 = reads_kmer[jj-1];
                ij2 = reads_kmer[kk];
              }
              else
              {
                jk = jk2 - num_lasthit;
              }
              if(reads_rcdx[reads_idex[jj-1]]==reads_rcdx[reads_idex[kk]])
              {
                if(abs(reads_loci[jk1]-reads_loci[jk2])==1)
                  reads_dupy[jk] = 0;
                else
                  reads_dupy[jk] = 1;
                rcdex1++;
              }
              else
              {
                if(abs(reads_loci[jk1]-((seq+ij2)->length-n_Sbase-reads_loci[jk2]))==1)
                  reads_dupy[jk] = 0;
                else
                  reads_dupy[jk] = 1;
                rcdex2++;
              }
            }
            else
            {
              if(reads_idex[kk] < num_lasthit)
              {
            //    reads_post[reads_kmer[kk]] = 0;
              }
            }
            kk = jj - 1;
         }
       }
       else
       {
         if(reads_rcdx[reads_idex[0]]==reads_rcdx[reads_idex[1]])
           rcdex1 = 1;
         else
           rcdex2 = 1;
       }
     }
     else
     {
       rcdex1 = 10;
       rcdex2 = 0;
     }

     ikk = 0;
     if(print_flag)
       printf("============================= %d %d\n",num_cons,num_lasthit);
     for(kk=0;kk<num_sect;kk++)
     {
        int low_index,uph_index,num_setcover = cfactor*num_cover;
        long kmer_next;

        idd = kmer_index[patch_head[i]+k+id_read[kk]];
        idk = patch_index[patch_head[i]+idd];
        ref_offset = (int)(patch_ofset[patch_head[i]+idd]>>kshift);
        ref_rcdex = (int)(patch_ofset[patch_head[i]+idd]&003);
        next_link = patch_nlink[patch_head[i]+idd];
        last_link = patch_plink[patch_head[i]+idd];
        reads_kmer[num_lasthit+ikk] = idk;
        reads_rcdx[num_lasthit+ikk] = ref_rcdex;
        reads_idex[num_lasthit+ikk] = num_lasthit+ikk;
        reads_loci[num_lasthit+ikk] = ref_offset;
        i_kmer2 = patch_kmer2[patch_head[i]+idd];
	reads_maps[reads_used] = idk;
	reads_used++;
	if(reads_mast[idk] == 0)
	  reads_mast[idk] = num_cons;
	if(hitnouniq[kk] == 0)
	{
	  if((num_cons-reads_mast[idk]) > (seq+idk)->length)
	  {
	    reads_mast[idk] = num_cons;
	    reads_post[idk] = 0;
	  }
          reads_post[idk]++;
	}
        if(print_flag)
//          printf("%d %d %-4d %-30s %-55s ",num_sect,ref_rcdex,reads_post[idk],(seq+idk)->name,m_align[id_read[kk]]);
          printf("%d %d %-4d %-30s %-55s ",num_sect,ref_rcdex,ref_offset,(seq+idk)->name,m_align[id_read[kk]]);
        if(num_pair[idk]==1)
        {
          if(num_cons>30)
          {
            int rd = mat_pair[idk];
            if((rpf_contig[rd] == i_contig)&&(rpf_offset[rd]>0)&&((num_cons-rpf_offset[rd])>50)&&((num_cons-rpf_offset[rd])<read_length)&&(num_sect<num_setcover))
            {
              rpf_contig[idk] = i_contig;
              rpf_offset[idk] = num_cons;
            }
          }
          if(num_sect<num_setcover)
          {
            rpf_contig[idk] = i_contig;
            if(rpf_offset[idk]==0)
              rpf_offset[idk] = num_cons;
          }
        }
        rcdex_nows = 0;
        if(rcdex1>=rcdex2)
        {
          if(rcdex_pres==0)
            rcdex_nows = 0;
          else
            rcdex_nows = 1;
        }
        else
        {
          if(rcdex_pres==0)
            rcdex_nows = 1;
          else
            rcdex_nows = 0;
        }
        if(rcdex_nows == 0)
        {
          if(ref_rcdex==1)
          {
            if(next_link>=0)
            {
	      int i_case = 1;
              move_dir = 1;
              Four_Cases(seq,m_align,m_qualy,patch_nlink,patch_plink,patch_array,patch_head,patch_2kmer,id_read,id_post,id_post2,id_rddex,id_nkmer,i,idd,idk,kk,max_read,i_case,n_Sbase,id_offset,ref_offset,kmer_long);
	      id_offset[hit_ikk-1] = ref_offset+1;
            }
            else
              if(print_flag) printf("\n");
          }
          else
          {
            if(last_link>=0)
            {
	      int i_case = 2;
              move_dir = 1;
              Four_Cases(seq,m_align,m_qualy,patch_nlink,patch_plink,patch_array,patch_head,patch_2kmer,id_read,id_post,id_post2,id_rddex,id_nkmer,i,idd,idk,kk,max_read,i_case,n_Sbase,id_offset,ref_offset,kmer_long);
	      id_offset[hit_ikk-1] = ref_offset+1;
            }
            else
              if(print_flag) printf("\n");
          }
        }
        else
        {
          if(ref_rcdex==2)
          {
            if(next_link>=0)
            {
	      int i_case = 3;
              move_dir = 2;
              Four_Cases(seq,m_align,m_qualy,patch_nlink,patch_plink,patch_array,patch_head,patch_2kmer,id_read,id_post,id_post2,id_rddex,id_nkmer,i,idd,idk,kk,max_read,i_case,n_Sbase,id_offset,ref_offset,kmer_long);
	      id_offset[hit_ikk-1] = (seq+idk)->length - ref_offset - n_Sbase + 1;
            }
            else
              if(print_flag) printf("\n");
          }
          else
          {
            if(last_link>=0)
            {
	      int i_case = 4;
              move_dir = 2;
              Four_Cases(seq,m_align,m_qualy,patch_nlink,patch_plink,patch_array,patch_head,patch_2kmer,id_read,id_post,id_post2,id_rddex,id_nkmer,i,idd,idk,kk,max_read,i_case,n_Sbase,id_offset,ref_offset,kmer_long);
	      id_offset[hit_ikk-1] = (seq+idk)->length - ref_offset - n_Sbase + 1;
            }
            else
              if(print_flag) printf("\n");
          }
        }
     }
     rcdex_pres = rcdex_nows;
     if(print_flag)
       printf("---------------------------------------------------------------------------------\n");
     if(move_dir==1)
     {
       long i_base = kmer_long >> (2*(n_Sbase-1));
       if(i_base==0)
         cons_base = 'A';
       else if(i_base==1)
         cons_base = 'C';
       else if(i_base==2)
         cons_base = 'G';
       else if(i_base==3)
         cons_base = 'T';
       else
         cons_base = 'C';
     }
     else
     {
       long i_base = kmer_long&003;
       if(i_base==0)
         cons_base = 'T';
       else if(i_base==1)
         cons_base = 'G';
       else if(i_base==2)
         cons_base = 'C';
       else if(i_base==3)
         cons_base = 'A';
       else
         cons_base = 'C';
     }
     if(print_flag)
       printf("%-7d %d %s   %-49s %d %d %c %c\n",num_cons,i_contig,(seq+n_list[i_contig])->name,kmer_cons,hit_ikk,i_node,cons_base,kmer_baseF[i_node]);
     cons_reads[i_contig][num_cons-1] = cons_base;
     if(kmer_baseF[i_node]=='n')
     {
       kmer_baseF[i_node] = cons_base;
       //kmer_nodes[i_prev] = i_node;
     }
//     hit_ikk = ikk;
}


/*   Subroutine to deal with four cases   */
/* =============================================  */
void Four_Cases(fasta *seq,char **m_align, int **m_qualy, long *patch_nlink, long *patch_plink, long *patch_array,long *patch_head,int *patch_2kmer,int *id_read,int *id_post,int *id_post2,int *id_rddex,long *id_nkmer,int i,int idd,int idk,int kk,int max_read,int i_case,int n_Sbase,int *id_offset,int offset,long kmer_long)
/* =============================================  */
{
     int m,ikk = hit_ikk;
     char *ss,*xx;
     int *sq,*xy,uph_index,low_index,i_kmer,num_sect2;
     long kmer_next;
     int jj,jk;

     if(i_case<=2)
     {
       ss = m_align[id_read[kk]]+max_read;
       sq = m_qualy[id_read[kk]]+max_read;
     }
     else
     {
       ss = m_align[id_read[kk]]+max_read-n_Sbase-1; 
       sq = m_qualy[id_read[kk]]+max_read-n_Sbase-1; 
     }
     xx = ss;
     xy = sq;
     if((i_case%2)==0)
     {
       uph_index = (int)(patch_plink[patch_head[i]+idd]>>nshift);
       low_index = (int)(patch_plink[patch_head[i]+idd]&hmask);
     }
     else
     {
       uph_index = (int)(patch_nlink[patch_head[i]+idd]>>nshift);
       low_index = (int)(patch_nlink[patch_head[i]+idd]&hmask);
     }
     i_kmer = patch_2kmer[patch_head[uph_index]+low_index];
     jj = kmer_idgrp[i_kmer];
     jk = kmer_locus[i_kmer];
     kmer_next = patch_array[patch_head[jj]+jk];
     num_sect2 = kmer_nhits[i_kmer];

     neib_qual[ikk] = *(sq);
     memset(neib_kmer[ikk],'\0',nkm_len);
     neib1_kmer[ikk] = *(ss);
     kmer2_hits[ikk] = num_sect2;
     for(m=0;m<nkm_len;m++)
     {
        if((*ss!='\0')&&(*sq>=qthresh))
          neib_kmer[ikk][m] = *(ss);
        else
          neib_kmer[ikk][m] = 'N';
	if(i_case<=2)
	{
          ss++;
          sq++;
	}
	else
	{
          ss--;
          sq--;
	}
     }
     if(print_flag)
     {
/*       if(i_case==1)
         printf("w1: %ld %ld %d %d %s\n",kmer_long,kmer_next,i_kmer,num_sect2,neib_kmer[ikk]);
       else if(i_case==2)
         printf("w2: %ld %ld %d %d %s\n",kmer_long,kmer_next,i_kmer,num_sect2,neib_kmer[ikk]);
       else if(i_case==3)
         printf("w3: %ld %ld %d %d %s\n",kmer_long,kmer_next,i_kmer,num_sect2,neib_kmer[ikk]);
       else if(i_case==4)
         printf("w4: %ld %ld %d %d %s\n",kmer_long,kmer_next,i_kmer,num_sect2,neib_kmer[ikk]);
           */
       if(i_case==1)
         printf("w1: %ld %c %d %d %s %d\n",kmer_long,*xx,i_kmer,num_sect2,neib_kmer[ikk],*xy);
       else if(i_case==2)
         printf("w2: %ld %c %d %d %s %d\n",kmer_long,*xx,i_kmer,num_sect2,neib_kmer[ikk],*xy);
       else if(i_case==3)
         printf("w3: %ld %c %d %d %s %d\n",kmer_long,*xx,i_kmer,num_sect2,neib_kmer[ikk],*xy);
       else if(i_case==4)
         printf("w4: %ld %c %d %d %s %d\n",kmer_long,*xx,i_kmer,num_sect2,neib_kmer[ikk],*xy);
     }
     id_read[ikk] = ikk;
     id_post[ikk] = i_kmer;
     id_post2[ikk] = i_kmer;
     id_rddex[ikk] = idk;
     id_nkmer[ikk] = kmer_next;
     ikk++;
     hit_ikk = ikk;
}

/*   Subroutine to sort the DNA sequences into a matrix table   */
/* =============================================  */
void Readpair_Stage(int nRead,char **argv, int argc, int args)
/* =============================================  */
{
     int i,j,k,nseq;
     FILE *fpMate;
     char *ptr,base[20],zero[20]={0},line2[500]={0},line[500]={0};
     int seqc;

     nseq = nRead;
     if((num_pair= (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - num_pair\n");
       exit(1);
     }
     if((mat_pair= (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - mat_pair\n");
       exit(1);
     }
     if((ins_pair= (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - ins_pair\n");
       exit(1);
     }
     if((dev_pair= (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - dev_pair\n");
       exit(1);
     }
     if((rpf_contig= (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - rpf_contig\n");
       exit(1);
     }
     if((rpf_offset= (int *)calloc(nseq,sizeof(int))) == NULL)
     {
       printf("fmate: calloc - rpf_offset\n");
       exit(1);
     }

     if((fpMate = fopen(argv[args],"r")) == NULL)
     {
       printf("Error fmate: mate input\n");
       exit(1);
     }
     seqc = 0;
     while(!feof(fpMate))
     {
       int id_pt=0,nPair;

       fgets(line,500,fpMate);
       if(feof(fpMate)) break;
       nPair=0;
       strcpy(line2,line);
       for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),nPair++)
       {
       }
       i=0;
       for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
       {
          if(i==0)
          {
            strcpy(base,zero);
            strcat(base,ptr);
            id_pt = atoi(base);
          }
          else if(i==1)
          {
            strcpy(base,zero);
            strcat(base,ptr);
            num_pair[id_pt] = atoi(base);
          }
          else if(i==2)
          {
            strcpy(base,zero);
            strcat(base,ptr);
            mat_pair[id_pt] = atoi(base);
          }
          else if((nPair>2)&&(i==(nPair-2)))
          {
            strcpy(base,zero);
            strcat(base,ptr);
            ins_pair[id_pt] = atoi(base);
          }
          else if((nPair>2)&&(i==(nPair-1)))
          {
            strcpy(base,zero);
            strcat(base,ptr);
            dev_pair[id_pt] = atoi(base);
          }
       }
       seqc++;
     }
     fclose(fpMate);
     printf("mates found: %d %d\n",nRead,seqc);
}

/*   get DNA bases from a given Kmer integer   */
/*   ===============================================  */
void Get_Kmerbase(long IntSeg,int n_Sbase,char *kmerch,char *kmerrc)
/*   ===============================================  */
{
     int i,j,k,IntBase;
     long mask,ns,nt;
     int nn[100];

     k=2*n_Sbase;
     for(k=0;k<2*n_Sbase;k++)
     {
        mask=(1L<<(k+1))-1;
        ns=IntSeg>>(k);
        nt=ns&01;
	nn[k]=ns&01;
     }
     k=0;
     while(k<=2*n_Sbase)
     {
         if(nn[k]==0&nn[k+1]==0)
	 {
	   kmerch[(2*n_Sbase-k)/2-1]='A';
	   kmerrc[k/2]='T';
	 }
	 else if (nn[k]==1&nn[k+1]==0)
	 {
	   kmerch[(2*n_Sbase-k)/2-1]='C';
	   kmerrc[k/2]='G';
	 }
	 else if (nn[k]==0&nn[k+1]==1)
	 {
	   kmerch[(2*n_Sbase-k)/2-1]='G';
	   kmerrc[k/2]='C';
	 }
	 else if (nn[k]==1&nn[k+1]==1)
         {
           kmerch[(2*n_Sbase-k)/2-1]='T';
	   kmerrc[k/2]='A';
         }
         k=k+2;
     }
}

/* creat an int matrix with subscript ange m[nrl...nrh][ncl...nch]  */
int  **imatrix(long nrl,long nrh,long ncl,long nch)
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        int  **m;

        /* allocate pointers to rows        */
        if((m=(int **)calloc(nrow,sizeof(int*)))==NULL)
        {
           printf("error imatrix: calloc error No. 1 \n");
           return(NULL);
        }
        m+=0;
        m-=nrl;

        /* allocate rows and set pointers to them        */
        if((m[nrl]=(int *)calloc(nrow*ncol,sizeof(int)))==NULL)
        {
           printf("error imatrix: calloc error No. 2 \n");
           return(NULL);
        }
        m[nrl]+=0;
        m[nrl]-=ncl;

        for(i=nrl+1;i<=nrh;i++)
           m[i]=m[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return m;
}


/* creat an int matrix with subscript ange m[nrl...nrh][ncl...nch]  */
unsigned int     **uimatrix(long nrl,long nrh,long ncl,long nch)
{
        long i, nrow=((nrh-nrl+1)/100 + 1)*100,ncol=nch-ncl+1;
        unsigned int  **m;
	long nri,nrn=nrow/100;

        /* allocate pointers to rows        */
        if((m=(unsigned int **)calloc(nrow,sizeof(int*)))==NULL)
        {
           printf("error imatrix: calloc error No. 1 \n");
           return(NULL);
        }

        /* allocate rows and set pointers to them        */
	/* allocate in 100 batches to use freed memory */
	nrl = 0;
	for(nri=0;nri<100;nri++,nrl+=nrn) {
           if((m[nrl]=(unsigned int *)calloc(nrn*ncol,sizeof(int)))==NULL)
           {
              printf("error imatrix: calloc error No. 2 \n");
              return(NULL);
           }

           for(i=1;i<nrn;i++)
              m[i+nrl]=m[i+nrl-1]+ncol;
	}
       /* return pointer to array of pointers to rows   */
        return m;
}

/* creat an int matrix with subscript ange m[nrl...nrh][ncl...nch]  */
long     **limatrix(long nrl,long nrh,long ncl,long nch)
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        long  **m;

        /* allocate pointers to rows        */
        if((m=(long **)calloc(nrow,sizeof(long*)))==NULL)
        {
           printf("error imatrix: calloc error No. 1 \n");
           return(NULL);
        }
        m+=0;
        m-=nrl;

        /* allocate rows and set pointers to them        */
        if((m[nrl]=(long *)calloc(nrow*ncol,sizeof(long)))==NULL)
        {
           printf("error imatrix: calloc error No. 2 \n");
           return(NULL);
        }
        m[nrl]+=0;
        m[nrl]-=ncl;

        for(i=nrl+1;i<=nrh;i++)
           m[i]=m[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return m;
}


/* creat char matrix with subscript ange cm[nrl...nrh][ncl...nch]  */
char    **cmatrix(long nrl,long nrh,long ncl,long nch)
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        char **cm;

        /* allocate pointers to rows        */
        if((cm=(char **)calloc(nrow,sizeof(char*)))==NULL)
        {
           printf("error cmatrix: calloc error No. 1 \n");
           return(NULL);
        }
        cm+=0;
        cm-=nrl;

        /* allocate rows and set pointers to them        */
        if((cm[nrl]=(char *)calloc(nrow*ncol,sizeof(char)))==NULL)
        {
           printf("error cmatrix: calloc error No. 2 \n");
           return(NULL);
        }
        cm[nrl]+=0;
        cm[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           cm[i]=cm[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return cm;
}

#define SWAP(a,b) temp=(a);(a)=b;(b)=temp;

/*   Subroutine to sort an array arr[0,...,n-1] into ascending order while
     making the corresponding reaarangement of the array brr[0,...,n-1]
     by the use of Quicksort (Sedgwick, R. 1978, Communications o fthe ACM,
     vol. 21, pp. 847-857) also see Numerical Recipes in C   

== =============================== */
void ArraySort_Mix(int n, long *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/* =============================== */
void ArraySort_Int2(int n, int *arr, int *brr)
/* =============================== */
{         
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/*   function to sort an array into a decreasing order:  a>b>c>....    */
/* =============================== */
void ArraySort2_Int2(int n, int *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]>=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]<arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]<arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]<arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]>a);
             do j--; while (arr[j]<a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}
 
