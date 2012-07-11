#!/bin/bash

usage()
{
cat << EOF

USAGE: $0 [options] 

Optimises parameters in Velvet.

OPTIONS:
   -e [int]		Expected coverage			[Default = 20]
   -f [fastq file]	Fastq file containing (forward/shuffled) reads
   -i [int]		Mean insert length 			[Default = 180]
   -k [int]		Starting kmer value 			[Odd number from 21 to 31. Default = 27]
   -n			Turn off scaffolding (no Ns in assembly)
   -p			Paired-end
   -r [fastq file]	Fastq file containing reverse reads	[Optional]
   -s [file name]	Name for shuffled fastq output file 	[Optional]
   -o ["string"]	Any other velvetg commands
   -h			Show this message

EXAMPLES:
   Single end reads starting from a kmer of 23:
   $0 -f bob.fastq -k 23
   Paired end reads in two fastq files with scaffolding turned off:
   $0 -p -n -i 180 -f bob_1.fastq -r bob_2.fastq -s bob_shuffled.fastq 
   Paired end reads in a single shuffled fastq file with minimum contig length set to 100 and read tracking on:
   $0 -p -i 180 -f bob_shuffled.fastq -o "-min_contig_lgth 100 -read_trkg yes"

   
EOF
}

FASTQDIR=
MININSERT=100
MAXINSERT=300
KMER=27
NAME=
RUNNAME=
PAIRED="y"
QUALITY=30
REF=
READLENGTH=54
SKIP=2


while getopts “f:hi:j:e:k:n:pq:r:s:o:” OPTION
do
     case $OPTION in
         e)
             EXPCOV=$OPTARG
             ;;
         f)
             FASTQDIR=$OPTARG
             ;;
         h)
             usage
             exit 1
             ;;
         i)
             MININSERT=$OPTARG
             ;; 
         j)
             MAXINSERT=$OPTARG
             ;;
         k)
             KMER=$OPTARG
             ;;
         l)
             READLENGTH=$OPTARG
             ;;
         n)
             NAME=$OPTARG
             ;;
         m)
             RUNNAME=$OPTARG
             ;;
         p)
             PAIRED=1
             ;;
         q)
             QUALITY=$OPTARG
             ;;
         r)
             REF=$OPTARG
             ;;
         t)
             RTYPE=$OPTARG
             ;;
         s)
             SKIP=$OPTARG
             ;;
         o)
             OTHERS=$OPTARG
             ;;
         ?)
             usage
             exit
             ;;
     esac
done



meaninsert=(($MAXINSERT+$MININSERT)/2)+$MININSERT


if [ PAIRED eq "n"]
then
	if [ RTYPE eq "454"]
		then /nfs/users/nfs_s/sh16/ssaha_2.2/ssaha2 -$RTYPE -seeds 5 -score $QUALITY -kmer 13 -skip 4 -diff 0 -output cigar $REF ${FASTQDIR}${NAME}.fastq > ${RUNNAME}/cigar1.tmp
		else
	/nfs/users/nfs_s/sh16/ssaha_2.2/ssaha2 -$RTYPE -score $QUALITY -skip 2 -diff 0 -output cigar $REF ${FASTQDIR}${NAME}.fastq > ${RUNNAME}/cigar1.tmp
	fi

	grep "^cigar" ${RUNNAME}/cigar1.tmp > ${RUNNAME}/cigar2.tmp
	mv ${RUNNAME}/cigar2.tmp ${RUNNAME}/cigar1.tmp
			
	perl /nfs/users/nfs_s/sh16/scripts/cigar00_2plot.pl ${RUNNAME}/cigar1.tmp ${RUNNAME}/allcoverage.plot $READLENGTH
	/nfs/users/nfs_t/tdo/bin/pileup_v0.4/ssaha_pileup/ssaha_pileup/ssaha_cigar ${RUNNAME}/cigar1.tmp ${RUNNAME}/cigar2.tmp
	awk '{print $2}' ${RUNNAME}/cigar2.tmp > ${RUNNAME}/readnames.tmp
	/nfs/users/nfs_t/tdo/bin/pileup_v0.4/ssaha_pileup/other_codes/get_seqreads/get_seqreads ${RUNNAME}/readnames.tmp ${FASTQDIR}${NAME}.fastq $	{RUNNAME}/fastq.tmp
		
	/nfs/users/nfs_s/sh16/scripts/exclude/get_excreads ${RUNNAME}/readnames.tmp ${FASTQDIR}${NAME}.fastq ${NAME}unmap.fastq"
	mv ${NAME}unmap.fastq ${RUNNAME}/unmap.fastq"
	if [ RTYPE eq "454" ]
		then
		/nfs/users/nfs_t/tdo/bin/pileup_v0.4/ssaha_pileup/ssaha_pileup/ssaha_pileup -solexa 0 -trans 0 ${RUNNAME}/cigar2.tmp $REF ${RUNNAME}/fastq.tmp > ${RUNNAME}/snp.tmp
		/nfs/users/nfs_t/tdo/bin/pileup_v0.4/ssaha_pileup/ssaha_pileup/ssaha_pileup -cons 1 -solexa 0 -trans 0 ${RUNNAME}/cigar2.tmp $REF ${RUNNAME}/fastq.tmp > ${RUNNAME}/pileup.tmp
	else
		/nfs/users/nfs_t/tdo/bin/pileup_v0.4/ssaha_pileup/ssaha_pileup/ssaha_pileup -solexa 1 -trans 0 ${RUNNAME}/cigar2.tmp $REF ${RUNNAME}/fastq.tmp > ${RUNNAME}/snp.tmp
		/nfs/users/nfs_t/tdo/bin/pileup_v0.4/ssaha_pileup/ssaha_pileup/ssaha_pileup -cons 1 -solexa 1 -trans 0 ${RUNNAME}/cigar2.tmp $REF ${RUNNAME}/fastq.tmp > ${RUNNAME}/pileup.tmp
	
	egrep ^SNP ${RUNNAME}/snp.tmp > ${RUNNAME}/all.snp"
	egrep ^cons ${RUNNAME}/pileup.tmp > ${RUNNAME}/all.pileup"
			
			
		#paired end
else
	if [ RTYPE eq '454' ]
		/nfs/users/nfs_s/sh16/ssaha_2.2/ssaha2 -$RTYPE -pair ${MININSERT},${MAXINSERT} -seeds 5 -score $QUALITY -kmer 13 -skip 4 -diff 0 -output cigar $REF ${FASTQDIR}${NAME}_1.fastq ${FASTQDIR}${NAME}_2.fastq > ${RUNNAME}/cigar1.tmp
			else
		/nfs/users/nfs_t/tdo//Bin/join.solexa2pairs.v2.pl ${FASTQDIR}${NAME}_1.fastq ${FASTQDIR}${NAME}_2.fastq ${RUNNAME}/shuffled.tmp
		/nfs/users/nfs_s/sh16/ssaha_2.2/ssaha2 -$RTYPE -score $QUALITY -skip 2 -diff 0 -output cigar $REF ${RUNNAME}/shuffled.tmp > ${RUNNAME}/cigar1.tmp

	grep "^cigar" ${RUNNAME}/cigar1.tmp > ${RUNNAME}/cigar2.tmp
	mv ${RUNNAME}/cigar2.tmp ${RUNNAME}/cigar1.tmp
	perl /nfs/users/nfs_s/sh16/scripts/cigar00_2plot.pl ${RUNNAME}/cigar1.tmp ${RUNNAME}/allcoverage.plot $READLENGTH#change 36 for different read lengths
				
	# -uniq 0 maps exactly identical repeat reads randomly
	/nfs/users/nfs_t/tdo/bin/pileup_v0.4/ssaha_pileup/ssaha_pileup/ssaha_cigar -uniq 0 ${RUNNAME}/cigar1.tmp ${RUNNAME}/cigar2.tmp
	perl /nfs/users/nfs_s/sh16/scripts/cigar00_2plot.pl ${RUNNAME}/cigar2.tmp ${RUNNAME}/random_exact_repeats_coverage.plot $READLENGTH
				
	/nfs/users/nfs_t/tdo/bin/pileup_v0.4/ssaha_pileup/ssaha_pileup/ssaha_pairs -insert $meaninsert ${RUNNAME}/cigar1.tmp ${RUNNAME}/cigar2unclean.tmp
	/nfs/users/nfs_t/tdo/bin/pileup_v0.4/ssaha_pileup/ssaha_pileup/ssaha_clean -insert $meaninsert ${RUNNAME}/cigar2unclean.tmp ${RUNNAME}/cigar2.tmp
	awk '{print $2}' ${RUNNAME}/cigar2.tmp > ${RUNNAME}/readnames.tmp
				
	/nfs/users/nfs_t/tdo/bin/pileup_v0.4/ssaha_pileup/other_codes/get_seqreads/get_seqreads ${RUNNAME}/readnames.tmp ${RUNNAME}/shuffled.tmp ${RUNNAME}/fastq.tmp
				
				
	/nfs/users/nfs_s/sh16/scripts/exclude/get_excreads ${RUNNAME}/readnames.tmp ${RUNNAME}/shuffled.tmp ${NAME}unmap.fastq"
	mv ${NAME}unmap.fastq ${RUNNAME}/unmap.fastq"
			
			
	if [ RTYPE eq "454" ]
		"/nfs/users/nfs_t/tdo/bin/pileup_v0.4/ssaha_pileup/ssaha_pileup/ssaha_pileup -solexa 0 -trans 0 ${RUNNAME}/cigar2.tmp $REF ${RUNNAME}/fastq.tmp > ${RUNNAME}/snp.tmp
		"/nfs/users/nfs_t/tdo/bin/pileup_v0.4/ssaha_pileup/ssaha_pileup/ssaha_pileup -cons 1 -solexa 0 -trans 0 ${RUNNAME}/cigar2.tmp $REF ${RUNNAME}/fastq.tmp > ${RUNNAME}/pileup.tmp
	else
		"/nfs/users/nfs_t/tdo/bin/pileup_v0.4/ssaha_pileup/ssaha_pileup/ssaha_pileup -solexa 1 -trans 0 ${RUNNAME}/cigar2.tmp $REF ${RUNNAME}/fastq.tmp > ${RUNNAME}/snp.tmp
		"/nfs/users/nfs_t/tdo/bin/pileup_v0.4/ssaha_pileup/ssaha_pileup/ssaha_pileup -cons 1 -solexa 1 -trans 0 ${RUNNAME}/cigar2.tmp $REF ${RUNNAME}/fastq.tmp > ${RUNNAME}/pileup.tmp
	
fi

egrep ^SNP ${RUNNAME}/snp.tmp > ${RUNNAME}/all.snp"
egrep ^cons ${RUNNAME}/pileup.tmp > ${RUNNAME}/all.pileup"
			
			
rm ${RUNNAME}/*.tmp
gzip ${RUNNAME}/*.plot



