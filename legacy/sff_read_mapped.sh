#!/bin/bash

usage()
{
cat << EOF

USAGE: $0 [options] <sff files>

Finds scaffolding errors (if ref is concatenated scaffolds) or genomic rearrangements (if ref is a similar strain).

OPTIONS:
   -d [int]		Standard deviation of insert sizes
   -i [int]		Mean insert length 			[Default = 1000]
   -r [file name]	Reference fastq file
   -h			Show this message

   
EOF
}

INSERT=3000
SDEV=1000
SFF=
REF=
while getopts “hi:d:s:r:” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         i)
             INSERT=$OPTARG
             ;;
         d)
             SDEV=$OPTARG
             ;;
         r)
             REF=$OPTARG
             ;;
         ?)
             usage
             exit
             ;;
     esac
done
shift $(($OPTIND - 1))
SFF=$*


FASTQ=${SFF%.sff}.fastq

#extract fasta and fasta.qual from sff using linker
~sh16/mira_2.9.37_dev_linux-gnu_x86_64/3rdparty/sff_extract -c -o $FASTQ -l ~tdo/work/linker.fasta -i "insert_size:$INSERT,insert_stdev:$SDEV" $SFF

PERL5LIB=$PERL5LIB:~tdo/bin/oldPerl

export PERL5LIB
#turn fasta+fasta.qual to fastq
perl -w -e "use AssemblyTools_unstable;AssemblyTools_unstable::fasta2fastq(  '${FASTQ}.fasta', '${FASTQ}' );"

#rename forward and reverse reads
cat $FASTQ | perl -nle 's/part1/F/g;s/part2/R/g; print' > ${FASTQ}.mate

#reverse complement forward reads
~tdo/Bin/revcompForwardRead.pl ${FASTQ}.mate


/nfs/pathdata/Streptococcus/pneumoniae/Transcriptome/bin/pileup_v0.4/pileup.csh -rtype 454 -paired 0 -cigar 1 ${FASTQ}.mate.revcompF $REF ${FASTQ}_revcomp

#to get rid of all same mapping read pairs
perl ~tdo/Bin/nico.helper.trivial.pl ${FASTQ}_revcomp.cigar > ${FASTQ}.cigar3

#to only select reads mapping with maximum quality
awk '{if ($1 == "cigar::50") print}' ${FASTQ}.cigar3 > ${FASTQ}50_output

#/nfs/team81/tdo/Bin/carma.Find.error.start.sh  <stuff.cigar3 file> <readLength> optional <maxInsert>

~tdo/Bin/carma.Find.error.start.sh  ${FASTQ}50_output 20 10000 