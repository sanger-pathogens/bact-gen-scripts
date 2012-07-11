#!/bin/bash

usage()
{
cat << EOF

USAGE: $0 [options] 

Runs abyss.

OPTIONS:
   -c [int] C option for Abyss [Default = 2]
   -d [directory]	Specify a directory for the abyss output [Default = new directory in current location using fastq prefix]
   -f [fastq file]	Fastq file containing (forward) reads
   -i [int]		Mean insert length 			[Default = 180]
   -k [int]		Starting kmer value 			[Default = 2/3 of read length]
   -n [int]		N option for Abyss [Default = 10]
   -r [fastq file]	Fastq file containing reverse reads	[Optional]
   -o ["string"]	Any other abyss commands
   -h			Show this message

NOTE:
   Please bsub this command. I if you expect it to use large amount of memory
   please include the relevant bsub options (see example below). You can change
   the number of processors Abyss will run on using the -n flag: 
   bsub -R 'select[mem>24000] rusage[mem=24000]' -e D7.err -o D7.out -M 30000000 -n 4 

   
EOF
}


SCRIPT_DIR="/nfs/users/nfs_s/sh16/scripts/"
ABYSS_DIR="/software/pathogen/external/applications/ABySS/ABySS/bin/"

KMER=0
FORWARD=
REVERSE=
OTHERS=
ABYSSDIR=
C=2
N=10
while getopts “hic::e:f:k:n:pr:o:d” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         i)
             INSERT=$OPTARG
             ;;
         c)
             C=$OPTARG
             ;;
         c)
             N=$OPTARG
             ;;
         d)
             ABYSSDIR=$OPTARG
             ;;
         e)
             EXPCOV=$OPTARG
             ;;
         f)
             FORWARD=$OPTARG
             ;;
         k)
             KMER=$OPTARG
             ;;
         r)
             REVERSE=$OPTARG
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


if [ -z "$FORWARD" ]
then
	echo
	echo ERROR: No fastq file selected
	echo
    usage
    exit 1
fi

readlength=$(head -n 2 $FORWARD | tail -n 1 | wc -m)
if [ "$KMER" -eq "0" ]
	then KMER=$[$readlength*2/3]
fi 

#echo $KMER

filename=${FORWARD##*/}
filename=${filename%_[12].fastq}
filename=${filename%_[12]_nonhuman.fastq}

if [ -z "$ABYSSDIR" ]
	then ABYSSDIR=${filename}_abyss/
fi

if [ -d $ABYSSDIR ]; then echo $ABYSSDIR found.
#do nothing
else
echo $ABYSSDIR not found. Creating...
mkdir $ABYSSDIR
fi



if [ -e $ABYSSDIR${filename}.fasta ]; then echo $ABYSSDIR${filename}.fasta found.
else
echo Creating fasta file...

${SCRIPT_DIR}fastq2fasta.pl $FORWARD $ABYSSDIR${filename}_f.fasta

if [ -e "$REVERSE" ]
	then ${SCRIPT_DIR}fastq2fasta.pl $REVERSE $ABYSSDIR${filename}_r.fasta
	cat $ABYSSDIR${filename}_f.fasta $ABYSSDIR${filename}_r.fasta > $ABYSSDIR${filename}.fasta
	rm $ABYSSDIR${filename}_r.fasta
fi

rm $ABYSSDIR${filename}_f.fasta
fi

cd $ABYSSDIR

echo Running Abyss...

#${ABYSS_DIR}abyss-pe k=$KMER l=$readlength c=$C n=$N in=${filename}.fasta lib=${filename}
/software/pathogen/projects/protocols/bin/my_abyss_pe_svn.pl -k $KMER -n $N -i ${filename}.fasta -o ${filename}
exit
awk '$2 >= 500 {print $0}' ${filename}-contigs.fa > ${filename}-contigs_gt500.headers

fgrep -f ${filename}-contigs_gt500.headers -A 1 ${filename}-contigs.fa > ${filename}_LargeContigs.fna

rm ${filename}-contigs_gt500.headers ${filename}.fasta

cd ..