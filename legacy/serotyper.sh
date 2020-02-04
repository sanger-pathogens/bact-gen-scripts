#! /usr/local/bin/bash

# Dealing with command line input

if [ $# -ge 1 ]; then
		echo "Reference sequence is $1";
	else
		echo "Need an input file containing the reference sequence in multifasta format and list of forward fastqs";
		exit
fi

#Indexing reference sequence

cat $1 | seqret -filter > unindexed_sequence.mfa;
grep -v ">" $1 | seqret -filter > unindexed_sequence.dna;
bwa index unindexed_sequence.dna 2> verbiage;

# Filling job array

declare -i N;
N=0;

for f in $@; do
	FILE=JobArray.$N;
	echo "${f}" > $FILE;
	N+=1;
done

N=$N-1;

# Running job array

echo bash '~nc3/Scripts/bwa_wrapper.sh unindexed_sequence.dna < JobArray.${LSB_JOBINDEX}' | bsub -M 2000000 -R 'select[mem>2000] rusage[mem=2000]' -J "BWAmap[1-$N]%25" -o verbiage -e verbiage;

echo "Strain		Cps locus	Coverage (%)" > serotypes.summary;

#bsub -w "ended(BWAmap[$N])" /nfs/pathogen/sh16_scripts/serotyping_tidyup.sh;
