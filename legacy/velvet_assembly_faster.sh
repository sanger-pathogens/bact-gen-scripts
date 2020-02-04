
if [ "$#" -lt "2" ]
then echo "Usage:"
echo "for single end reads:"
echo "velvet_assembly.sh <fastq file> <minimum expected coverage>"
echo "for paired end reads:"
echo "velvet_assembly.sh <name to give shuffled fastq file> <forward read fastq file> <reverse read fastq file> <mean insert length> <minimum expected coverage>"
echo "or if you already have the shuffled sequences:"
echo "velvet_assembly.sh <name of shuffled fastq file> <mean insert length> <minimum expected coverage>"
exit
fi

lastcov=20
eval covlimit=\$$#
mincov=$[$covlimit/2] #change this to change the minimum coverage that can be accepted
covcutoff=mincov
kmer=27
insert=180
readtype="short"

logfile=${1%.fastq}.log

echo velvet_assembly_faster.sh $@ > $logfile

if [ "$#" -gt "2" ]
	then readtype="shortPaired"
	insert=$2
fi
if [ "$#" -gt "3" ]
	then echo Shuffling sequences...
	/nfs/pathogen/sh16_scripts/shufflefastqSequences.pl $2 $3 $1
fi

echo Running velveth with kmer of $kmer

velveth ${1%.fastq}_velvet $kmer -fastq -$readtype $1 > ${1%.fastq}_velvet.log

echo Running velvetg

if [ "$#" -gt "2" ]
	then velvetg ${1%.fastq}_velvet -ins_length $insert > ${1%.fastq}_velvet.log
fi
if [ "$#" -lt "3" ]
	then velvetg ${1%.fastq}_velvet > ${1%.fastq}_velvet.log
fi
expcov=$(/nfs/pathogen/sh16_scripts/velvet_stats_2_av_cov.py ${1%.fastq}_velvet/stats.txt $mincov)

echo Found expected coverage of $expcov
echo kmer=$kmer expected coverage=$expcov  >> $logfile
grep "n50" ${1%.fastq}_velvet.log | head -n 1
grep "n50" ${1%.fastq}_velvet.log | head -n 1 >> $logfile

if [ "$expcov" -gt $[$covlimit-1] ]

	then kmer=$[$kmer+2]
	while [ "$kmer" -lt "32" ]
	do
		echo Expected coverage above $covlimit, increasing kmer
		echo Running velveth with kmer of $kmer
		velveth ${1%.fastq}_velvet $kmer -fastq -$readtype $1 > ${1%.fastq}_velvet.log
		echo Running velvetg
		if [  "$#" -gt "2" ]
			then velvetg ${1%.fastq}_velvet -ins_length $insert  > ${1%.fastq}_velvet.log
		fi
		if [ "$#" -lt "3" ]
			then velvetg ${1%.fastq}_velvet  > ${1%.fastq}_velvet.log
		fi
		expcov=$(/nfs/pathogen/sh16_scripts/velvet_stats_2_av_cov.py ${1%.fastq}_velvet/stats.txt $mincov)
		echo Found expected coverage of $expcov
		echo kmer=$kmer expected coverage=$expcov  >> $logfile
		grep "n50" ${1%.fastq}_velvet.log | head -n 1
		grep "n50" ${1%.fastq}_velvet.log | head -n 1 >> $logfile
		if  [ "$expcov" -lt $covlimit ]
			then echo Found best combination:
			if  [ "$kmer" -gt "22" ]
				then kmer=$[$kmer-2]
				expcov=$lastcov
			fi
			echo kmer=$kmer expected coverage=$expcov
			break
		fi
		lastcov=$expcov
		kmer=$[$kmer+2]
	done

else kmer=$[$kmer-2]
	while [ "$kmer" -gt "20" ]
	do
		echo Expected coverage below $covlimit, decreasing kmer
		velveth ${1%.fastq}_velvet $kmer -fastq -$readtype $1 > ${1%.fastq}_velvet.log
		echo Running velveth with kmer of $kmer
		if [ "$#" -gt "2" ]
			then velvetg ${1%.fastq}_velvet -ins_length $insert > ${1%.fastq}_velvet.log
		fi
		if [ "$#" -lt "3" ]
			then velvetg ${1%.fastq}_velvet > ${1%.fastq}_velvet.log
		fi
		expcov=$(/nfs/pathogen/sh16_scripts/velvet_stats_2_av_cov.py ${1%.fastq}_velvet/stats.txt $mincov)
		echo Found expected coverage of $expcov
		echo kmer=$kmer expected coverage=$expcov  >> $logfile
		grep "n50" ${1%.fastq}_velvet.log | head -n 1
		grep "n50" ${1%.fastq}_velvet.log | head -n 1 >> $logfile
		if  [ "$expcov" -gt $[$covlimit-1] ]
			then echo Found best combination:
			echo kmer=$kmer expected coverage=$expcov
			break
		fi
		kmer=$[$kmer-2]
	done
fi

if [ "$kmer" -lt "20" ]
	then kmer=21
fi
if [ "$kmer" -gt "31" ]
        then kmer=31
fi

covcutoff=$[$expcov/2]

echo Final run parameters: >> $logfile
echo kmer=$kmer expected coverage=$expcov coverage cutoff=$covcutoff  >> $logfile 

echo Running final run with kmer of $kmer, expected coverage of $expcov and coverage cutoff of $covcutoff

if [ "$#" -gt "2" ]
	then velvetg ${1%.fastq}_velvet -ins_length $insert -min_contig_lgth 100 -exp_cov $expcov -cov_cutoff $covcutoff > ${1%.fastq}_velvet.log
fi
if [ "$#" -lt "3" ]
	then velvetg ${1%.fastq}_velvet -min_contig_lgth 100 -exp_cov $expcov -cov_cutoff $covcutoff > ${1%.fastq}_velvet.log
fi

grep "n50" ${1%.fastq}_velvet.log | head -n 1
grep "n50" ${1%.fastq}_velvet.log | head -n 1 >> $logfile

rm ${1%.fastq}_velvet.log

echo "Assembly finished"
