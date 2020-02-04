#!/bin/tcsh
# by zn1@sanger.ac.uk
# 11/07/2008
#

set locPath = `pwd`

cd $locPath/fuzzy; make
cd $locPath/fastq2fasta; make

cd $locPath
awk -v path=$locPath '/^# set fuzzyProgramPath/{printf "set prog=%s/ ",path}//{print $0}' fuzzypath.csh_src >! fuzzypath.csh
chmod 755 fuzzypath.csh

