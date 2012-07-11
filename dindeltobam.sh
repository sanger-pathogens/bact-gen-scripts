#!/bin/bash 
samtools view $1 >> $2.sam 
echo $3 $4 $5 >> $2.windows 
rm $1
