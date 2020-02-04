#!/bin/bash

echo file %C %G %A %T %GC

for FILE in "$@"
do

C=$(CHAR="C"; cat $FILE | grep -v "[^ACGTN]" | tr '\n' 'R' | tr $CHAR "\n" | wc -l)
G=$(CHAR="G"; cat $FILE | grep -v "[^ACGTN]" | tr '\n' 'R' | tr $CHAR "\n" | wc -l)
A=$(CHAR="A"; cat $FILE | grep -v "[^ACGTN]" | tr '\n' 'R' | tr $CHAR "\n" | wc -l)
T=$(CHAR="T"; cat $FILE | grep -v "[^ACGTN]" | tr '\n' 'R' | tr $CHAR "\n" | wc -l)

TOTAL=$[$A+$G+$C+$T]

Cpercent=$(echo "scale=4; 100*$C/$TOTAL" | bc)
Gpercent=$(echo "scale=4; 100*$G/$TOTAL" | bc)
Apercent=$(echo "scale=4; 100*$A/$TOTAL" | bc)
Tpercent=$(echo "scale=4; 100*$T/$TOTAL" | bc)

GC=$(echo "scale=4; $Cpercent+$Gpercent" | bc)


echo $FILE $Cpercent $Gpercent $Apercent $Tpercent $GC

done