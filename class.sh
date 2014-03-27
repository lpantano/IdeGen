#!/bin/sh
set -e 
if [ -e $1.bam ]
then
        FILE="$1.bam" 
else
        FILE="$1.sam"
fi
samtools depth $FILE > $1.depth
junc $FILE > $1.splice 
if [ $# -eq 2 ]
then
        class $1 --set-cover -F 0.01  > $2
else
        class $1 --set-cover -F 0.01   
fi
rm $1.depth
rm $1.splice
