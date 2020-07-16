#!/bin/bash

VAR1="ERR35848"
VAR2="_1.fastq.gz"
VAR3="_2.fastq.gz"
VAR4="_1.trimmed.fastq.gz"
VAR5="_2.trimmed.fastq.gz"

for i in 5 6 7 8; do

	java -jar /usr/share/java/trimmomatic.jar PE -phred33 "$VAR1$i$VAR2" "$VAR1$i$VAR3" "$VAR1$i$VAR4" /dev/null "$VAR1$i$VAR5" /dev/null ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:60

done
