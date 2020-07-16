#!/bin/bash

PATH_TO_OUTPUT=$HOME/RNA-seq/Licencjat/Kallisto/Output/
PATH_TO_INDEX=$HOME/RNA-seq/Licencjat/Kallisto/Index/Homo_sapiens.GRCh38.cdna.ncrna.index

for f in *; 
do
	kallisto quant -i "$PATH_TO_INDEX" -o "$PATH_TO_OUTPUT$f"0 -b 100 $f $f;
done
