#!/usr/sh

file=$1

echo $file; 
basename=$(basename $file _R1.bw); 
echo "track type=bigWig name=${basename} visibility=full smoothingWindow=4 autoScale=on description=${basename} bigDataUrl=http://bone.bmc.lu.se/Public/$file" >> UCSC.bw.txt ; 
echo $basename; 

