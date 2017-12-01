#!/bin/bash

bam=false

## sam or bam?
while true; do
    read -p "do you wish to use BAM or SAM files? (bam/sam): " bs
    case $bs in
	"bam" ) echo "BAM files chosen."; bam=true; break;;
	"sam" ) echo "SAM files chosen."; bam=false; break;;
	* ) echo "Please answer \"bam\" or \"sam\".";;
    esac
done

## if bam, convert to sam
if $bam; then
    module load  GCC/5.4.0-2.26  OpenMPI/1.10.3 SAMtools/1.4
    for bamfile in *bam;
    do 	samfile=$(basename $bamfile .bam).sam
	echo "Converting BAM: $bamfile"
	echo "To SAM: $samfile"
	samtools view -h $bamfile > $samfile
    done    
fi


## calculate statistics for each sam file
echo ">>> Iterate through all *sam files in directory, calculating average and sd of insert size"
touch insertSize.stats.txt
for file in *sam
do echo "> Calculating average insert size for $file.."
    mean=$(cat $file | awk '{if ($9 > 0 ) {S+=$9;T+=1}}END{print "Mean: " S/T}')
    echo "$file $mean" >> insertSize.stats.txt
done
    
    
