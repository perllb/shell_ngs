#!/bin/sh
ml foss/2016b
ml Python/3.5.2
ml GCC/5.4.0-2.26  OpenMPI/1.10.3
ml SAMtools/1.4

bamfiles=$1
echo "The current bamfiles will be processed:"
echo $bamfiles

for bam in $bamfiles

do echo "Current bamfile: "
echo $bam
echo "..."


# get sample name
sample=$(basename $bam .bam)
echo "Sample: $sample"
bamsort=$sample.sortedCoord.bam

if [ -f $bamsort ]
then 
    echo "-- Sorted bamfile exists.."
else
    echo ">>Sorting bam file.."
    samtools sort -o $bamsort $bam
    echo "-- Bamfile sorted.."
fi



## index bam
if [ -f $bamsort.bai ] 
then
    echo "Index file for $bam exists"
else 
    echo ">>Creating index file for $bam.."
    samtools index -b $bamsort -o $bamsort.bai
fi 


# make bigwig
echo ">>Generating bigwig file from $bam.."
bamCoverage -b $bamsort -o $sample.bw

echo ">>Indexing and bigwig-generation completed for $bam..!"
echo " ++++++++++++++++++++++++++++++++++++++++++++++++++++ "

done ## processing done!



