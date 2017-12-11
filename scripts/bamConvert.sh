#!/bin/sh
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -A lsens2017-3-2
#SBATCH -p dell
#SBATCH -t 05:00:00
#SBATCH -C mem256GB
#SBATCH -J bamToBw.%j
#SBATCH -o bamToBw.%j.out
#SBATCH -e bamToBw.%j.err

ml foss/2016b
ml Python/3.5.2
ml GCC/5.4.0-2.26  OpenMPI/1.10.3
ml SAMtools/1.4

bam=$1
echo ">> The current bam will be processed:"
echo $bam

## Ask user if sorted needed
sorted=true
while true; do
    read -p ">> Is bamfile sorted by coord? (y/n): " yn
    case $yn in
	"y" ) echo "> Bam file is sorted."; sorted=true; break;;
	"n" ) echo "> Bam is not sorted, so sorting needed.."; sorted=false; break;;
	* ) echo "Please answer \"y\" or \"n\".";;
    esac
done

# get sample name
sample=$(basename $bam .bam)
echo "Sample: $sample"


if [ sorted==true ]
then 
    bamsort=$bam
else
    bamsort=${sample}.sorted.bam
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
    samtools index -b $bamsort $bamsort.bai
fi 

# make bigwig
if [ -f $sample.bw ]
then
    echo "> Bigwig file for $sample exist.. Nothing done"
else
    echo ">>Generating bigwig file from $bam.."
    bamCoverage -b $bamsort -o $sample.bw
fi
echo ">>Indexing and bigwig-generation completed for $bam..!"
echo " ++++++++++++++++++++++++++++++++++++++++++++++++++++ "

# make bed

