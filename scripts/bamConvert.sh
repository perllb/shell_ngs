#!/bin/sh
ml foss/2016b
ml Python/3.5.2
ml GCC/5.4.0-2.26  OpenMPI/1.10.3
ml SAMtools/1.4

bam=$1
sort=$2

### If sort argument is empty, prompt user 
if [ -z "$2" ]
  then
  
    echo "> No 'sort' argument input. If sort needed or not can be specified with 'sort' or 'nosort' arguments"
  ## Ask user if sorted needed
  sorted=true
  while true; do
      read -p ">> Is bamfile sorted by coord? (y/n): " yn
      case $yn in
	  "y" ) echo "> Bam file is sorted."; sorted=true; break;;
	  "n" ) echo "> Bam is not sorted, so sorting needed.."; sorted=false; break;;
	  * ) echo "> Please answer \"y\" or \"n\".";;
      esac
  done

## If sort argument is nonsort, set $sorted to true
elif [ "$sort" == 'nosort' ]
   then
   sorted=true

## If sort argument is not empty and not sort or nonsort, then prompt user for input
elif [ "$sort" != 'nosort' ] && [ "$sort" != 'sort' ]
then
  echo "> Sort argument unclear.."
  ## Ask user if sorted needed
  sorted=true
  while true; do
      read -p ">> Is bamfile sorted by coord? (y/n): " yn
      case $yn in
	  "y" ) echo "> Bam file is sorted."; sorted=true; break;;
	  "n" ) echo "> Bam is not sorted, so sorting needed.."; sorted=false; break;;
	  * ) echo "> Please answer \"y\" or \"n\".";;
      esac
  done

## If sort argument is 'sort' then set $sorted to false
elif [ "$sort" == 'sort' ]
then 
    echo "> Bamfile will be sorted.."
    sorted=false
fi

echo "Sorted: $sorted"
echo '++++++++++++++++++'

#echo $sort
echo ">> The current bam will be processed:"
echo $bam
# get sample name
sample=$(basename $bam .bam)
echo "Sample: $sample"

if [ $sorted = true ]
then 
    bamsort=$bam
else
    bamsort=${sample}.sorted.bam
    echo ">>Sorting bam file.."
    samtools sort -o $bamsort $bam
    echo "> Bamfile sorted.."
fi

## index bam
if [ -f $bamsort.bai ] 
then
    echo "> Index file for $bam exists"
else 
    echo ">> Creating index file for $bam.."
    samtools index -b $bamsort $bamsort.bai
fi 

# make bigwig
if [ -f $sample.bw ]
then
    echo "> Bigwig file for $sample exist.. Nothing done"
else
    echo ">> Generating bigwig file from $bam.."
    bamCoverage -b $bamsort -o $sample.bw
fi

echo ">> Indexing and bigwig-generation completed for $bam..!"
echo " ++++++++++++++++++++++++++++++++++++++++++++++++++++ "






