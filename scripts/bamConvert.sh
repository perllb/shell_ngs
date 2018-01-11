#!/bin/sh

echo "> Loading modules..."
ml foss/2016b
ml Python/3.5.2
ml GCC/5.4.0-2.26  OpenMPI/1.10.3
ml SAMtools/1.4
ml ucsc-tools/3.4.3
ml BEDTools/2.26.0


usage="Usage:

> $(basename "$0") [-h] [-i <bamfile>] [-s 'y'/'n'] [-b] -- program to convert bam files to bw and/or bed

where:
    -h  show this help text
    -i <bamfile>  name of bamfile to convert
    -s <y/n>  'y' if sorting needed, 'n' if already sorted!
    -b  if set, bedfile will be generated

    "

seed=42
while getopts ':hi:s:b' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    i)
    if [ -r $bam ]
      then
        bam=$OPTARG
    elif [ ! -r $bam ]
      then
        echo "> ERROR: Bamfile $OPTARG does not exists or is not readable!"
        exit 1
    fi
       ;;
    s) if [ $OPTARG == 'y' ] || [ $OPTARG == 'n' ]
	  then sort=$OPTARG
       else
      echo "> ERROR: Sort argument invalid!"
      echo "$usage"
    fi
       ;;
    b)
      bed='y'
      ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done
shift $((OPTIND - 1))

if [ -z "$bam" ]
  then
    printf "> Error: bamfile not entered..\n\n"
    echo "$usage" >&2
    exit 1
fi
if [ -z "$sort" ]
  then
    printf "> Error: Sort not specified..\n\n"
    echo "$usage" >&2
    exit 1
fi

echo ">> The current bam will be processed:"
printf "$bam\n"

# get sample name
sample=$(basename $bam .bam)
bamsort=${sample}.sorted.bam
if [ $sort = 'y' ]
then
    echo ">> Sorting bam file.."
    samtools sort -o $bamsort $bam
    echo "> Bamfile sorted.."
# if sorted not desired, then check if the bamsort file with .sorted.bam suffix exists.
# If not, it might imply that the original bamfile is sorted and should then be used.
else
    if [ ! -f $bamsort ]
    then
	     bamsort=$bam
    fi
fi


## index bam
if [ -f $bamsort.bai ]
then
    echo "> Index file for $bam exists"
else
    echo ">> Creating index file for $bamsort.."
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

## Bed conversion?
if [ $bed == 'y' ]
  then
    # make bed from bam
    if [ -f $sample.bed ]
    then
        echo "> Bed file for $sample exist.. Nothing done"
    else
        echo "> Generating bed file from $bam.."
        bedtools bamtobed -i $bam > $sample.bed
        echo "> Bedfile generated: $sample.bed"
    fi
fi

echo ">> Indexing and bigwig-generation completed for $bam..!"
echo " ++++++++++++++++++++++++++++++++++++++++++++++++++++ "
