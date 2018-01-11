#!/bin/sh
ml foss/2016b
ml Python/3.5.2
ml GCC/5.4.0-2.26  OpenMPI/1.10.3
ml SAMtools/1.4
ml ucsc-tools/3.4.3
ml BEDTools/2.26.0


usage="$(basename "$0") [-h] [-b <bamfile>] [-s 'y'/'n'] [-bed] -- program to convert bam files to bw and/or bed

where:
    -h  show this help text
    -b <bamfile>  name of bamfile to convert
    -s <y/n>  'y' if sorting needed, 'n' if already sorted!
    -bed  if set, bedfile will be generated

    "

seed=42
while getopts ':hb:s:bed:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    b)
    if [ -r $bam ]
      then
        echo "> Bamfile $OPTARG exists and is readable!.. proceeding..";
        bam=$OPTARG
    elif [ ! -r $bam ]
      then
        echo "> Bamfile $OPTARG does not exists or is not readable!"
        exit 1
    fi
       ;;
    s) if [ $OPTARG == 'y' || $OPTARG == 'n' ]
      then sort=$OPTARG
    else
      echo "> ERROR: Sort argument invalid!"
      echo "$usage"
    fi
       ;;
    bed)
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

echo "Sort?: $sort"
echo '++++++++++++++++++'

echo ">> The current bam will be processed:"
echo $bam
# get sample name
sample=$(basename $bam .bam)
echo "Sample: $sample"

if [ $sort = 'n' ]
then
    bamsort=$bam
else
    bamsort=${sample}.sorted.bam
    echo ">> Sorting bam file.."
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

## Bed conversion?
if [ $bed == 'y']
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
