#########################################
# Generate scripts for STAR alignments
#########################################

usage="

Usage:


> $(basename "$0") -p <project folder> -f <fastq.folder> -g <ref.genome (mm10/hg38)> -m <mapping strategy (multi/unique)> -z <zipped> -s <y/n>
where:
    -p: Main project directory (where STAR outputfolder will be placed)
    -f: Folder with all fastq files to be mapped
    -g: reference genome (either 'mm10' or 'hg38')
    -m: mapping strategy. 'multi' for multimapping (--outFilterMultimapNmax 10) and 'unique' for (--outFilterMultimapNmax 1)
    -z: set if fastq files are zipped!
    -s: if using spliceJunction database (gencode) annotation to assist mapping (y/n)
    -h help
-- Program to generate scripts for mapping to STAR

    "

go=0
seed=42
zipped=0
sjdb=0

while getopts ':hp:f:g:m:zs' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
# -p
    p)
    path=$OPTARG
    if [ ! -d $path ]; then
        echo ">> ERROR: project directory $path does not exists!"
        echo "$usage" >$2
        exit 1
    fi
     ;;

# -f
    f)
    fastq=$OPTARG

    # check if dir exist
    if [ ! -d $fastq ]; then
        echo "$usage" >&2
        echo ">> ERROR: fastq directory $fastq does not exists!"
        exit 1
    fi
    ;;


# -g
    g)
    genome=$OPTARG
    if [ $genome != 'hg38' ] && [ $genome != 'mm10' ]; then
      echo "$usage" >%2
      echo ">> ERROR: Genome -g must be set to 'hg38' or 'mm10'.."
      exit 1
    fi
    ;;

# -m
    m)
    mapping=$OPTARG
    if [ $mapping == 'multi' ]; then
      multiMapMax=10
    elif [ $mapping == 'unique' ]; then
      multiMapMax=1
    else
      echo "$usage" >%2
      echo ">> ERROR: -m must be set to either 'unique' or 'multi'"
      exit 1
    fi
    ;;

# -z
    z)
    zipped=1
    ;;
# s
    s)
    sjdb=1
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

# check if all mandatory variables set
if [ ! -z "$path" ] && [ ! -z "$fastq" ] && [ ! -z "$mapping" ] && [ ! -z "$genome" ]; then

  # if -z not set, then prompt user to make sure not zipped
  if [ $zipped == 0 ]; then
  a=0
  while [ $a == 0 ]
  do
    printf " -z not set. Are fastq files zipped? (y/n): "
    read arg
    if [ $arg == 'y' ]; then
      zipped=1
      a=1
    elif [ $arg == 'n' ]; then
      zipped=0
      a=1
    else
      echo "> Error..  Please answer 'y' (fastq are zipped *fastq.gz) or 'n' (fastq are not zipped *fastq) "
    fi
  done
  fi

  # if -s (sjdb) not set, then prompt user to make no sjdb should be used
  if [ $sjdb == 0 ]; then
  a=0
  while [ $a == 0 ]
  do
    printf " -s (sjdb) not set. Do you want to use gencode annotation for SJ assistance during alignment? (y/n): "
    read arg
    if [ $arg == 'y' ]; then
      printf "$gencode will be used as SJDB in alignment.. \n"
      sjdb=1
      a=1
    elif [ $arg == 'n' ]; then
      printf "No SJDB will be used in alignment.. \n"
      sjdb=0
      a=1
    else
      echo "> Error.. Please answer 'y' (Use gencode as sjdb) or 'n' (do not use gencode as sjdb) "
    fi
  done
  fi

  # check if there are fastq files in the fir
  if [ $zipped == '1' ]; then
    count=$(ls -1 $fastq/*.fastq.gz 2>/dev/null | wc -l)
  else
    count=$(ls -1 $fastq/*.fastq 2>/dev/null | wc -l)
  fi
  if [ $count == 0 ]; then
    echo "$usage" >&2
    echo ">> ERROR: no fastq files in fastq dir $fastq.."
    echo ">> Or they are zipped! in that case, add -z option!"
    exit 1
  fi

  ## set ref genome folder and gencode annot.
  if [ $genome == 'mm10' ]; then
      genomedir='/projects/fs1/common/genome/lunarc/indicies/star/mouse/mm10'
      gencode='/projects/fs1/medpvb/no_backup/genomicData/mm10/gencode/gencode.vM16/gencode.vM16.annotation.gtf'

  elif [ $genome == 'hg38' ]; then
       genomedir='/projects/fs1/common/genome/lunarc/indicies/star/human/hg38'
       gencode='/projects/fs1/medpvb/no_backup/genomicData/hg38/gencode/gencode.v27/gencode.v27.annotation.gtf'
  else
      echo ">ERROR: Set correct reference genome!"
      echo "$usage" >&2
      exit 1
  fi

  ## Check if gnome folder exist
  if [ ! -d $genomedir ]; then
    echo ">> ERROR: Genome dir: $genomedir does not exist! "
    echo "$usage" >&2
    exit 1
  fi
  echo "> Genome dir:   $genomedir"

  ## if sjdb -s set
  if [ $sjdb == 1 ]; then
    if [ ! -f $gencode ]; then
      echo ">> ERROR: Gencode annotation file $gencode does not exist! "
      echo "$usage" >&2
      exit 1
    else
      echo "> Gencode sjdb: $gencode"
      gencode="--sjdbGTFfile $gencode"
    fi
  else
    echo "> No SJDB used!"
  fi

  # star output dir
  starout=$path/Aligned_${genome}_STAR_$mapping
  if [ ! -d $starout ]; then
    mkdir $starout
  fi

  # script output dir
  scrout=$path/scripts/Align_${genome}.$mapping/
  if [ ! -d $scrout ]; then
    mkdir $scrout
  fi

  # get input fastq files
  cd $fastq
  # grep for read mate 1 (R1) to get samples
  files1=$(ls *R1.fastq.gz)
  samples=$(echo $files1 | sed 's/.R1.fastq.gz//g')
  echo $samples

  # generate scripts
  for sample in $samples;

  do echo "SAMPLE: $sample"

  echo "#!/bin/bash
#SBATCH -n 16
#SBATCH -N 1
#SBATCH -A lsens2017-3-2
#SBATCH -p dell
#SBATCH -t 01:00:00
#SBATCH -J %j.STARr_$sample
#SBATCH -o %j.STARr_$sample.out
#SBATCH -e %j.STARr_$sample.err
ml GCC/4.9.3-binutils-2.25 STAR/2.5.0a
STAR --genomeDir $genomedir  \
 --readFilesIn $fastq/${sample}.R1.fastq $fastq/${sample}.R2.fastq \
 --readFilesCommand gunzip -c \
 --outFilterMismatchNoverLmax 0.03 \
 --runThreadN 16 --outSAMattributes All \
 --outSAMtype BAM Unsorted \
  $gencode \
 --outFilterMultimapNmax $multiMapMax \
 --outFileNamePrefix $starout/$genome.$mapping.$sample" > $scrout/Map_$mapping_$genome.$sample.sh

  done

  echo ""
  echo ""
  echo ""
  echo ""
  echo "> Scripts are put in $scrout. Go there to sbatch them! "


# If not all variables set
else
  echo "$usage" >&2
  if [ -z $path ]
    then
      echo "ERROR: Missing -p <project folder>!"
  fi
  if [ -z $fastq ]
    then
      echo "ERROR: Missing -f <fastq-path>!"
  fi
  if [ -z $genome ]
    then
      echo "ERROR: Missing -g <reference genome>!"
  fi
  if [ -z $mapping ]
    then
      echo "ERROR: Missing -m <mapping strategy>!"
  fi
fi
