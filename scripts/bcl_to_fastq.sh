#!/bin/bash
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -A lsens2017-3-2
#SBATCH -p dell
#SBATCH -t 03:00:00
#SBATCH -J bcl2fastq
#SBATCH -o bcl2fastq.out
#SBATCH -e bcl2fastq.err

## Run script in folder output from illumina sequencer.
## SampleSheet.csv has to be in this folder!

path=$1 #e.g. '/projects/fs1/medpvb/backup/projects/NextSeq500/Seq004/171215_NB502004_0010_AHN7KJAFXX'
cd $path

ml  GCCcore/6.3.0  bcl2fastq/2.19.1

bcl2fastq
