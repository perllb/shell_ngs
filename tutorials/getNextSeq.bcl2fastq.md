1. Get NextSeq500 run output folder from bone

```shell
cd <path to where you want to put the run output folder>
ml ncftp # load module for transfer
ncftpget -R -u <username> bone.bmc.lu.se . /nextseq/<run folder (e.g. Seq001)>  
```

1. Get SampleSheet from NGS folder in dropbox

```shell
cd <path to the NGS folder and the Seq0XX folder that should contain the SampleSheet for the run>
sftp -P 22220 <username>@aurora-ls2.lunarc.lu.se
 > put <samplesheet.csv>
```
     1. Change name of the sample sheet to SampleSheet.csv

 ```shell
 mv <old sample sheet name> SampleSheet.csv
 ```

1. Make script to send bcl2fastq to computing node
### Bcl2fastq must be run from the project folder (e.g. Seq004/171215_NB502004_0010_AHN7KJAFXX ), and this folder must contain the SampleSheet.csv!
```shell
#!/bin/bash

#SBATCH -n 8
#SBATCH -N 1
#SBATCH -A lsens2017-3-2
#SBATCH -p dell
#SBATCH -t 03:00:00
#SBATCH -J bcl2fastq
#SBATCH -o bcl2fastq.out
#SBATCH -e bcl2fastq.err

cd <NextSeq500 Run folder> # eg cd /projects/fs1/medpvb/backup/projects/NextSeq500/Seq004/171215_NB502004_0010_AHN7KJAFXX

 # load modules
ml  GCCcore/6.3.0  bcl2fastq/2.19.1

 # run bcl2fastq
bcl2fastq

```
