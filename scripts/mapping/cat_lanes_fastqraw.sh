#!/usr/bin
###########################################################################
# Script to fuse raw fastq of different lanes for each sample / read mate
###########################################################################

cd <main folder (with each sample having a folder in which the raw fastq files are)>
#E.g cd /projects/fs1/medpvb/backup/projects/hEmbryo-dissection/PE_150bp/fastq/
## /fastq (<- Main folder!)
##      /1.sample
##        /1samp.L001.R1.fastq.gz
##        /1samp.L002.R1.fastq.gz
##        /1samp.L003.R1.fastq.gz
##        /1samp.L001.R2.fastq.gz
##        /1samp.L002.R2.fastq.gz
##        /1samp.L003.R2.fastq.gz
##      /2.sample
##        /2samp.L001.R1.fastq.gz
##        /2samp.L002.R1.fastq.gz
##        /2samp.L003.R1.fastq.gz
##        /2samp.L001.R2.fastq.gz
##        /2samp.L002.R2.fastq.gz
##        /2samp.L003.R2.fastq.gz
##

## Get all sample names, by listing folders in current directory
samples=$(ls -d 1*/)
echo $samples

## For each sample, cat files of for each read mate in pair
for sample in $samples;
do echo $sample;

    # cd into sample's dir
    samp=$(echo $sample | sed 's/\///g')
    cd $samp
    currdir=$(pwd)
    echo "now in $currdir.. \n"

    # remove potentially catted files already there.. (careful)
    rm $samp.R1.fastq.gz
    rm $samp.R2.fastq.gz

    # create empty files to store the catted files
    touch $samp.R1.fastq.gz
    touch $samp.R2.fastq.gz

    # cat all files of each mate
    cat *R1_001.fastq.gz > $samp.R1.fastq.gz
    cat *R2_001.fastq.gz > $samp.R2.fastq.gz

    # move the ready files back to the main folder
    mv $samp.R1.fastq.gz ..
    mv $samp.R2.fastq.gz ..

    # move back to main folder
    cd ..

    currdir=$(pwd)
    echo "now in $currdir.. \n"

done
