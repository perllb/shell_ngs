# shellRNAseq
Common scripts for NGS data: from fastq to bam


Shell commands:

## Extract sequences from fastq file
`awk '(NR%4==2)' file.fq > seq.txt`

## Extract subset (shuffle) of rows in file
`shuf -n 1000 file` 
