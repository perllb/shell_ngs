# shellRNAseq
Common scripts for NGS data: from fastq to bam


Shell commands:

## Extract sequence from fastq file
`awk '(NR%4==2)' file.fq > seq.txt`
