##################################################################
# This scripts finds L1s that overlap a TSS of gencode genes
# 1. get TSS coordinates for all gencode.v27 transcripts
# 2. get all l1 that overlap exon
# 3. keep those l1-exon fragments that also is part of TSS
##################################################################
# get TSS from gencode.v27 annotation

cd ~/genomicData/gencode/hg38

## 1. Get TSS coordinates of gencode transcripts
# get transcript rows only
grep transcript -w ~/genomicData/gencode/hg38/gencode.v27.annotation.gtf > gencode.v27.transcript.gtf
# get TSS - write to bed format
awk -F "\t" ' { if ($7 == "+" ) print $1,$4,$4+1,$9,".",$7; else print $1,$5-1,$5,$9,".",$7; } ' OFS='\t' gencode.v27.transcript.gtf > gencode.v27.TSS.bed

## 2. overlap L1HS-PA4 with exons
##  - make bedfile of L1s that overlap a gencode exon - and keep the part that overlap.
# get  exon rows only
grep exon -w ~/genomicData/gencode/hg38/gencode.v27.annotation.gtf > gencode.v27.exon.gtf
# convert exon gtf to bed
python ~/tmp/python_ngs/scripts/gtfToBed.py -i gencode.v27.exon.gtf -f transcript_name
# make bedfile with the coordinates of L1-exon overlap
# pipe bedtools output, and get max left and min right coordinate to get the part that overlaps
bedtools intersect -a ~/genomicData/Repeats/hg38/L1/L1.hg38.uniq.bed -b gencode.v27.exon.transcript_name.bed -f 1E-9 -wo | awk ' { if ( $2 > $8 )
{ if ($3 > $9) {print $0,$1,$2,$9} else {print $0,$1,$2,$3}}
else
{ if ($3 > $9) {print $0,$1,$8,$9} else {print $0,$1,$8,$3}}} ' OFS="\t" > gencode.v27.exon_L1.txt
# make bedfile with the part that overlap (extend coordinates with 5 bases on each end to avoid annotation bias)
awk ' { if ( $6 == $12 )
{ print $1,$15-5,$16+5,$4":"$10,".",$12; }
else
{ print $1,$15-5,$16+5,$4":AS:"$10,".",$12; }
} ' OFS="\t" gencode.v27.exon_L1.txt > gencode.v27.exon_L1.bed

## 3. get those l1-exon overlaps that overlap with TSS
bedtools intersect -a gencode.v27.exon_L1.bed -b gencode.v27.TSS.bed -wa > gencode.v27.TSS_L1.bed
## Make gtf
python ~/tmp/python_ngs/scripts/bedToGtf.py -i gencode.v27.TSS_L1.bed -a genc.v27.TSS_L1 -l TSS_fusion -o gencode.v27.TSS_L1.gtf
