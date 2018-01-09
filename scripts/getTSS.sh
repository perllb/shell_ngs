# get TSS from gencode.v27 annotation:
# 1. Get transcript entries 
# 2. first (if +) or last (-) coordinate should be TSS

cd ~/genomicData/gencode/hg38

# get transcript rows only
grep transcript -w ~/genomicData/gencode/hg38/gencode.v27.annotation.gtf > gencode.v27.transcript.gtf

# get  exon rows only
grep exon -w ~/genomicData/gencode/hg38/gencode.v27.annotation.gtf > gencode.v27.exon.gtf

# convert exon gtf to bed
python ~/tmp/python_ngs/scripts/gtfToBed.py -i gencode.v27.exon.gtf -f transcript_name 

# get TSS - write to bed format
awk -F "\t" ' { if ($7 == "+" ) print $1,$4,$4+1,$9,".",$7; else print $1,$5-1,$5,$9,".",$7; } ' OFS='\t' gencode.v27.transcript.gtf > gencode.v27.TSS.bed

# overlap L1HS-PA4 with exons
# pipe bedtools output, and get max left and min right coordinate to get the part that overlaps
bedtools intersect -a ~/genomicData/Repeats/hg38/L1/L1HS.L1PA2-4.hg38.uniq.bed -b gencode.v27.exon.transcript_name.bed -f 1E-9 -wo | awk ' { if ( $2 > $8 ) 
{ if ($3 > $9) {print $0,$1,$2,$9} else {print $0,$1,$2,$3}} 
else 
{ if ($3 > $9) {print $0,$1,$8,$9} else {print $0,$1,$8,$3}}} ' OFS="\t" > gencode.v27.exon_L1hs.pa4.txt

# make bedfile with the part that overlap
awk ' { if ( $6 == $12 ) 
{ print $1,$15,$16,$4":"$10,".",$12; } 
else 
{ print $1,$15,$16,$4":AS:"$10,".",$12; }
} ' OFS="\t" gencode.v27.exon_L1hs.pa4.txt > gencode.v27.exon_L1hs.pa4.bed

## Overlap with TSS
bedtools intersect -a gencode.v27.exon_L1hs.pa4.bed -b gencode.v27.TSS.bed -wa > gencode.v27.TSS_L1hs.pa4.bed

## Make gtf
python ~/tmp/python_ngs/scripts/bedToGtf.py -i gencode.v27.TSS_L1hs.pa4.bed -a genc.v27.TSS_L1hs.pa4 -l TSS_fusion -o gencode.v27.TSS_L1hs.pa4.gtf

mkdir TSS/L1
mv *L1* TSS/L1/
