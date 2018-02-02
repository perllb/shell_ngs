#!/bin/bash
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -A lsens2017-3-2
#SBATCH -p dell
#SBATCH -t 03:00:00
#SBATCH -J %j.stringTie.pipe
#SBATCH -o %j.stringTie.pipe.out
#SBATCH -e %j.stringTie.pipe.err

## 1. make folder to output gtf and bedfiles
if [ ! -d assembly ]
then
  mkdir assembly
fi

guide="/projects/fs1/medpvb/no_backup/genomicData/hg38/gencode/gencode.v27/gencode.v27.annotation.gtf"
## 2. Run stringtie
# sample 107
echo "> Run stringtie assembly .. "
sh stringtie.run.guide.sh ../Aligned_hg38_STAR_unique/hg38.unique.P5660_107_S47Aligned.sortedByCoord.out.bam 107 genc.v27 $guide
sh stringtie.run.guide.sh ../Aligned_hg38_STAR_unique/hg38.unique.P5660_110_S50Aligned.sortedByCoord.out.bam 110  genc.v27 $guide
sh stringtie.run.guide.sh ../Aligned_hg38_STAR_unique/hg38.unique.P5660_111_S51Aligned.sortedByCoord.out.bam 111 genc.v27 $guide
sh stringtie.run.guide.sh ../Aligned_hg38_STAR_unique/hg38.unique.P5660_112_S52Aligned.sortedByCoord.out.bam 112 genc.v27 $guide
sh stringtie.run.guide.sh ../Aligned_hg38_STAR_unique/hg38.unique.P5660_108_S48Aligned.sortedByCoord.out.bam 108 genc.v27 $guide
sh stringtie.run.guide.sh ../Aligned_hg38_STAR_unique/hg38.unique.P5660_109_S49Aligned.sortedByCoord.out.bam 109 genc.v27 $guide

## 3. Merge transcripts
# First, make mergelist
ls assembly/*.genc.*gtf > assembly/mergelist.txt
merged=assembly/stringtie_merged.gencGuide.gtf
echo "> Merge transcripts"
sh stringtie_merge.sh

# get number of transcripts
#cat $merged  | grep -v "^#" | awk '$3=="transcript" {print}' | wc -l

# compare the assembled transcripts to known transcripts
echo "> Compare assembled to known.."
if [ ! -f $merged ]
then
  gffcompare -r $guide -G -o assembly/merged assembly/stringtie_merged.gencGuide.gtf
fi
cat assembly/merged.stats

## 4. Estimate abuncances
echo "> Estimate abundances of merged transcripts.."
sh stringtie.run.merged.sh ../Aligned_hg38_STAR_unique/hg38.unique.P5660_110_S50Aligned.sortedByCoord.out.bam 110 $merged
sh stringtie.run.merged.sh ../Aligned_hg38_STAR_unique/hg38.unique.P5660_107_S47Aligned.sortedByCoord.out.bam 107 $merged
sh stringtie.run.merged.sh ../Aligned_hg38_STAR_unique/hg38.unique.P5660_108_S48Aligned.sortedByCoord.out.bam 108 $merged
sh stringtie.run.merged.sh ../Aligned_hg38_STAR_unique/hg38.unique.P5660_109_S49Aligned.sortedByCoord.out.bam 109 $merged
sh stringtie.run.merged.sh ../Aligned_hg38_STAR_unique/hg38.unique.P5660_111_S51Aligned.sortedByCoord.out.bam 111 $merged
sh stringtie.run.merged.sh ../Aligned_hg38_STAR_unique/hg38.unique.P5660_112_S52Aligned.sortedByCoord.out.bam 112 $merged
