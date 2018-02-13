## Run generate.STAR.align.sh to generate scripts to send for STAR alignment

Usage:

```shell
> sh generate.STAR.align.sh -p <project folder> -f <fastq.folder> -g <ref.genome (mm10/hg38)> -m <mapping strategy (multi/unique)> -z <zipped> -s <y/n>
where:
    -p: Main project directory (where STAR outputfolder will be placed)
    -f: Path to folder with all fastq files to be mapped
    -g: reference genome (either 'mm10' or 'hg38')
    -m: mapping strategy. (either 'multi' or 'unique'): 'multi' for multimapping (--outFilterMultimapNmax 10) and 'unique' for (--outFilterMultimapNmax 1)
    -z: set if fastq files are zipped!
    -s: if using spliceJunction database (gencode) annotation to assist mapping (y/n)
    -h help
-- Program to generate scripts for mapping to STAR

```

## Examples:

### This script will generate STAR alignment scripts that 
  - Aligns to mm10 reference genome (-g mm10)
  - Allow for multimapping (-m multi)
  - Take zipped fastq files (-z)
  - Use gencode annotation for SJDB assistance upon mapping (-s)
  
```shell
> sh generate.STAR.align.sh -p /projects/fs1/medpvb/backup/projects/monocytesinAD/ \
-f /projects/fs1/medpvb/backup/projects/monocytesinAD/data/fastq/ \
-g mm10 \
-m multi \
-z -s  
```

### This script will generate STAR alignment scripts that 
  - Aligns to hg38 reference genome (-g mm10)
  - Keep only uniquely aligned reads (-m unique)
  - Does NOT take zipped fastq files ( NO -z argument)
  - Does NOT use gencode annotation for SJDB assistance upon mapping ( NO -s argument)
  
```shell
> sh generate.STAR.align.sh -p /projects/fs1/medpvb/backup/projects/DNMT1/ \
-f /projects/fs1/medpvb/backup/projects/DNMT1/data/fastq/ \
-g hg38 \
-m unique \

```
