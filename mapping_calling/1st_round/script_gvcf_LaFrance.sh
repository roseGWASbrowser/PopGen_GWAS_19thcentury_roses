#!/bin/bash
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=12G
#SBATCH -c 8
#SBATCH -p unlimitq
module load bioinfo/gatk-4.2.2.0
module load bioinfo/samtools-1.14
module load bioinfo/picard-2.20.7
ind=$(echo "LaFrance")

# add read group to bam
java -Xmx4g -jar $PICARD AddOrReplaceReadGroups \
I="$ind".dedup.bam \
O="$ind".dedupRG.bam  \
RGLB=lib1 \
RGPL=ILLUMINA \
RGPU=barcode \
RGSM="$ind" \

rm "$ind".dedup.bam
# index
samtools index "$ind".dedupRG.bam 

# gatk index ref & gvcf
#gatk --java-options "-Xmx4G" CreateSequenceDictionary -R ~/work/Rose/fastaref/OBDH_1.0_formated.fasta
#samtools faidx ~/work/Rose/fastaref/OBDH_1.0_formated.fasta
gatk --java-options "-Xmx4G" HaplotypeCaller \
--sample-name $ind \
-R ~/work/Rose/fastaref/OBDH_1.0_formated.fasta \
-I "$ind".dedupRG.bam \
--native-pair-hmm-threads 8 \
-ERC GVCF \
-O "$ind".g.vcf
