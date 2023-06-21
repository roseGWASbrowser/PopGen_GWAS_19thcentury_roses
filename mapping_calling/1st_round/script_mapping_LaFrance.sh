#!/bin/bash
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=30G
module load bioinfo/bwa-mem2-2.2
module load bioinfo/samtools-1.14
module load bioinfo/picard-2.20.7
ind=$(echo "LaFrance")
cd /home/tleroy/work/Rose/
# sorting
#bwa-mem2 index ./fastaref/OBDH_1.0_formated.fasta
# mapping
bwa-mem2 mem ./fastaref/OBDH_1.0_formated.fasta /home/tleroy/work/Rose/rawdata_LaFrance/SRR6175521_1_cleaned.fastq.gz /home/tleroy/work/Rose/rawdata_LaFrance/SRR6175521_2_cleaned.fastq.gz | samtools view -bS - > "$ind"_paired.bam
bwa-mem2 mem ./fastaref/OBDH_1.0_formated.fasta /home/tleroy/work/Rose/rawdata_LaFrance/SRR6175521_1_cleaned_unpaired.fastq.gz | samtools view -bS - > "$ind"_unpaired1.bam
bwa-mem2 mem ./fastaref/OBDH_1.0_formated.fasta /home/tleroy/work/Rose/rawdata_LaFrance/SRR6175521_2_cleaned_unpaired.fastq.gz | samtools view -bS - > "$ind"_unpaired2.bam
# merging bams
samtools merge "$ind".bam "$ind"_paired.bam "$ind"_unpaired1.bam "$ind"_unpaired2.bam
# samtools sort 
samtools sort "$ind".bam > "$ind".sorted.bam
rm "$ind".bam
# samtools index 
samtools index "$ind".sorted.bam
# picard remove duplicates
java -Xmx4g -jar $PICARD MarkDuplicates  \
        I="$ind".sorted.bam \
        O="$ind".dedup.bam \
        M="$ind".MarkDuplicates.metrics.txt \
        REMOVE_DUPLICATES=true
# first round = no BQSR => need to generate a first round of variant calling
rm "$ind".sorted.bam
