#!/bin/bash
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=12G
#SBATCH -c 8
#SBATCH -p unlimitq
module load bioinfo/gatk-4.2.2.0
module load bioinfo/samtools-1.14
module load bioinfo/picard-2.20.7
ind=$(echo "OB2n")

cd /work/project/irhs/Rose/roses_BQSRed_bam 
samtools index "$ind".dedupRG.BQSRed.bam

# gatk index ref & gvcf
#gatk --java-options ""-Xmx4G"" CreateSequenceDictionary -R ~/work/Rose/fastaref/OBDH_1.0_formated.fasta
#samtools faidx ~/work/Rose/fastaref/OBDH_1.0_formated.fasta
gatk --java-options ""-Xmx4G"" HaplotypeCaller \
--sample-name $ind \
-R ~/work/Rose/fastaref/OBDH_1.0_formated.fasta \
-I ""$ind"".dedupRG.BQSRed.bam \
--native-pair-hmm-threads 8 \
-ERC GVCF \
-O ""$ind"".BQSRed.g.vcf
