#!/bin/bash
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=40G
#SBATCH -c 10
#SBATCH -p unlimitq
module load bioinfo/gatk-4.2.2.0
module load bioinfo/samtools-1.14
module load bioinfo/picard-2.20.7

# gatk index ref & gvcf
#gatk --java-options "-Xmx4G" CreateSequenceDictionary -R ~/work/Rose/fastaref/OBDH_1.0_formated.fasta
#samtools faidx ~/work/Rose/fastaref/OBDH_1.0_formated.fasta

#for i in *g.vcf.gz; do
#	gatk --java-options "-Xmx4G" IndexFeatureFile -I $i
#done

cd /work/project/irhs/Rose/roses_BQSRed_vcf

gatk --java-options "-Xmx40G" CombineGVCFs \
   -R /home/tleroy/work/Rose/fastaref/OBDH_1.0_formated.fasta \
   --variant 175581.4n.BQSRed.g.vcf \
   --variant 175583.4n.BQSRed.g.vcf \
   --variant 175584.BQSRed.g.vcf \
   --variant 175585.4n.BQSRed.g.vcf \
   --variant 175586.4n.BQSRed.g.vcf --variant 175587.4n.BQSRed.g.vcf \
   --variant 175623.4n.BQSRed.g.vcf \
   --variant 175624.4n.BQSRed.g.vcf --variant 175626.4n.BQSRed.g.vcf \
   --variant 175627.4n.BQSRed.g.vcf \
   --variant 175628.4n.BQSRed.g.vcf --variant 177142.4n.BQSRed.g.vcf \
   --variant 177143.4n.BQSRed.g.vcf \
   --variant OB2n.BQSRed.g.vcf \
   --variant Rchisan.BQSRed.g.vcf \
   --variant Rchimut.BQSRed.g.vcf \
   --variant Rgaloff.4n.BQSRed.g.vcf \
   --variant Rchispo.BQSRed.g.vcf \
   --variant Rogig.BQSRed.g.vcf \
   -O multigvcf_round2.g.vcf.gz
