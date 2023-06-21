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

gatk --java-options "-Xmx40G" CombineGVCFs \
   -R /home/tleroy/work/Rose/fastaref/OBDH_1.0_formated.fasta \
   --variant 175580.g.vcf \
   --variant 175581.g.vcf \
   --variant 175583.g.vcf \
   --variant 175584.g.vcf \
   --variant 175585.g.vcf \
   --variant 175586.g.vcf --variant 175587.g.vcf \
   --variant 175588.g.vcf \
   --variant 175623.g.vcf \
   --variant 175624.g.vcf --variant 175626.g.vcf \
   --variant 175627.g.vcf \
   --variant 175628.g.vcf --variant 177142.g.vcf \
   --variant 177143.g.vcf \
   --variant OB2n.g.vcf --variant LaFrance.g.vcf \
   --variant HumeBlush.g.vcf \
   --variant Rchisan.g.vcf \
   --variant Rchimut.g.vcf \
   --variant Rgaloff.g.vcf \
   --variant Rxdama.g.vcf \
   --variant Rmosch.g.vcf \
   --variant Rchispo.g.vcf \
   --variant Rxanth.g.vcf \
   --variant YellowIsland.g.vcf \
   --variant Rpend.g.vcf \
   --variant Rwich.g.vcf \
   --variant Rmaja.g.vcf \
   --variant Rarve.g.vcf \
   --variant Rogig.g.vcf \
   --variant Rminut.g.vcf \
   -O multigvcf.g.vcf.gz
