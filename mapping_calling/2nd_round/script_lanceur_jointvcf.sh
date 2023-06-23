#!/bin/bash
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=80G
#SBATCH -c 10
#SBATCH -p unlimitq
module load bioinfo/gatk-4.2.2.0
module load bioinfo/samtools-1.14
module load bioinfo/picard-2.20.7

cd /work/project/irhs/Rose/roses_BQSRed_vcf/
# gatk index ref & gvcf
#gatk --java-options "-Xmx4G" CreateSequenceDictionary -R ~/work/Rose/fastaref/OBDH_1.0_formated.fasta
#samtools faidx ~/work/Rose/fastaref/OBDH_1.0_formated.fasta

#for i in *g.vcf.gz; do
#	gatk --java-options "-Xmx4G" IndexFeatureFile -I $i
#done

gatk --java-options "-Xmx80G" GenotypeGVCFs \
   -R /home/tleroy/work/Rose/fastaref/OBDH_1.0_formated.fasta \
   -V multigvcf_round2.g.vcf.gz \
   -O rose_round2.jointvcf.gz \
   --all-sites true
