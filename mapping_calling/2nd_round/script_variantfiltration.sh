#!/bin/bash
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=150G
#SBATCH -c 10
#SBATCH -p unlimitq
module load bioinfo/gatk-4.2.2.0
module load bioinfo/samtools-1.14
module load bioinfo/picard-2.20.7

cd /work/project/irhs/Rose/roses_BQSRed_vcf

# gatk index ref & gvcf
#gatk --java-options "-Xmx4G" CreateSequenceDictionary -R ~/work/Rose/fastaref/OBDH_1.0_formated.fasta
#samtools faidx ~/work/Rose/fastaref/OBDH_1.0_formated.fasta

#for i in *g.vcf.gz; do
#	gatk --java-options "-Xmx4G" IndexFeatureFile -I $i
#done
gatk --java-options "-Xmx150G" \
    VariantFiltration \
    --QUIET true \
    -verbosity ERROR \
    -O rose_round2.jointvcf.clean.gz \
    -V rose_round2.jointvcf.gz \
    --filter-name "FAILED_QUAL" --filter-expression "QUAL < 0" \
    --filter-name "FAILED_SOR"  --filter-expression "SOR > 4.000"\
    --filter-name "FAILED_MQ"   --filter-expression "MQ < 30.00" \
    --filter-name "FAILED_QD"   --filter-expression "QD < 2.00" \
    --filter-name "FAILED_FS"   --filter-expression "FS > 60.000" \
    --filter-name "FAILED_MQRS" --filter-expression "MQRankSum < -20.000" \
    --filter-name "FAILED_RPR"  --filter-expression "ReadPosRankSum < -10.000 || ReadPosRankSum > 10.000"
