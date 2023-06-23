#!/bin/bash
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=50G
#SBATCH -c 2

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

#cd /work/project/irhs/Rose

#gzip rose32.jointvcf.clean.PASSonly.varonly.vcf 

gatk --java-options "-Xmx20G" \
IndexFeatureFile \
     -I /work/project/irhs/Rose/rose32.jointvcf.clean.gz

gatk --java-options "-Xmx20G" \
    SelectVariants \
    -V /work/project/irhs/Rose/rose32.jointvcf.clean.gz \
    --select-type-to-include SNP \
    -O /work/project/irhs/Rose/rose32.jointvcf.clean.SNPSelectVariants.vcf.gz

gatk --java-options "-Xmx50G" \
    VariantsToTable \
    -V /work/project/irhs/Rose/rose32.jointvcf.clean.SNPSelectVariants.vcf.gz \
    -F CHROM -F POS -F ID  -F REF -F ALT -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
    -O /work/project/irhs/Rose/rose32.jointvcf.clean.SNPSelectVariants.vcf.gz.table 
