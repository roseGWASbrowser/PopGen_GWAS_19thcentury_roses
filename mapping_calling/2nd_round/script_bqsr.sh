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
BaseRecalibrator \
	-I /work/project/irhs/Rose/roses_dedupRG_bam/175626.dedupRG.bam \
	-R /home/tleroy/work/Rose/fastaref/OBDH_1.0_formated.fasta \
	--known-sites /work/project/irhs/Rose/rose32.jointvcf.clean.SNPSelectVariants.vcf.gz \
	-O /work/project/irhs/Rose/rose32.jointvcf.clean.SNPSelectVariants.vcf.gz.175626recal_data.table

gatk --java-options "-Xmx20G" \
ApplyBQSR \
   -R /home/tleroy/work/Rose/fastaref/OBDH_1.0_formated.fasta \
   -I /work/project/irhs/Rose/roses_dedupRG_bam/175626.dedupRG.bam \
   --bqsr-recal-file /work/project/irhs/Rose/rose32.jointvcf.clean.SNPSelectVariants.vcf.gz.175626recal_data.table \
   -O /work/project/irhs/Rose/roses_BQSRed_bam/175626.dedupRG.BQSRed.bam
