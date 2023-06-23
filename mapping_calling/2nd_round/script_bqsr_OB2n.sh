#!/bin/bash
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=50G
#SBATCH -c 2

#SBATCH -p unlimitq
module load bioinfo/gatk-4.2.2.0
module load bioinfo/samtools-1.14
module load bioinfo/picard-2.20.7

ind=$(echo "OB2n")

gatk --java-options "-Xmx20G" \
BaseRecalibrator \
	-I /work/project/irhs/Rose/roses_dedupRG_bam/"$ind".dedupRG.bam \
	-R /home/tleroy/work/Rose/fastaref/OBDH_1.0_formated.fasta \
	--known-sites /work/project/irhs/Rose/rose32.jointvcf.clean.SNPSelectVariants.vcf.gz \
	-O /work/project/irhs/Rose/rose32.jointvcf.clean.SNPSelectVariants.vcf.gz."$ind"recal_data.table

gatk --java-options "-Xmx20G" \
ApplyBQSR \
   -R /home/tleroy/work/Rose/fastaref/OBDH_1.0_formated.fasta \
   -I /work/project/irhs/Rose/roses_dedupRG_bam/"$ind".dedupRG.bam \
   --bqsr-recal-file /work/project/irhs/Rose/rose32.jointvcf.clean.SNPSelectVariants.vcf.gz."$ind"recal_data.table \
   -O /work/project/irhs/Rose/roses_BQSRed_bam/"$ind".dedupRG.BQSRed.bam
