#!/bin/bash
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=12G
#SBATCH -c 1
#SBATCH -p workq
module load bioinfo/gatk-4.2.2.0
module load bioinfo/samtools-1.14
module load bioinfo/picard-2.20.7
ind=$(echo "LaFrance")

# index
samtools depth -aa -H -o "$ind".dedupRG.bam.cov "$ind".dedupRG.bam 
awk 'NR % 100 == 0' "$ind".dedupRG.bam.cov  > "$ind".dedupRG.bam.cov.1pc
gzip "$ind".dedupRG.bam.cov.1pc
rm "$ind".dedupRG.bam.cov
