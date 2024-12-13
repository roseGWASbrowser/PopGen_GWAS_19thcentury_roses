#!/bin/bash
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=thibault.leroy@inrae.fr
#SBATCH --cpus-per-task=20
#SBATCH --mem=250G
# TL - 031124
module load statistics/R/4.1.1
module load devel/python/Python-3.7.9
cd /work/genphyse/cytogen/Thibault/beegenomics_2023disk/beegenomics/Rose
#python script_generateallelecounts.py
#gunzip rose32.jointvcf.clean.PASSonly.varonly.vcf.gz
Rscript script_poolfstat_PCA_031124_frommpythonmatrix_nvxtest.R
