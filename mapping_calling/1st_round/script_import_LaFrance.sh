#!/bin/bash
#SBATCH --mail-type=BEGIN,END,FAIL
mkdir rawdata_LaFrance
cd rawdata_LaFrance 
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR617/001/SRR6175521/SRR6175521_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR617/001/SRR6175521/SRR6175521_2.fastq.gz

