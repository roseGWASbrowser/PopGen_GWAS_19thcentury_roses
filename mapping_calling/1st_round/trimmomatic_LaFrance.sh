#!/bin/bash
#SBATCH --mail-type=BEGIN,END,FAIL
cd /home/tleroy/work/Rose/rawdata_LaFrance
module load bioinfo/Trimmomatic-0.38
#wget -c --no-check-certificate --auth-no-challenge --user 'thibault.leroy' --password 'CHILD|LAUGHED|equal|SUCCESS|19' https://ngs.vbcf.ac.at/filemanager/byurl/lSJSNPaost-HKKK2DRXY_1_R12489_20211118.tar.gz
#tar -xvzf lSJSNPaost-HKKK2DRXY_1_R12489_20211118.tar.gz
java -jar /usr/local/bioinfo/src/Trimmomatic/Trimmomatic-0.38/trimmomatic.jar PE -threads 1 -phred33 SRR6175521_1.fastq.gz SRR6175521_2.fastq.gz SRR6175521_1_cleaned.fastq.gz SRR6175521_1_cleaned_unpaired.fastq.gz SRR6175521_2_cleaned.fastq.gz SRR6175521_2_cleaned_unpaired.fastq.gz ILLUMINACLIP:/usr/local/bioinfo/src/Trimmomatic/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
