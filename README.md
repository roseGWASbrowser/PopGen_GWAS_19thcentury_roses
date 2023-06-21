# PopGen_GWAS_19thcentury_roses

This repository contains the scripts needed to redo the analyses done by Thibault Leroy et al. (preprint available on bioRxiv) to reconstruct the history of rose breeding during the 19th century, including the analysis of the SNP array dataset, the mapping and SNP calling of the whole-genome sequences, as well as all subsequent analyses and figures. Given the number of analyses performed, some scripts (and more broadly information) could be missing. Please send me an email.

Preprint: Leroy T., Albert E., Thouroude T., Baudino S., Caissard J-C., Chastellier A., Chameau J., Jeauffre J., Loubert T., Paramita S.N., Pernet A., Soufflet-Freslon A., Oghina-Pavie C., Foucher C., Hibrand-Saint Oyant L, Clotault J.
*Dark side of the honeymoon: reconstructing the Asian x European rose breeding history through the lens of genomics*, bioRxiv

Last update: 21/06/23 (under progress) 

### 1/ SNP array dataset (see ./SNParray)

### 2/ Whole-genome sequences (from the raw sequencing to the final SNP set, see ./mapping_calling)

Main softwares: bwa-mem2, picard & GATK4
For the reanalysis of publicly available raw sequncing reads, we downloaded the sequences using wget from the SRA archive (e.g. see *script_import_LaFrance.sh*). All other sequences were de novo sequenced and therefore directly downloaded from the sequencing facility. Then, we then used Trimmomatic (e.g. see *script_trimmomatic_LaFrance*). We then used fastqc and multiqc to check the quality of the sequencing data. Then, we used bwa-mem2 to map the paired-end reads, as well as the single-end reads after trimming (see *script_mapping_LaFrance.sh*, note that several samples available from SRA were sequenced with a SE strategy: the botanical roses from Hibrand Saint-Oyant and the Yellow Island cultivar). Mapping results from PE and SE reads were then merged with samtools and picard was then used to remove duplicates. <br>


### 3/ Population structure, kinship and diversity estimates (see ./popgenomics)

### 4/ GWAS analyses (see ./GWAS)

### 5/ Scripts and files to generate the main figures (see ./main_figures)

### 6/ Website (see ./website)

