# PopGenomics & GWAS of the 19th century roses <br>
<br>
This repository contains the scripts needed to redo the analyses done by Thibault Leroy et al. (preprint available on bioRxiv) to reconstruct the history of rose breeding during the 19th century, including the analysis of the SNP array dataset, the mapping and SNP calling of the whole-genome sequences, as well as all subsequent analyses and figures. Given the number of analyses performed, some scripts (and more broadly information) could be missing. Please send me an email.<br>

```
Preprint: Leroy T., Albert E., Thouroude T., Baudino S., Caissard J-C., Chastellier A., Chameau J., Jeauffre J., Loubert T., Paramita S.N., Pernet A., Soufflet-Freslon A., Oghina-Pavie C., Foucher C., Hibrand-Saint Oyant L, Clotault J.
Dark side of the honeymoon: reconstructing the Asian x European rose breeding history through the lens of genomics, bioRxiv 2023.06.22.546162 (last update (v2): April 07, 2024); doi: https://doi.org/10.1101/2023.06.22.546162
```
<i>Dedication: This work is dedicated to the memory of our colleague, Laurence Hibrand-Saint Oyant, who passed away in November 2024 (see the dedication in the Supplementary Information of our latest bioRxiv version).</i><br>
<br>
Author: Thibault Leroy (thibault.leroy_at_inrae.fr)<br>
Last update: 13/12/24 <br><br>

### 1/ SNP array dataset (see ./SNParray)

Main software: [fitPoly](https://cran.r-project.org/web/packages/fitPoly/index.html) <br>
In this section, we have made the raw data available (CEL files corresponding to three 96-well plates), as well as the final dataset after the filtering of the SNPs and the clone-correction (204 genotypes). The R codes to generate this matrix are also provided. In bried, we filtered the low accurate SNP based on the raw signal prior to the use of fitpoly (see *script_QC_filtering_renaming_before_fitPoly_270_samples.R*), before to launch fitpoly (*script_launch_genotyping_fitPoly_270samplesCARO.R*). Note that fitpoly was run on a specific computing cluster given the large computational requirements (here: 30 CPUs used). We then filtered polymorphic markers and investigated the ploidy (e.g. *script_diagnostic_ploidy_plot_270_perchromosome.R*). Note that the final dataset for the SNP array is available (you have to download and decompress the archive: *Dataset_SNParray_matrix_fitPoly_270_scores_formated_only_passing_probes.204CloneCorrectedsamples.txt.gz*).

### 2/ Whole-genome sequences (from the raw sequencing to the final SNP set, see ./mapping_calling)

Main softwares: [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2), [picard](https://broadinstitute.github.io/picard/) & [GATK4](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4)<br> <br>
**A/ First round of calling (see ./mapping_calling/1st_round/):**<br>
Below, we provide the scripts for one of the rose cultivar (La France), but the same strategy was used for all individuals. The only notable exception is for individuals that come from the literature and that were only sequenced with single-ends reads. For these individuals the pipeline was edited accordingly (trimming considering SE reads, mapping of SE reads), but the rest of the pipeline is similar. We also removed PCR duplicates for SE, which can be considered as a quite conservative strategy. <br>
For the reanalysis of publicly available raw sequncing reads, we downloaded the sequences using wget from the SRA archive (e.g. see *script_import_LaFrance.sh*). All other sequences were de novo sequenced and therefore directly downloaded from the sequencing facility. Then, we then used Trimmomatic (e.g. see *script_trimmomatic_LaFrance*). We then used fastqc and multiqc to check the quality of the sequencing data. Then, we used bwa-mem2 to map the paired-end reads, as well as the single-end reads after trimming (see *script_mapping_LaFrance.sh*, note that several samples available from SRA were sequenced with a SE strategy: the botanical roses from Hibrand Saint-Oyant and the Yellow Island cultivar). Mapping results from PE and SE reads were then merged with samtools and picard was then used to remove duplicates. We then used GATK4 to perform the variant calling (see the following scripts; 1/ generate a gvcf: script_gvcf_LaFrance.sh*; 2/ generate a multivcf:*script_multigvcf.sh*; 3/ perform the jointcalling:*script_lanceur_jointvcf.sh*; and 4/ filter the variants: *script_variantfiltration.sh*).<br>

**B/ Second round of calling (see ./mapping_calling/2nd_round/):** <br>
Thanks to this first round of analysis, we were able to generate an extended catalog of SNPs for roses, in order to perform a BQSR, as well as infer the ploidy level thanks to the data. For this analysis, we only considered 19 non-botanical varieties, also excluding roses with a first-degree relationship (see the analysis of kinship in section 3/). We first generated a SNP catalog based on the first SNP set (see *script_tablebeforebqsr.sh*), before to perform the BQSR for each individual (here an example based on 'Old Blush', see *script_bqsr_OB2n.sh*), before to redo the variant calling (see the following scripts; 1/ generate a gvcf: script_gvcf_OB2n.sh*; 2/ generate a multivcf:*script_multigvcf_round2.sh*; 3/ perform the jointcalling:*script_lanceur_jointvcf.sh*; and 4/ filter the variants: *script_variantfiltration.sh*). <br>

### 3/ Population structure, kinship and diversity estimates (see ./popgenomics)

Main softwares: [KING](https://www.kingrelatedness.com/), [fastStructure](https://github.com/rajanil/fastStructure)  <br> 
This section covers different analyses perform on both datasets, including the inference of the population structure based on the SNP array and a subset of the WGS data*, kinship analysis and local ancestry estimates based on the WGS data and ploidy level inferences. More detail is given within each directory (a dedicated README is available for each directory).<br>
This section also includes the PCA generated with poolFstat based on allele counts (no ploidy assumption, see ./)<br>

### 4/ GWAS analyses (see ./GWAS)

Main software: [GWASpoly](https://github.com/jendelman/GWASpoly), [circlize](https://jokergoo.github.io/circlize_book/book/) <br>
The scripts needed for the GWAS analysis, as well as the R scripts to generate the circlize plots (except Fig.4, see below), are available in this subsection. We used the R script *1_script_PCA_PCOA_on_markers_204cc_tetra_TL111021.R* to first filter the SNP set and keep only informative SNPs. Then, we used GWASpoly to perform the GWAS analysis *2_script_test_GWAS_poly_204cc_CARO_121021.R*. This latter script requires the genetic and phenotypic data, which corresponds to the TableS1 file and the Table S2 files. Then, we considered a series of R scripts to generate the circlize plots that are shown on the webste, see the section "6/ Website" below. <br><br>
Importantly, this directory contains also all the results of our GWAS, including linear Manhattan plots, qqplots and raw files used to generate the circlize plots shown on the website (p-values, q-values, etc): see ./GWAS/ALL_GWAS_RESULTS/ . Concretely, the organization of result directories mirrors the structure of the website making it easy to locate the zip files. These zip files contain three types of data (Manhattan plots, qqplots and raw files), organized into three separate folders. For any questions regarding the results, please send me an email. <br>

### 5/ Scripts and files to generate the main figures (see ./main_figures)

Main software: [circlize](https://jokergoo.github.io/circlize_book/book/) <br>
The R scripts - as well as the corresponding files containing the results - needed to generate the main figures of the article are made available on this subsection. <br>

### 6/ Website (see ./website)

This github account hosts the [rose GWAS browser](https://rosegwasbrowser.github.io/)!<br>
The current implementation of the website is based on a series of HTML webpages and a CSS style sheet. The objective of the website was to develop a user-friendly interface to explore the catalog, as a result of the elevated number of Manhattan plots generated (currently available on the website: 84 phenotypes in total * 8 models = 672 circlize plots !). Thanks to different drop-down menus, the user can select several options and the corresponding PDF is then (expected to be) loaded. All the webpages are generated by a single bash script (*script_generate_GWASbrowserwebpages.sh*) made available in the ./website directory (see also: [the github repository hosting the roseGWASbrowser website](https://github.com/roseGWASbrowser/rosegwasbrowser.github.io)).<br>
