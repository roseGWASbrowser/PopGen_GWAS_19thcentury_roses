This repository includes the scripts that were used to generate the input file for the PCA based on allele counts proposed by Mathieu Gautier (our PCI Evol Biol recommender).<br>

I first extracted allele counts from the AD column of the VCF with the script ("script_generateallelecounts.py"). <br>

Then, the idea was to generate the largest possible dataset to perform the PCA based on the allele read counts. We successfully applied the analysis with a success on 65 million SNPs of rose SNPs on the genotoul cluster. Note that this step required >200Gb of RAM. 
<code>175580  175581  175583 ... (32 individuals in total)</code><br>
<code>6,0     13,0    1,0 ... </code><br>
<code>6,0     12,0    1,0 ... </code><br>
<code>2,0     0,4     2,0 ... </code><br>
<code>0,0     0,4     1,0 ... </code><br><br>

Then the R script (PCA_poolfstat_131224.R) was used to plot the results and the generate the PDF (SupFigureX_PCA_allelecounts_vs_PCAcallsFig2B_131224.pdf).<br>
