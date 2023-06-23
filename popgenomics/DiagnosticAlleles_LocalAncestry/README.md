To identify the diagnostic alleles, we first computed the allele frequency at all PASS SNPs with the python script *script_freqgroup4localancestryprop.py*.<br>
Then, we filtered SNPs with fixed different alleles between the ancient Asian and the ancient European. <br>
The file containing the diagnostic SNPs is provided: see *rose_round2.jointvcf.clean.AllFreqRefPerGroup.diagnostics.txt*<br><br>
Then, we computed local ancestry by considering average or median allele frequecy across sliding windows spanning the genome, see the R script for details *script_R_localancestry_300322.R*
