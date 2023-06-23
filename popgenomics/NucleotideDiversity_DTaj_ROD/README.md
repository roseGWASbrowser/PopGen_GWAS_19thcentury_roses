Regarding the reconstruction of the fasta sequences from the vcf, this section has been extensively described in Leroy et al. 2021 Current Biology and Leroy et al. 2021 Peer Community Journal. <br>
Here, the main change is associated to a reconstruction taking into account the ploidy level (note that it works for a up to 4 ploidy level, see the python script *VCF2Fasta_fast_withcovquantiles_polyploid_outputfilenames.py*). <br><br>
Importantly, the reconstructed sequences are **unphased**. This is appropriate to generate nucleotide diversity estimates for instance, but these sequences **should not be used for methods that requires phased data** (e.g. MSMC. )!<br><br>
The rest of the files can be used to explore the results. It includes the values of pi, Tajima's D and the reduction of diversity (RoD index). To have an overview of the results, see script_R_diversity_pi_ROD_280322.R<br>

