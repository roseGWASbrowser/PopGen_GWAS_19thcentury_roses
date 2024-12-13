#!/bin/python
import pysam
import pandas as pd

# Open the VCF file
vcf_file = "rose32.jointvcf.clean.PASSonly.varonly.vcf.withheader"
vcf = pysam.VariantFile(vcf_file)

# Initialize a list to store the allele counts
allele_counts = ""

# Process each record in the VCF file
for record in vcf:
    # Get the position and alleles
    chrom = record.chrom
    pos = record.pos
    ref = record.ref
    alts = record.alts
    #print(chrom,pos,ref,alts)
    # Loop through each sample and extract AD if available
    if len(ref) == 1 and len(alts) == 1:
        # initiate allele_counts
        line_allele_counts=chrom+"\t"+str(pos)+"\t"+ref+"\t"+alts[0]
        sample_nb= 0
        for sample_name in record.samples:
            sample = record.samples[sample_name]
            #print("incl.",chrom,"-",pos,"/",ref,":",alts,"|",sample)
            # Extract the AD field (allelic depths for ref and alt alleles)
            if 'AD' in sample:
                sample_nb+=1
                ad_values = sample['AD']
                #print(ad_values)
                # Append reference and alternate allele counts
                if ad_values is not None:
                    # Add counts to the table as [chromosome, position, ref_count, alt_count]
                    #allele_counts.append([chrom, pos,sample_name,ref, alts[0], ad_values[0], ad_values[1]])
                    line_allele_counts=line_allele_counts+"\t"+str(ad_values[0])+"\t"+str(ad_values[1])
                else:
                    #allele_counts.append([chrom, pos, sample_name,ref, alts[0], None, None])
                    line_allele_counts=line_allele_counts+"\t"+"NA"+"\t"+"NA"
    if sample_nb==32:
        #print("incl. 32 samples at this position car ",sample_nb)
        line_allele_counts=line_allele_counts+"\n"
        allele_counts=allele_counts+line_allele_counts
    else:
        print("excluded != 32 samples at this position EXCLUDED car ",sample_nb)

#    else:
#        print("excl.",chrom,"-",pos,"/",ref,":",alts,"EXCLUDED")

print(allele_counts)

# Convert the list to a DataFrame for better handling
#df = pd.DataFrame(
#    allele_counts, 
#    columns=["Chromosome", "Position", "Ref", "Alt", "Ref_Count_sample1", "Alt_Count_sample1", "Ref_Count_sample2", "Alt_Count_sample2", "Ref_Count_sample3", "Alt_Count_sample3", "Ref_Count_sample4", "Alt_Count_sample4", "Ref_Count_sample5", "Alt_Count_sample5", "Ref_Count_sample6", "Alt_Count_sample6", "Ref_Count_sample7", "Alt_Count_sample7", "Ref_Count_sample8", "Alt_Count_sample8", "Ref_Count_sample9", "Alt_Count_sample9", "Ref_Count_sample10", "Alt_Count_sample10", "Ref_Count_sample11", "Alt_Count_sample11", "Ref_Count_sample12", "Alt_Count_sample12", "Ref_Count_sample13", "Alt_Count_sample13", "Ref_Count_sample14", "Alt_Count_sample14", "Ref_Count_sample15", "Alt_Count_sample15", "Ref_Count_sample16", "Alt_Count_sample16","Ref_Count_sample17", "Alt_Count_sample17", "Ref_Count_sample18", "Alt_Count_sample18", "Ref_Count_sample19", "Alt_Count_sample19", "Ref_Count_sample20", "Alt_Count_sample20", "Ref_Count_sample21", "Alt_Count_sample21", "Ref_Count_sample22", "Alt_Count_sample22", "Ref_Count_sample23", "Alt_Count_sample23", "Ref_Count_sample24", "Alt_Count_sample24", "Ref_Count_sample25", "Alt_Count_sample25", "Ref_Count_sample26", "Alt_Count_sample26", "Ref_Count_sample27", "Alt_Count_sample27", "Ref_Count_sample28", "Alt_Count_sample28", "Ref_Count_sample29", "Alt_Count_sample29", "Ref_Count_sample30", "Alt_Count_sample30", "Ref_Count_sample31", "Alt_Count_sample31", "Ref_Count_sample32", "Alt_Count_sample32"]
#)

# Save the DataFrame to a TSV file
#output_file = "rose32.jointvcf.clean.PASSonly.varonly.vcf.allele_counts_table.tsv"
#df.to_csv(output_file, sep="\t", index=False, na_rep="NA")

#print(f"Allele counts saved to {output_file}")
