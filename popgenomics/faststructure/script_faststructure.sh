# WITHOUT CROSS VALIDATION
for i in {1..6}; do echo starting simulation for K=$i; python ~/Software/faststructure/structure.py -K $i --input rose32.jointvcf.clean.PASSonly.varonly.biallelic50kwithheader --output rose32.jointvcf.clean.PASSonly.varonly.biallelic50kwithheader --full --format bed; done

# WITH CROSS VALIDATION BUT AN INFINITATE NUMBER OF FAILED OBSERVED WEIRD
#for i in {1..6}; do echo starting simulation for K=$i; python ~/Software/faststructure/structure.py -K $i --input rose32.jointvcf.clean.PASSonly.varonly.biallelic50kwithheader --output rose32.jointvcf.clean.PASSonly.varonly.biallelic50kwithheader --full --cv 10 --format bed; done     
