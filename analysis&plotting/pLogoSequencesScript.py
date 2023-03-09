#!/usr/bin/env python3

# pLogoSequencesScript
#
# Generate the kmers necessary for the pLogo plotting
#
# Bogdan Fedeles, Nov 2022.




import MutLib as ml


ref_gpt = "EG10_rgc_Corrected.fasta"
ref_twnstr = "mm10_twst-probes50.fa"

#mutpos_file = "Dev/79_FINAL2.mutpos"
#mut_file = "Dev/8217_mgmt.mut"


#mutpos files
mutfilespath = "../MutposUnique/"
name_suff1 = "_FINAL2_unique.mutpos"
name_suff2 = "_FINAL_unique.mutpos"

out_path = "../pLogoFiles/"
GA_suff = "GtoA_15base_"
AG_suff = "AtoG_15base_"

filerange1 = range(74, 84)
filerange2 = range(349, 359)

GAmut = "G>A"
AGmut = "A>G"

interval = 7


#extract_mutpos_contexts(mutpos_file, ref_gpt, "sample79_19base_G-A.txt", "G>A", 9)


for f in filerange1:

    mut_file = mutfilespath+str(f)+name_suff1
    GA_file = out_path+GA_suff+str(f)+'.txt'
    AG_file = out_path+AG_suff+str(f)+'.txt'

    print('Proceesing file '+mut_file)
    
    ml.extract_mutpos_contexts(mut_file, ref_gpt, GA_file, GAmut, interval) 
    ml.extract_mutpos_contexts(mut_file, ref_gpt, AG_file, AGmut, interval)


for f in filerange2:

    mut_file = mutfilespath+str(f)+name_suff2
    GA_file = out_path+GA_suff+str(f)+'.txt'
    AG_file = out_path+AG_suff+str(f)+'.txt'

    print('Proceesing file '+mut_file)
    
    ml.extract_mutpos_contexts(mut_file, ref_gpt, GA_file, GAmut, interval) 
    ml.extract_mutpos_contexts(mut_file, ref_gpt, AG_file, AGmut, interval)

