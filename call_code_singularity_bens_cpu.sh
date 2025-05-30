#!/bin/bash
cd ~/singularity-ce-4.0.1

#Scripts run on a small computer which take less computational power
#Note:  I have symbolic links outside of singularity for these two directories as alphafold_symlink and alphafold_github since it doesn't seem to be able to overwrite the symbolic link with the bound directory.  
#It won't work whether you just call the symbolic link or the actual directory itself.  
#Until I load singularity nothing goes into those directories
#echo here

#Needed to mount google drive before running this
#sudo mount -t drvfs G: /mnt/g  

export SINGULARITY_BIND='/mnt/g/My Drive/Crick_LMS/projects/diverse_yeasts/alphafold':/home/heineike/alphafold,'/mnt/c/Documents and Settings/bheineike/Documents/GitHub':/home/heineike/github


#singularity exec codeml.sif /home/heineike_wsl2/github_s/diverse_yeast/20221206_struct_align_dnds_preps.sh
#singularity exec codeml.sif /home/heineike_wsl2/github_s/diverse_yeast/20221206_struct_align_dnds_msa_bens_cpu.sh

#DN/DS
#singularity exec codeml.sif /home/heineike/github/diverse_yeast/20230914_run_m0_subset.sh
#singularity exec codeml.sif /home/heineike_wsl2/github_s/diverse_yeast/20221220_test_m0.sh
#singularity exec codeml.sif /home/heineike_wsl2/github_s/diverse_yeast/20230301_run_yn00_test_cleandata.sh

#Test yn00 and m0 while removing sequences
#singularity exec codeml.sif /home/heineike_wsl2/github_s/diverse_yeast/20230301_struct_align_dnds_msas_yn00.sh
#singularity exec codeml.sif /home/heineike_wsl2/github_s/diverse_yeast/20230307_struct_align_dnds_msas_m0.sh
#singularity exec codeml.sif /home/heineike_wsl2/github_s/diverse_yeast/20230301_run_yn00_filterseqs.sh


#Check M0 for fully trimmed alignments for example proteins
#singularity exec codeml.sif /home/heineike/github/diverse_yeast/20231110_struct_align_dnds_msas_full_trim.sh
#singularity exec codeml.sif /home/heineike/github/diverse_yeast/20231110_run_m0_complete_trim.sh

#Calculate tree for example orthogroup
#singularity exec codeml.sif /home/heineike/github/diverse_yeast/20240215_full_og_tree.sh

#Run evcoupling
#singularity exec codeml.sif /home/heineike/github/diverse_yeast/20250122_evcoupling.sh