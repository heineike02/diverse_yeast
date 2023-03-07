#!/bin/bash
cd ~/singularity-ce-3.10.4
export SINGULARITY_BIND='/mnt/g/My Drive/Crick_LMS/projects/diverse_yeasts/alphafold':/home/heineike_wsl2/alphafold,/mnt/c/Users/heineib/Documents/Github:/home/heineike_wsl2/github_s

#singularity exec codeml.sif ls /home/heineike_wsl2/github_s/diverse_yeast/20221206_struct_align_dnds_preps.sh

#singularity exec codeml.sif /home/heineike_wsl2/github_s/diverse_yeast/20230301_run_yn00_test_cleandata.sh
#singularity exec codeml.sif /home/heineike_wsl2/github_s/diverse_yeast/20221206_struct_align_dnds_preps.sh
#singularity exec codeml.sif /home/heineike_wsl2/github_s/diverse_yeast/20221207_run_yn00.sh
#singularity exec codeml.sif /home/heineike_wsl2/github_s/diverse_yeast/20221214_run_m0.sh
#singularity exec codeml.sif /home/heineike_wsl2/github_s/diverse_yeast/20221220_test_m0.sh
#singularity exec codeml.sif /home/heineike_wsl2/github_s/diverse_yeast/20230110_seq_align.sh

#Test yn00 while removing sequences
#singularity exec codeml.sif /home/heineike_wsl2/github_s/diverse_yeast/20220301_filter_sparse_seqs.sh
singularity exec codeml.sif /home/heineike_wsl2/github_s/diverse_yeast/20230301_run_yn00_filterseqs.sh