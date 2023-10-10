#!/bin/bash
cd ~/singularity
#export SINGULARITY_BIND='/mnt/g/My Drive/Crick_LMS/projects/diverse_yeasts/alphafold':/home/heineike_wsl2/alphafold,/mnt/c/Users/heineib/Documents/Github:/home/heineike_wsl2/github_s

#singularity exec codeml.sif /home/heineike_wsl2/github_s/diverse_yeast/20221207_run_yn00.sh
#singularity exec codeml.sif /home/heineikeb/github/diverse_yeast/20221206_struct_align_dnds_trees.sh
#singularity exec codeml.sif /home/heineikeb/github/diverse_yeast/20221206_struct_align_dnds_msas.sh
#singularity exec codeml.sif /home/heineikeb/github/diverse_yeast/20221214_run_m0.sh
#singularity exec codeml.sif /home/heineikeb/github/diverse_yeast/20230119_run_m1.sh
singularity exec codeml.sif /home/heineikeb/github/diverse_yeast/20221214_run_BS.sh
#singularity exec codeml.sif /home/heineikeb/github/diverse_yeast/20230309_run_m0_troubleshoot.sh
#singularity exec codeml.sif /home/heineikeb/github/diverse_yeast/20230110_seq_align.sh
#singularity exec codeml.sif /home/heineikeb/github/diverse_yeast/20230324_seq_align_dnds_trees.sh
#singularity exec codeml.sif /home/heineikeb/github/diverse_yeast/20230326_calculate_tree_quantities.sh
