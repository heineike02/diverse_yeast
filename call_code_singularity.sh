#!/bin/bash

cd ~/singularity
#export SINGULARITY_BIND='/mnt/g/My Drive/Crick_LMS/projects/diverse_yeasts/alphafold':/home/heineike_wsl2/alphafold,/mnt/c/Users/heineib/Documents/Github:/home/heineike_wsl2/github_s

#sequence based alignments
#singularity exec codeml.sif /home/heineikeb/github/diverse_yeast/20230110_seq_align.sh
#singularity exec codeml.sif /home/heineikeb/github/diverse_yeast/20230324_seq_align_dnds_trees.sh

#structure based alignments
#singularity exec codeml.sif /home/heineikeb/github/diverse_yeast/20221206_struct_align_dnds_trees.sh
#singularity exec codeml.sif /home/heineikeb/github/diverse_yeast/20221206_struct_align_dnds_msas.sh

#Trees for mitochondrial encoded proteins
#singularity exec codeml.sif /home/heineikeb/github/diverse_yeast/20241023_struct_align_dnds_msas_mito.sh
#singularity exec codeml.sif /home/heineikeb/github/diverse_yeast/20241023_struct_align_dnds_trees_mito.sh

#Tree calculations (for comparing alignments)
#singularity exec codeml.sif /home/heineikeb/github/diverse_yeast/20230326_calculate_tree_quantities.sh

#DN/DS
#singularity exec codeml.sif /home/heineikeb/github/diverse_yeast/20221207_run_yn00.sh
#singularity exec codeml.sif /home/heineikeb/github/diverse_yeast/20221214_run_m0.sh
#singularity exec codeml.sif /home/heineikeb/github/diverse_yeast/20230309_run_m0_troubleshoot.sh
#singularity exec codeml.sif /home/heineikeb/github/diverse_yeast/20231010_run_site_model.sh
#singularity exec codeml.sif /home/heineikeb/github/diverse_yeast/20230926_call_hyphy_busted.sh
#singularity exec codeml.sif /home/heineikeb/github/diverse_yeast/20230119_run_m1.sh
#singularity exec codeml.sif /home/heineikeb/github/diverse_yeast/20221214_run_BS.sh

#Feature Breakdown
#singularity exec codeml.sif /home/heineikeb/github/diverse_yeast/20230913_struct_align_dnds_msas_binding_site.sh
#singularity exec codeml.sif /home/heineikeb/github/diverse_yeast/20240109_struct_align_dnds_msas_uniprot_binding_site.sh
#singularity exec codeml.sif /home/heineikeb/github/diverse_yeast/20240916_struct_align_dnds_msas_core_surface.sh
#singularity exec codeml.sif /home/heineikeb/github/diverse_yeast/20240109_run_m0_features.sh

