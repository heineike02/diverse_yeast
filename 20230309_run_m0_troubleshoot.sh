#!/bin/bash

. /opt/conda/etc/profile.d/conda.sh
conda activate diverse_yeast_env

BASE_DIR=/home/heineikeb/alphafold/
#/home/heineike_wsl2/alphafold/
BASE_M0=${BASE_DIR}selection_calculations/m0_yn00_comparison/m0/


for OG_BASE in OG1273_REF_Scer_AF-P40012-F1-model_v2 
#OG3104_REF_Scer_AF-P53954-F1-model_v2 OG1364_REF_Scer_AF-P27796-F1-model_v2 OG1230_REF_Scer_AF-P29465-F1-model_v2

#OG4*.tm.fasta.clipkit.renamed
#${BASE_DIR}msas/structural/tm_align/trees/
#for A in /home/heineikeb/alphafold/msas/structural/tm_align/trees/OG4*.tm.fasta.clipkit.treefile.renamed; 
#do 
#    echo $A
#done
#for TREE_FILE in /home/heineikeb/alphafold/msas/structural/tm_align/trees/OG*.tm.fasta.clipkit.treefile.renamed
do 
    #DIR_OG_BASE=$(echo $TREE_FILE | cut -d '.' -f 1)
    #OG_BASE=$(echo $DIR_OG_BASE | cut -d '/' -f 9)
    #echo $TREE_FILE
    #echo $DIR_OG_BASE
    echo $OG_BASE
    #Make directory
    CALC_DIR=${BASE_M0}${OG_BASE}/
    echo $CALC_DIR    
    mkdir $CALC_DIR
    
    #Change into directory and call codeML
    cd $CALC_DIR

    ## Copy in tree
    cp ${BASE_M0}/alignments/${OG_BASE}.filtered ${CALC_DIR}tree.treefile
    
    # Copy in phy
    cp ${BASE_M0}/alignments/${OG_BASE}.m0_filtered.codeML.phy ${CALC_DIR}aln.phy
    
    # Copy in control file
    cp ${BASE_M0}m0.ctl ${CALC_DIR}m0.ctl
    
    # Run codeML
    /var/lib/paml/paml-4.10.6/bin/codeml ${CALC_DIR}m0.ctl


    # remove data
    # for now skipping it
    # paml_results_out = base_dir + os.sep + os.path.normpath('selection_calculations/yn00/' + og_base + '.yn00.out')
    # shutil.copy(calc_dir + og_base + '.yn00.out', paml_results_out)


    #python /home/heineike_wsl2/github_s/diverse_yeast/20221206_struct_align_yn00.py

done
