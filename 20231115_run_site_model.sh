#!/bin/bash

. /opt/conda/etc/profile.d/conda.sh
conda activate diverse_yeast_env

BASE_DIR=/home/heineikeb/alphafold/
#/home/heineike_wsl2/alphafold/
BASE_SEL_CALC=${BASE_DIR}selection_calculations/site_model/


#OG_BASE=OG1004_REF_Scer_AF-P40459-F1-model_v2
for OG_BASE in OG3575_REF_Scer_AF-P08067-F1-model_v2 OG4352_REF_Scer_AF-P00127-F1-model_v2 OG2714_REF_Scer_AF-P00427-F1-model_v2 OG4346_REF_Scer_AF-P10174-F1-model_v2 OG4751_REF_Scer_AF-P04039-F1-model_v2 OG1122_REF_Scer_AF-P13711-F1-model_v2

#Examples: 

#Complete
#OG4316_REF_Scer_AF-P00424-F1-model_v2 

#Failed to complete: 
#OG2248_REF_Scer_AF-P07143-F1-model_v2


#for OG_BASE in OG4150_REF_Scer_AF-P07256-F1-model_v2 OG2603_REF_Scer_AF-P50076-F1-model_v2 OG2845_REF_Scer_AF-P43577-F1-model_v2 OG3677_REF_Scer_AF-P47125-F1-model_v2 OG1299_REF_Scer_AF-P00549-F1-model_v2 
#Also ran 2006 QCR 2 before 

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
    #echo $OG_BASE
    #Make directory
    CALC_DIR=${BASE_SEL_CALC}${OG_BASE}/
    echo $CALC_DIR    
    mkdir $CALC_DIR
    
    #Change into directory and call codeML
    cd $CALC_DIR

    ## Copy in tree
    cp ${BASE_DIR}msas/structural/tm_align/trees/${OG_BASE}.tm.fasta.clipkit.treefile.renamed ${CALC_DIR}tree.treefile
    
    # Copy in phy
    cp ${BASE_DIR}msas/structural/tm_align/cds_trim_strict/${OG_BASE}.tm.fasta.clipkit.cds.renamed.codeML.phy ${CALC_DIR}aln.phy
    
    # Copy in control file
    cp ${BASE_SEL_CALC}site_model.ctl ${CALC_DIR}site_model.ctl
    #cp ${BASE_M0}m0.ctl ${CALC_DIR}m0.ctl

    # Run codeML
    /var/lib/paml/paml-4.10.6/bin/codeml ${CALC_DIR}site_model.ctl
    #/var/lib/paml/paml-4.10.6/bin/codeml ${CALC_DIR}m0.ctl

    # remove data
    # for now skipping it
    # paml_results_out = base_dir + os.sep + os.path.normpath('selection_calculations/yn00/' + og_base + '.yn00.out')
    # shutil.copy(calc_dir + og_base + '.yn00.out', paml_results_out)


    #python /home/heineike_wsl2/github_s/diverse_yeast/20221206_struct_align_yn00.py

done
