#!/bin/bash

. /opt/conda/etc/profile.d/conda.sh
conda activate diverse_yeast_env

BASE_DIR=/home/heineike_wsl2/alphafold/
BASE_BS=${BASE_DIR}selection_calculations/branch_site/


for OG_BASE in OG4150_REF_Scer_AF-P07256-F1-model_v2 
#OG2603_REF_Scer_AF-P50076-F1-model_v2 OG2845_REF_Scer_AF-P43577-F1-model_v2 OG3677_REF_Scer_AF-P47125-F1-model_v2 OG1299_REF_Scer_AF-P00549-F1-model_v2

#for ALN_FILE in ${BASE_DIR}msas/structural/tm_align/cds_trim_strict/OG*.tm.fasta.clipkit.cds
do 
#    DIR_OG_BASE=$(echo $ALN_FILE | cut -d '.' -f 1)
#    OG_BASE=$(echo $DIR_OG_BASE | cut -d '/' -f 9)
    echo $OG_BASE
    ##Make directory
    CALC_DIR=${BASE_BS}${OG_BASE}/
        
    #mkdir $CALC_DIR
    
    #Change into directory and call codeML
    cd $CALC_DIR

    ## Copy in tree
    #cp ${BASE_DIR}msas/structural/tm_align/trees/${OG_BASE}.tm.fasta.clipkit.treefile.renamed ${CALC_DIR}tree.treefile
    
    ## Copy in phy
    #cp ${BASE_DIR}msas/structural/tm_align/cds_trim_strict/${OG_BASE}.tm.fasta.clipkit.cds.renamed.codeML.phy ${CALC_DIR}aln.phy
    
    ## Copy in control file
    #cp ${BASE_M0}m0.ctl ${CALC_DIR}m0.ctl
    
    ## Run codeML Branch Site model
    #/var/lib/paml/paml-4.10.6/bin/codeml ${CALC_DIR}codemlBS_A.ctl
    
    # Run codeML Branch Site Fixed Omega model (background)
    /var/lib/paml/paml-4.10.6/bin/codeml ${CALC_DIR}codemlBS_Afix1.ctl


    # remove data
    # for now skipping it
    # paml_results_out = base_dir + os.sep + os.path.normpath('selection_calculations/yn00/' + og_base + '.yn00.out')
    # shutil.copy(calc_dir + og_base + '.yn00.out', paml_results_out)


    #python /home/heineike_wsl2/github_s/diverse_yeast/20221206_struct_align_yn00.py

done
