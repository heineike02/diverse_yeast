#!/bin/bash

. /opt/conda/etc/profile.d/conda.sh
conda activate diverse_yeast_env

BASE_DIR=/home/heineike_wsl2/alphafold/
BASE_YN00=${BASE_DIR}selection_calculations/m0_yn00_comparison/yn00/


#OG4150_REF_Scer_AF-P07256-F1-model_v2 OG2603_REF_Scer_AF-P50076-F1-model_v2 OG2845_REF_Scer_AF-P43577-F1-model_v2 OG3677_REF_Scer_AF-P47125-F1-model_v2 OG1299_REF_Scer_AF-P00549-F1-model_v2

#for ALN_FILE in ${BASE_DIR}msas/structural/tm_align/cds_trim_strict/OG*.tm.fasta.clipkit.cds
for OG_BASE in OG1230_REF_Scer_AF-P29465-F1-model_v2 OG1364_REF_Scer_AF-P27796-F1-model_v2
do 
    #DIR_OG_BASE=$(echo $ALN_FILE | cut -d '.' -f 1)
    #OG_BASE=$(echo $DIR_OG_BASE | cut -d '/' -f 9)
    echo $OG_BASE
    #Make directory
    CALC_DIR=${BASE_YN00}${OG_BASE}_filtered/
        
    mkdir $CALC_DIR
    
    #Change into directory and call codeML
    cd $CALC_DIR
    
    # Copy in alignment file phy
    cp ${BASE_YN00}alignments/${OG_BASE}.yn00_filtered.codeML.phy ${CALC_DIR}aln.phy
    
    # Copy in control file
    cp ${BASE_YN00}yn00.ctl ${CALC_DIR}yn00.ctl
    
    # Run codeML
    /var/lib/paml/paml-4.10.6/bin/yn00 ${CALC_DIR}yn00.ctl

    # remove data
    # for now skipping it
    # paml_results_out = base_dir + os.sep + os.path.normpath('selection_calculations/yn00/' + og_base + '.yn00.out')
    # shutil.copy(calc_dir + og_base + '.yn00.out', paml_results_out)
    
    #was doing this in python originally
    #python /home/heineike_wsl2/github_s/diverse_yeast/20221206_struct_align_yn00.py

done
