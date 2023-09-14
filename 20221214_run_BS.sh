#!/bin/bash

. /opt/conda/etc/profile.d/conda.sh
conda activate diverse_yeast_env

BASE_DIR=/home/heineikeb/alphafold/
BASE_BS=${BASE_DIR}selection_calculations/branch_site/


#for OG_BASE in OG4552_REF_Scer_AF-P37299-F1-model_v2 OG4555_REF_Scer_AF-P32799-F1-model_v2 
for OG_BASE in OG1299_REF_Scer_AF-P00549-F1-model_v2 OG1390_REF_Scer_AF-P47164-F1-model_v2 
#OG1254_REF_Scer_AF-Q01519-F1-model_v2 OG1275_REF_Scer_AF-P47052-F1-model_v2 OG1287_REF_Scer_AF-P37298-F1-model_v2 OG2006_REF_Scer_AF-P07257-F1-model_v2 OG2248_REF_Scer_AF-P07143-F1-model_v2 OG2704_REF_Scer_AF-P21801-F1-model_v2 OG2714_REF_Scer_AF-P00427-F1-model_v2 OG3208_REF_Scer_AF-P00128-F1-model_v2 OG3505_REF_Scer_AF-P04037-F1-model_v2 OG3575_REF_Scer_AF-P08067-F1-model_v2 OG4118_REF_Scer_AF-P07255-F1-model_v2 OG4352_REF_Scer_AF-P00127-F1-model_v2 OG4751_REF_Scer_AF-P04039-F1-model_v2 OG5490_REF_Scer_AF-P32340-F1-model_v2
#OG4555_REF_Scer_AF-P32799-F1-model_v2 OG4316_REF_Scer_AF-P00424-F1-model_v2 OG4150_REF_Scer_AF-P07256-F1-model_v2 OG4346_REF_Scer_AF-P10174-F1-model_v2 OG4360_REF_Scer_AF-P08525-F1-model_v2 OG2112_REF_Scer_AF-P33421-F1-model_v2 OG4744_REF_Scer_AF-P22289-F1-model_v2 OG4552_REF_Scer_AF-P37299-F1-model_v2 
#OG1122_REF_Scer_AF-P13711-F1-model_v2 
#OG4552_REF_Scer_AF-P37299-F1-model_v2 
#OG4150_REF_Scer_AF-P07256-F1-model_v2 
#OG4150_REF_Scer_AF-P07256-F1-model_v2 
#OG2603_REF_Scer_AF-P50076-F1-model_v2 OG2845_REF_Scer_AF-P43577-F1-model_v2 OG3677_REF_Scer_AF-P47125-F1-model_v2 OG1299_REF_Scer_AF-P00549-F1-model_v2

#for ALN_FILE in ${BASE_DIR}msas/structural/tm_align/cds_trim_strict/OG*.tm.fasta.clipkit.cds
do 
#    DIR_OG_BASE=$(echo $ALN_FILE | cut -d '.' -f 1)
#    OG_BASE=$(echo $DIR_OG_BASE | cut -d '/' -f 9)
    echo $OG_BASE
    ##Make directory
    CALC_DIR=${BASE_BS}${OG_BASE}/
        
    mkdir $CALC_DIR
    
    #Change into directory and call codeML
    cd $CALC_DIR

    ## Copy in tree
    #cp ${BASE_DIR}msas/structural/tm_align/trees/${OG_BASE}.tm.fasta.clipkit.treefile.renamed ${CALC_DIR}tree.treefile
    cp ${BASE_BS}trees/${OG_BASE}.bs_trees ${CALC_DIR}tree.treefile

    ## Copy in phy
    cp ${BASE_DIR}msas/structural/tm_align/cds_trim_strict/${OG_BASE}.tm.fasta.clipkit.cds.renamed.codeML.phy ${CALC_DIR}aln.phy
    
    ## Copy in control files
    cp ${BASE_BS}codemlBS_A.ctl ${CALC_DIR}codemlBS_A.ctl
    cp ${BASE_BS}codemlBS_Afix1.ctl ${CALC_DIR}codemlBS_Afix1.ctl
    
    ## Run codeML Branch Site model
    /var/lib/paml/paml-4.10.6/bin/codeml ${CALC_DIR}codemlBS_A.ctl
    
    mkdir ${CALC_DIR}A
    
    for OUTPUT_FILE in 2NG.dN 2NG.dS 2NG.t 4fold.nuc BSA.out lnf rst rst1 rub
    do
       mv ${CALC_DIR}${OUTPUT_FILE} ${CALC_DIR}A/${OUTPUT_FILE} 
    done

    # Run codeML Branch Site Fixed Omega model (background)
    /var/lib/paml/paml-4.10.6/bin/codeml ${CALC_DIR}codemlBS_Afix1.ctl
   
    mkdir ${CALC_DIR}Afix1
    
    for OUTPUT_FILE in 2NG.dN 2NG.dS 2NG.t 4fold.nuc BSAfix1.out lnf rst rst1 rub
    do
       mv ${CALC_DIR}${OUTPUT_FILE} ${CALC_DIR}Afix1/${OUTPUT_FILE} 
    done

done
