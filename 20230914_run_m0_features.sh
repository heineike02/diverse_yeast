#!/bin/bash

. /opt/conda/etc/profile.d/conda.sh
conda activate diverse_yeast_env

BASE_DIR=/home/heineikeb/alphafold/
FEATURE=binding_site
#/home/heineikeb/alphafold/
#Before running, need to have made the feature directory and also added a control file into it called m0.ctl
BASE_M0=${BASE_DIR}selection_calculations/m0_features/${FEATURE}/

for OG_BASE in OG4352_REF_Scer_AF-P00127-F1-model_v2
#OG1316_REF_Scer_AF-P19097-F1-model_v2
#OG4150_REF_Scer_AF-P07256-F1-model_v2 OG2603_REF_Scer_AF-P50076-F1-model_v2 OG2845_REF_Scer_AF-P43577-F1-model_v2 OG3677_REF_Scer_AF-P47125-F1-model_v2 OG1299_REF_Scer_AF-P00549-F1-model_v2

#OG4*.tm.fasta.clipkit.renamed
#${BASE_DIR}msas/structural/tm_align/trees/
#for A in /home/heineikeb/alphafold/msas/structural/tm_align/trees/OG4*.tm.fasta.clipkit.treefile.renamed; 
#do 
#    echo $A
#done
#for TREE_FILE in /home/heineikeb/alphafold/msas/structural/tm_align/trees/OG*.tm.fasta.clipkit.treefile.renamed
#for ALN_FILE in ${BASE_DIR}msas/structural/tm_align/feature_subsets/${FEATURE}/fasta_renamed/OG*.tm.fasta
do 
    #Use these first four lines when iterating through the tree files

	#Remove the file extensions
	#DIR_OG_BASE=$(echo $ALN_FILE | cut -d '.' -f 1)
	
	#Remove the directory and leave just the filename without the extension
    #OG_BASE=$(echo $DIR_OG_BASE | cut -d '/' -f 11)
    #echo $ALN_FILE
    #echo $DIR_OG_BASE
    echo $OG_BASE
    #Make directory
    CALC_DIR=${BASE_M0}${OG_BASE}/
    echo $CALC_DIR    
    mkdir $CALC_DIR
    
    #Change into directory and call codeML
    cd $CALC_DIR

    ## Copy in tree
    cp ${BASE_DIR}msas/structural/tm_align/trees/${OG_BASE}.tm.fasta.clipkit.treefile.renamed ${CALC_DIR}tree.treefile
    
    # Copy in phy
    cp ${BASE_DIR}msas/structural/tm_align/feature_subsets/${FEATURE}/cds_trim_strict/${OG_BASE}.tm.fasta.clipkit.cds.renamed.codeML.phy ${CALC_DIR}aln.phy
    
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
