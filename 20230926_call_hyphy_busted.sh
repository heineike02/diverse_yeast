#!/bin/bash

. /opt/conda/etc/profile.d/conda.sh
conda activate diverse_yeast_env

BASE_DIR=/home/heineike/alphafold/


for ALN_FILE in /home/heineike/alphafold/msas/structural/tm_align/fasta_renamed/OG*.tm.fasta
do 
	#Remove the file extensions
	DIR_OG_BASE=$(echo $ALN_FILE | cut -d '.' -f 1)
	
	#Remove the directory and leave just the filename without the extension
    #OG_BASE=OG1299_REF_Scer_AF-P00549-F1-model_v2
	OG_BASE=$(echo $DIR_OG_BASE | cut -d '/' -f 11)
    echo $ALN_FILE
    echo $DIR_OG_BASE
    echo $OG_BASE
    #Make directory
    CALC_DIR=${BASE_DIR}selection_calculations/hyphy_busted${OG_BASE}
    echo $CALC_DIR    
    mkdir $CALC_DIR

	cd $CALC_DIR

	hyphy busted --alignment ${BASE_DIR}msas/structural/tm_align/cds_trim_strict/${OG_BASE}.tm.fasta.clipkit.cds.renamed.codeML.phy --tree ${BASE_DIR}msas/structural/tm_align/trees/${OG_BASE}.tm.fasta.clipkit.treefile.renamed --output ${CALC_DIR} -m
