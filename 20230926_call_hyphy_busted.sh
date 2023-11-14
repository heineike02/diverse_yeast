#!/bin/bash

. /opt/conda/etc/profile.d/conda.sh
conda activate diverse_yeast_env

BASE_DIR=/home/heineikeb/alphafold/
BASE_HYPHY=${BASE_DIR}selection_calculations/hyphy_busted/


for TREE_FILE in /home/heineikeb/alphafold/msas/structural/tm_align/trees/OG*.tm.fasta.clipkit.treefile.renamed
do 
    DIR_OG_BASE=$(echo $TREE_FILE | cut -d '.' -f 1)
    OG_BASE=$(echo $DIR_OG_BASE | cut -d '/' -f 9)
    echo $TREE_FILE
    echo $DIR_OG_BASE
    echo $OG_BASE
    #Make directory
    CALC_DIR=${BASE_HYPHY}${OG_BASE}/
    echo $CALC_DIR    
    mkdir $CALC_DIR
    
    #Change into directory and call hyphy
    cd $CALC_DIR

	hyphy busted --alignment ${BASE_DIR}msas/structural/tm_align/cds_trim_strict/${OG_BASE}.tm.fasta.clipkit.cds.renamed.codeML.phy --tree ${BASE_DIR}msas/structural/tm_align/trees/${OG_BASE}.tm.fasta.clipkit.treefile.renamed --output ${CALC_DIR}hyphy_busted.json -m
done
