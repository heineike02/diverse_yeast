#!/bin/bash

. /opt/conda/etc/profile.d/conda.sh
conda activate diverse_yeast_env

BASE_DIR=/home/heineike/alphafold/
OG_BASE=OG1299_REF_Scer_AF-P00549-F1-model_v2

cd ${BASE_DIR}selection_calculations/hyphy_busted

mkdir ${OG_BASE}

hyphy busted --alignment ${BASE_DIR}msas/structural/tm_align/cds_trim_strict/${OG_BASE}.tm.fasta.clipkit.cds.renamed.codeML.phy --tree ${BASE_DIR}msas/structural/tm_align/trees/${OG_BASE}.tm.fasta.clipkit.treefile.renamed --output ${BASE_DIR}/selection_calculations/hyphy_busted/${OG_BASE} -m
