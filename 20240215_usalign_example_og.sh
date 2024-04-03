#!/bin/bash


PDB_DIR=/home/heineike/alphafold_ln/selected_proteins/pdbs/OG1022/
PDB_LIST=${PDB_DIR}pdb_list.txt
/home/heineike/usalign/USalign -dir $PDB_DIR $PDB_LIST -suffix .pdb -mm 4 > ${PDB_DIR}us_align.fasta
