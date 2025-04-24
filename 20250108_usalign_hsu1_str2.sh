#!/bin/bash

#Needed to mount google drive before running this
#sudo mount -t drvfs G: /mnt/g  

PDB_BASE=/home/heineike/alphafold_ln/examples/hsu1_str2/

PDB_DIR=${PDB_BASE}/pdbs/
PDB_LIST=${PDB_BASE}/pdb_list_hsu1.txt
/home/heineike/usalign/USalign -dir $PDB_DIR $PDB_LIST -suffix .pdb -mm 4 > ${PDB_BASE}/us_align_hsu1_str2.fasta
